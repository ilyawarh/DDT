#!/usr/bin/env nextflow

// Binning module for metagenomics pipeline

process METAQUAST {
    tag { "assembly quality assessment" }
    publishDir "${params.outdir}/binning/quast", mode: 'copy'

    input:
    path assembly

    output:
    path "quast_results/*", emit: results
    path "quast_results/report.tsv", emit: report

    script:
    """
    metaquast.py \
        ${assembly} \
        --threads ${task.cpus} \
        --max-ref-number 0 \
        --gene-finding \
        ${params.quast_opts}
    """
}

process MINIMAP2_INDEX {
    tag { "indexing assembly" }
    publishDir "${params.outdir}/binning/mapping", mode: 'copy'

    input:
    path assembly

    output:
    path "*.mmi", emit: index

    script:
    """
    minimap2 -d ${assembly.baseName}.mmi ${assembly}
    """
}

process MINIMAP2_ALIGN {
    tag { "mapping reads to assembly" }
    publishDir "${params.outdir}/binning/mapping", mode: 'copy'

    input:
    path assembly
    path index
    path reads

    output:
    path "*.bam", emit: bam
    path "*.bam.bai", emit: bai
    path "*.depth", emit: depth

    script:
    """
    # Map reads
    minimap2 -ax map-ont -t ${task.cpus} ${index} ${reads} | \
        samtools sort -@ ${task.cpus} -o ${assembly.baseName}.bam

    # Index BAM
    samtools index ${assembly.baseName}.bam

    # Calculate depth
    samtools depth ${assembly.baseName}.bam > ${assembly.baseName}.depth
    """
}

process SEMIBIN2 {
    tag { "binning with SemiBin2" }
    publishDir "${params.outdir}/binning/semibin2", mode: 'copy'

    when:
    !params.poly_binning

    input:
    path assembly
    path bam

    output:
    path "bins/*.fa", emit: bins
    path "*.summary", emit: summary

    script:
    """
    semibin2 single_easy_bin \
        --input-fa ${assembly} \
        --input-bam ${bam} \
        --output bins \
        --threads ${task.cpus} \
        ${params.semibin_opts}
    """
}

process METABINNER {
    tag { "binning with MetaBinner" }
    publishDir "${params.outdir}/binning/metabinner", mode: 'copy'

    when:
    params.poly_binning

    input:
    path assembly
    path bam

    output:
    path "metabinner_bins/*.fa", emit: bins

    script:
    """
    metabinner.sh \
        -a ${assembly} \
        -o metabinner_bins \
        -t ${task.cpus} \
        ${params.metabinner_opts}
    """
}

process COMEBIN {
    tag { "binning with ComeBin" }
    publishDir "${params.outdir}/binning/comebin", mode: 'copy'

    when:
    params.poly_binning

    input:
    path assembly
    path bam

    output:
    path "comebin_bins/*.fa", emit: bins

    script:
    """
    comebin.py \
        --contigs ${assembly} \
        --bam ${bam} \
        --output comebin_bins \
        --threads ${task.cpus} \
        ${params.comebin_opts}
    """
}

process METAWRAP_REFINE {
    tag { "refining bins with metaWRAP" }
    publishDir "${params.outdir}/binning/metawrap", mode: 'copy'

    when:
    params.poly_binning

    input:
    path '*'

    output:
    path "refined_bins/*.fa", emit: bins
    path "stats/*", emit: stats

    script:
    """
    metawrap bin_refinement \
        -o . \
        -A metabinner_bins \
        -B comebin_bins \
        -c 50 \
        -x 10 \
        -t ${task.cpus}
    """
}

process AGGREGATE_BINS {
    tag { "aggregating binning results" }
    publishDir "${params.outdir}/binning", mode: 'copy'

    input:
    path '*'

    output:
    path "final_bins/*.fa", emit: bins
    path "binning_summary.tsv", emit: summary

    script:
    """
    mkdir -p final_bins
    cp -L */*.fa final_bins/

    # Generate summary
    echo "Bin\tSize\tCompleteness\tContamination" > binning_summary.tsv
    for bin in final_bins/*.fa; do
        size=\$(grep -v ">" \$bin | tr -d '\\n' | wc -c)
        echo "\$(basename \$bin)\t\$size\t-\t-" >> binning_summary.tsv
    done
    """
}

workflow {
    // Input validation
    if (!params.input_ont) {
        error "Please provide ONT reads (--input_ont)"
    }

    // Main workflow
    assembly = Channel.fromPath(params.assembly)
    reads = Channel.fromPath(params.input_ont)

    // Quality assessment
    quast_results = METAQUAST(assembly)

    // Mapping
    assembly_index = MINIMAP2_INDEX(assembly)
    mapped_reads = MINIMAP2_ALIGN(assembly, assembly_index.index, reads)

    // Binning
    if (params.poly_binning) {
        // Multiple binning tools + refinement
        metabinner_bins = METABINNER(assembly, mapped_reads.bam)
        comebin_bins = COMEBIN(assembly, mapped_reads.bam)
        refined_bins = METAWRAP_REFINE(metabinner_bins.bins.mix(comebin_bins.bins).collect())
        final_bins = refined_bins.bins
    } else {
        // Single binning tool
        semibin_results = SEMIBIN2(assembly, mapped_reads.bam)
        final_bins = semibin_results.bins
    }

    // Aggregate results
    aggregated = AGGREGATE_BINS(final_bins.collect())

    emit:
    bins = aggregated.bins
    summary = aggregated.summary
    quast = quast_results.report
    bam = mapped_reads.bam
} 