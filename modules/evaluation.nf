#!/usr/bin/env nextflow

// Evaluation module for metagenomics pipeline

process CHECKM2 {
    tag { "genome quality assessment" }
    publishDir "${params.outdir}/evaluation/checkm2", mode: 'copy'

    input:
    path bins

    output:
    path "checkm2_results/*", emit: results
    path "checkm2_results/quality_report.tsv", emit: report

    script:
    """
    export CHECKM2DB=${params.databases}/checkm2
    
    checkm2 predict \
        --threads ${task.cpus} \
        --input ${bins} \
        --output-directory checkm2_results \
        ${params.checkm2_opts}
    """
}

process COVERM {
    tag { "genome abundance estimation" }
    publishDir "${params.outdir}/evaluation/coverm", mode: 'copy'

    input:
    path bins
    path bam

    output:
    path "*.tsv", emit: coverage

    script:
    """
    coverm genome \
        --bam-files ${bam} \
        --genome-fasta-directory ${bins} \
        --min-covered-fraction 0 \
        --methods mean trimmed_mean covered_fraction variance length count rpkm \
        -t ${task.cpus} \
        ${params.coverm_opts} \
        > coverage_stats.tsv
    """
}

process MERQURY {
    tag { "assembly quality assessment" }
    publishDir "${params.outdir}/evaluation/merqury", mode: 'copy'

    when:
    params.hybrid

    input:
    path assembly
    tuple path(short_1), path(short_2)

    output:
    path "merqury_out/*", emit: results

    script:
    """
    # Build Illumina k-mer database
    meryl k=21 count output merqury.meryl ${short_1} ${short_2}

    # Run Merqury
    merqury.sh merqury.meryl ${assembly} merqury_out
    """
}

process REAPR {
    tag { "assembly error assessment" }
    publishDir "${params.outdir}/evaluation/reapr", mode: 'copy'

    when:
    params.hybrid

    input:
    path assembly
    tuple path(short_1), path(short_2)

    output:
    path "reapr_out/*", emit: results

    script:
    """
    # Map reads
    reapr smaltmap ${assembly} ${short_1} ${short_2} mapped.bam

    # Run REAPR
    reapr pipeline ${assembly} mapped.bam reapr_out
    """
}

process AGGREGATE_METRICS {
    tag { "aggregating quality metrics" }
    publishDir "${params.outdir}/evaluation", mode: 'copy'

    input:
    path '*'

    output:
    path "final_evaluation.tsv", emit: report
    path "detailed_metrics/*", emit: details

    script:
    """
    mkdir -p detailed_metrics
    cp -r * detailed_metrics/

    # Generate comprehensive report
    echo -e "Bin\tCompleteness\tContamination\tSize\tCoverage\tQuality_Score" > final_evaluation.tsv
    
    # Parse CheckM2 results
    while IFS=\$'\\t' read -r bin completeness contamination; do
        if [[ \$bin != "Name" ]]; then
            size=\$(stat -L -c%s "detailed_metrics/\${bin}")
            coverage=\$(grep "\${bin}" detailed_metrics/coverage_stats.tsv | cut -f2)
            quality=\$(echo "scale=2; \$completeness - 5*\$contamination" | bc)
            echo -e "\${bin}\t\${completeness}\t\${contamination}\t\${size}\t\${coverage}\t\${quality}" >> final_evaluation.tsv
        fi
    done < detailed_metrics/checkm2_results/quality_report.tsv
    """
}

workflow {
    // Input validation
    if (!params.bins) {
        error "Please provide the path to genome bins (--bins)"
    }

    // Main workflow
    bins = Channel.fromPath(params.bins)
    bam = Channel.fromPath(params.bam)

    // Core evaluation
    checkm2_results = CHECKM2(bins)
    coverm_results = COVERM(bins, bam)

    // Additional evaluation for hybrid assemblies
    if (params.hybrid && params.short_1 && params.short_2) {
        assembly = Channel.fromPath(params.assembly)
        short_reads = tuple(file(params.short_1), file(params.short_2))
        
        merqury_results = MERQURY(assembly, short_reads)
        reapr_results = REAPR(assembly, short_reads)
        
        // Aggregate all results
        all_results = checkm2_results.report
            .mix(coverm_results.coverage)
            .mix(merqury_results.results)
            .mix(reapr_results.results)
            .collect()
    } else {
        // Aggregate core results only
        all_results = checkm2_results.report
            .mix(coverm_results.coverage)
            .collect()
    }

    // Generate final report
    final_report = AGGREGATE_METRICS(all_results)

    emit:
    report = final_report.report
    details = final_report.details
    quality = checkm2_results.report
    coverage = coverm_results.coverage
} 