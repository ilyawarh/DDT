#!/usr/bin/env nextflow

// Annotation module for metagenomics pipeline

process GTDBTK {
    tag { "taxonomic classification with GTDB-Tk" }
    publishDir "${params.outdir}/annotation/gtdbtk", mode: 'copy'

    input:
    path bins

    output:
    path "gtdbtk_out/*.summary.tsv", emit: summary
    path "gtdbtk_out/*", emit: all

    script:
    """
    export GTDBTK_DATA_PATH=${params.databases}/gtdbtk
    
    gtdbtk classify_wf \
        --genome_dir ${bins} \
        --out_dir gtdbtk_out \
        --extension .fa \
        --cpus ${task.cpus} \
        ${params.gtdbtk_opts}
    """
}

process ABRICATE {
    tag { "virulence and resistance genes search" }
    publishDir "${params.outdir}/annotation/abricate", mode: 'copy'

    input:
    path bins

    output:
    path "*.tsv", emit: reports

    script:
    """
    # Search against multiple databases
    for DB in card vfdb resfinder; do
        abricate --db \$DB --threads ${task.cpus} ${bins}/*.fa > \${DB}_hits.tsv
    done

    # Summarize
    abricate --summary *.tsv > summary.tsv
    """
}

process PROKKA {
    tag { "functional annotation" }
    publishDir "${params.outdir}/annotation/prokka", mode: 'copy'

    when:
    params.annotate_more

    input:
    path bin

    output:
    path "${bin.baseName}/*", emit: annotation

    script:
    """
    prokka \
        --outdir ${bin.baseName} \
        --prefix ${bin.baseName} \
        --cpus ${task.cpus} \
        ${params.prokka_opts} \
        ${bin}
    """
}

process GENOMAD {
    tag { "viral and plasmid detection" }
    publishDir "${params.outdir}/annotation/genomad", mode: 'copy'

    when:
    params.annotate_more

    input:
    path bins

    output:
    path "genomad_out/*", emit: results

    script:
    """
    genomad \
        --threads ${task.cpus} \
        all \
        ${bins} \
        genomad_out \
        ${params.databases}/genomad
    """
}

process TAXPASTA {
    tag { "taxonomy aggregation" }
    publishDir "${params.outdir}/annotation/taxpasta", mode: 'copy'

    input:
    path gtdb_summary
    path read_taxonomy

    output:
    path "taxpasta_merged.tsv", emit: merged

    script:
    """
    taxpasta merge \
        -p ${read_taxonomy} \
        -g ${gtdb_summary} \
        -o taxpasta_merged.tsv \
        ${params.taxpasta_opts}
    """
}

process KEGG_ANNOTATION {
    tag { "KEGG pathway annotation" }
    publishDir "${params.outdir}/annotation/kegg", mode: 'copy'

    when:
    params.annotate_more

    input:
    path prokka_results

    output:
    path "kegg/*", emit: annotation

    script:
    """
    mkdir -p kegg
    
    # Extract protein sequences from Prokka results
    cp ${prokka_results}/*.faa proteins.faa

    # Run DIAMOND against KEGG
    diamond blastp \
        --db ${params.databases}/kegg/genes.dmnd \
        --query proteins.faa \
        --out kegg/diamond.txt \
        --threads ${task.cpus} \
        --outfmt 6 \
        --max-target-seqs 1 \
        --evalue 1e-5

    # Map to KEGG pathways (placeholder - implement actual mapping)
    echo "KEGG pathway mapping" > kegg/pathways.txt
    """
}

process AGGREGATE_ANNOTATIONS {
    tag { "aggregating annotations" }
    publishDir "${params.outdir}/annotation", mode: 'copy'

    input:
    path '*'

    output:
    path "final_report.tsv", emit: report
    path "detailed_reports/*", emit: details

    script:
    """
    mkdir -p detailed_reports
    cp -r * detailed_reports/

    # Generate comprehensive report
    echo -e "Bin\tTaxonomy\tCompleteness\tVirulence\tResistance\tPathways" > final_report.tsv
    
    # Combine all annotations (placeholder - implement actual combining logic)
    for bin in detailed_reports/gtdbtk_out/*.summary.tsv; do
        bin_name=\$(basename \$bin)
        echo -e "\${bin_name}\tNA\tNA\tNA\tNA\tNA" >> final_report.tsv
    done
    """
}

workflow {
    // Input validation
    if (!params.databases) {
        error "Please provide the path to annotation databases (--databases)"
    }

    // Main workflow
    bins = Channel.fromPath(params.bins)
    read_taxonomy = Channel.fromPath(params.read_taxonomy)

    // Core annotation
    gtdb_results = GTDBTK(bins)
    abricate_results = ABRICATE(bins)
    
    // Merge taxonomies
    taxpasta_results = TAXPASTA(gtdb_results.summary, read_taxonomy)

    // Extended annotation if requested
    if (params.annotate_more) {
        // Process each bin individually for Prokka
        bins.flatten()
            .map { bin -> 
                def binFile = file(bin)
                return binFile
            }
            .set { individual_bins }

        prokka_results = PROKKA(individual_bins)
        genomad_results = GENOMAD(bins)
        kegg_results = KEGG_ANNOTATION(prokka_results.annotation.collect())
        
        // Aggregate all results
        all_results = gtdb_results.summary
            .mix(abricate_results.reports)
            .mix(taxpasta_results.merged)
            .mix(prokka_results.annotation)
            .mix(genomad_results.results)
            .mix(kegg_results.annotation)
            .collect()
    } else {
        // Aggregate core results only
        all_results = gtdb_results.summary
            .mix(abricate_results.reports)
            .mix(taxpasta_results.merged)
            .collect()
    }

    // Generate final report
    final_report = AGGREGATE_ANNOTATIONS(all_results)

    emit:
    report = final_report.report
    details = final_report.details
    taxonomy = taxpasta_results.merged
} 