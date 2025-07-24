#!/usr/bin/env nextflow

/*
 * main.nf: Main Nextflow pipeline for metagenomic analysis (Nanopore/Illumina)
 * Modular, profile-driven, with externalized parameters and help flag
 */

nextflow.enable.dsl=2

// Print help and exit if --help is passed
if (params.help) {
    println file('help.txt').text
    exit 0
}

// Import modules
include { DORADO; PORECHOP; NANOFILT; FASTP } from './modules/preprocess.nf' 
include { KRAKENUNIQ; KAIJU; CENTRIFUGE; DIAMOND; MEGAN_LR; MERGE_TAXONOMY } from './modules/read_based.nf'
include { METAFLYE; RAVEN; UNICYCLER; MINIMAP2_ALIGN as ASSEMBLY_MAP; RACON; MEDAKA; PILON; HOMOPOLISH } from './modules/assembly.nf'
include { METAQUAST; MINIMAP2_INDEX; MINIMAP2_ALIGN as BINNING_MAP; SEMIBIN2; METABINNER; COMEBIN; METAWRAP_REFINE; AGGREGATE_BINS } from './modules/binning.nf'
include { GTDBTK; ABRICATE; PROKKA; GENOMAD; TAXPASTA; KEGG_ANNOTATION; AGGREGATE_ANNOTATIONS } from './modules/annotation.nf'
include { CHECKM2; COVERM; MERQURY; REAPR; AGGREGATE_METRICS } from './modules/evaluation.nf'

// Function to generate final report
def generateFinalReport(taxonomyReport, binningReport, annotationReport, evaluationReport) {
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    # Read reports
    taxonomy = pd.read_csv("${taxonomyReport}", sep='\t')
    binning = pd.read_csv("${binningReport}", sep='\t')
    annotation = pd.read_csv("${annotationReport}", sep='\t')
    evaluation = pd.read_csv("${evaluationReport}", sep='\t')
    
    # Merge information
    final_report = pd.merge(
        binning,
        evaluation,
        on='Bin',
        how='outer'
    )
    
    # Add taxonomic and functional information
    final_report = pd.merge(
        final_report,
        annotation[['Bin', 'Taxonomy', 'Virulence', 'Resistance']],
        on='Bin',
        how='left'
    )
    
    # Generate summary visualizations
    plt.figure(figsize=(12, 8))
    sns.scatterplot(data=final_report, x='Completeness', y='Contamination', size='Size', hue='Quality_Score')
    plt.savefig('bin_quality_plot.png')
    
    # Save reports
    final_report.to_csv('comprehensive_report.tsv', sep='\t', index=False)
    
    # Generate summary statistics
    with open('summary_stats.txt', 'w') as f:
        f.write("=== Metagenome Analysis Summary ===\\n")
        f.write(f"Total Bins: {len(final_report)}\\n")
        f.write(f"High Quality Bins (>90% complete, <5% contamination): {len(final_report[(final_report['Completeness'] > 90) & (final_report['Contamination'] < 5)])}\\n")
        f.write(f"Medium Quality Bins (>70% complete, <10% contamination): {len(final_report[(final_report['Completeness'] > 70) & (final_report['Contamination'] < 10)])}\\n")
        f.write("\\nTaxonomic Distribution:\\n")
        f.write(taxonomy.groupby('Rank')['Count'].sum().to_string())
        f.write("\\n\\nVirulence Factors Found: {len(annotation[annotation['Virulence'].notna()])}\\n")
    """
}

workflow {
    // Input validation
    if (!params.input_ont && !params.input_fast5) {
        error "Please provide either ONT FASTQ (--input_ont) or FAST5 (--input_fast5) files"
    }

    // Preprocessing
    if (params.input_fast5) {
        raw_reads = DORADO(params.input_fast5)
    } else {
        raw_reads = Channel.fromPath(params.input_ont)
    }

    trimmed_reads = PORECHOP(raw_reads)
    filtered_reads = NANOFILT(trimmed_reads)

    if (params.hybrid) {
        short_reads = FASTP(
            Channel.fromPath(params.short_1),
            Channel.fromPath(params.short_2)
        )
    }

    // Read-based taxonomy
    if (params.read_expand) {
        centrifuge_results = CENTRIFUGE(filtered_reads)
        diamond_results = DIAMOND(filtered_reads)
        megan_results = MEGAN_LR(diamond_results.daa)
        tax_results = centrifuge_results.mix(megan_results.classifications)
    } else {
        krakenuniq_results = KRAKENUNIQ(filtered_reads)
        kaiju_results = KAIJU(filtered_reads)
        tax_results = krakenuniq_results.mix(kaiju_results.report)
    }
    
    merged_taxonomy = MERGE_TAXONOMY(tax_results.collect())

    // Assembly
    if (params.hybrid) {
        assembly = UNICYCLER(filtered_reads, short_reads)
    } else if (params.raven) {
        assembly = RAVEN(filtered_reads)
    } else {
        assembly = METAFLYE(filtered_reads)
    }

    // Assembly polishing
    mapped_reads = ASSEMBLY_MAP(assembly.assembly, filtered_reads)
    racon_polished = RACON(assembly.assembly, filtered_reads, mapped_reads.bam)
    medaka_polished = MEDAKA(racon_polished.polished, filtered_reads)

    if (params.hybrid) {
        final_assembly = PILON(medaka_polished.polished, short_reads)
    } else {
        final_assembly = HOMOPOLISH(medaka_polished.polished)
    }

    // Binning
    quast_results = METAQUAST(final_assembly.polished)
    assembly_index = MINIMAP2_INDEX(final_assembly.polished)
    mapped_reads_binning = BINNING_MAP(final_assembly.polished, assembly_index.index, filtered_reads)

    if (params.poly_binning) {
        metabinner_bins = METABINNER(final_assembly.polished, mapped_reads_binning.bam)
        comebin_bins = COMEBIN(final_assembly.polished, mapped_reads_binning.bam)
        refined_bins = METAWRAP_REFINE(metabinner_bins.bins.mix(comebin_bins.bins).collect())
        final_bins = refined_bins.bins
    } else {
        semibin_results = SEMIBIN2(final_assembly.polished, mapped_reads_binning.bam)
        final_bins = semibin_results.bins
    }

    binning_summary = AGGREGATE_BINS(final_bins.collect())

    // Annotation
    gtdb_results = GTDBTK(binning_summary.bins)
    abricate_results = ABRICATE(binning_summary.bins)
    taxpasta_results = TAXPASTA(gtdb_results.summary, merged_taxonomy.merged)

    if (params.annotate_more) {
        prokka_results = PROKKA(binning_summary.bins)
        genomad_results = GENOMAD(binning_summary.bins)
        kegg_results = KEGG_ANNOTATION(prokka_results.annotation.collect())
        
        annotation_results = AGGREGATE_ANNOTATIONS(
            gtdb_results.summary
                .mix(abricate_results.reports)
                .mix(taxpasta_results.merged)
                .mix(prokka_results.annotation)
                .mix(genomad_results.results)
                .mix(kegg_results.annotation)
                .collect()
        )
    } else {
        annotation_results = AGGREGATE_ANNOTATIONS(
            gtdb_results.summary
                .mix(abricate_results.reports)
                .mix(taxpasta_results.merged)
                .collect()
        )
    }

    // Evaluation
    checkm2_results = CHECKM2(binning_summary.bins)
    coverm_results = COVERM(binning_summary.bins, mapped_reads_binning.bam)

    if (params.hybrid) {
        merqury_results = MERQURY(final_assembly.polished, short_reads)
        reapr_results = REAPR(final_assembly.polished, short_reads)
        
        evaluation_results = AGGREGATE_METRICS(
            checkm2_results.report
                .mix(coverm_results.coverage)
                .mix(merqury_results.results)
                .mix(reapr_results.results)
                .collect()
        )
    } else {
        evaluation_results = AGGREGATE_METRICS(
            checkm2_results.report
                .mix(coverm_results.coverage)
                .collect()
        )
    }

    // Generate final report
    process FINAL_REPORT {
        publishDir "${params.outdir}/final_report", mode: 'copy'

        input:
        path taxonomy from merged_taxonomy.merged
        path binning from binning_summary.summary
        path annotation from annotation_results.report
        path evaluation from evaluation_results.report

        output:
        path "comprehensive_report.tsv"
        path "summary_stats.txt"
        path "bin_quality_plot.png"

        script:
        generateFinalReport(taxonomy, binning, annotation, evaluation)
    }
}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
} 