#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def helpMessage() {
    println """
    ===============================
    Metagenomics QC Pipeline
    ===============================
    
    A standalone quality control workflow for metagenomic reads.
    
    Usage:
        nextflow run qc.nf [options]
    
    Options:
        --input_ont      Directory containing ONT FASTQ files
        --input_short    Directory containing paired Illumina FASTQ files
        --pattern_short  Glob pattern for paired short reads (default: *{1,2}.{fastq,fq,fastq.gz,fq.gz})
        --outdir         Output directory (default: reports/QC)
        --help          Show this message
    
    Example:
        nextflow run qc.nf --input_ont reads/ont --input_short reads/illumina
    """.stripIndent()
}

process NANOPLOT {
    tag { "ONT read quality metrics" }
    publishDir "${params.outdir}/long/nanoplot", mode: 'copy'
    maxForks 1  // Run one instance that processes all files

    input:
    path reads

    output:
    path "*"

    when:
    params.input_ont

    script:
    """
    NanoPlot \
        --fastq ${reads} \
        --outdir . \
        --threads ${task.cpus} \
        --plots hex dot
    """
}

process FASTQC_LONG {
    tag { "long read quality assessment" }
    publishDir "${params.outdir}/long", mode: 'copy'
    maxForks 1

    input:
    path reads

    output:
    path "*_fastqc.{zip,html}"

    script:
    """
    fastqc \
        -o . \
        -t ${task.cpus} \
        ${reads}
    """
}

process FASTQC_SHORT {
    tag { "short read quality assessment" }
    publishDir "${params.outdir}/short", mode: 'copy'
    maxForks 1

    input:
    path reads

    output:
    path "*_fastqc.{zip,html}"

    script:
    """
    fastqc \
        -o . \
        -t ${task.cpus} \
        ${reads}
    """
}

process MULTIQC_LONG {
    tag { "aggregate long read QC reports" }
    publishDir params.outdir, mode: 'copy'

    input:
    path '*'

    output:
    path "multiqc_long.html"
    path "multiqc_long_data"

    script:
    """
    multiqc . \
        --force \
        --filename multiqc_long \
        --title "Long Reads QC Report"
    """
}

process MULTIQC_SHORT {
    tag { "aggregate short read QC reports" }
    publishDir params.outdir, mode: 'copy'

    input:
    path '*'

    output:
    path "multiqc_short.html"
    path "multiqc_short_data"

    script:
    """
    multiqc . \
        --force \
        --filename multiqc_short \
        --title "Short Reads QC Report"
    """
}

workflow {
    // Show help message if requested
    if (params.help) {
        helpMessage()
        exit 0
    }

    // Input validation
    if (!params.input_ont && !params.input_short) {
        error "Please provide either ONT (--input_ont) or Illumina (--input_short) reads"
    }

    // Create channels for input reads
    //ont_reads = params.input_ont ? Channel.fromPath("${params.input_ont}/*.{fastq,fq,fastq.gz,fq.gz}") : Channel.empty()
    all_ont_reads = params.input_ont ? Channel.fromPath("${params.input_ont}/*.{fastq,fq,fastq.gz,fq.gz}").collect() : Channel.empty()
    
    // Create channel for paired Illumina reads
    illumina_reads = params.input_short ? Channel
        .fromFilePairs("${params.input_short}/${params.pattern_short}", checkIfExists: true)
        .map { id, reads -> reads }
        .flatten()
        .collect() : Channel.empty()

    // Long read QC
    if (params.input_ont) {
        nanoplot_results = NANOPLOT(all_ont_reads)
        fastqc_long = FASTQC_LONG(all_ont_reads)
        
        // Aggregate long read QC results
        MULTIQC_LONG(
            nanoplot_results.mix(
                fastqc_long
            ).collect()
        )
    }

    // Short read QC
    if (params.input_short) {
        fastqc_short = FASTQC_SHORT(illumina_reads)
        
        // Aggregate short read QC results
        MULTIQC_SHORT(fastqc_short.collect())
    }
}

