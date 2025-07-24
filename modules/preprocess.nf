#!/usr/bin/env nextflow

// Preprocessing module for metagenomics pipeline

process DORADO {
    tag { "basecalling_and_processing" }
    publishDir "${params.outdir}/basecalled", mode: 'copy'
    
    input:
    path fast5_dir

    output:
    path "*.fastq.gz", emit: fastq
    path "sequencing_summary.txt", optional: true, emit: summary

    when:
    params.input_fast5

    script:
    def demux_cmd = params.demultiplex ? "--kit ${params.barcode_kit} --trim_barcodes" : ""
    def mod_cmd = params.detect_modifications ? "--modified-bases ${params.modification_model}" : ""
    def duplex_cmd = params.duplex ? "--duplex" : ""
    """
    # Basecalling with optional demultiplexing and modification detection
    dorado basecaller ${params.dorado_model} \
        ${demux_cmd} \
        ${mod_cmd} \
        ${duplex_cmd} \
        ${params.dorado_opts} \
        ${fast5_dir} > basecalled.fastq

    # If demultiplexing was performed, process each barcode separately
    if [ "${params.demultiplex}" = "true" ]; then
        mkdir -p barcodes
        cat basecalled.fastq | awk -v FS=" " '{
            if (substr(\$0,1,1)==">") {
                split(\$2,a,"=")
                file="barcodes/barcode_"a[2]".fastq"
            }
            print \$0 > file
        }'
        for f in barcodes/*.fastq; do
            gzip \$f
        done
    else
        gzip basecalled.fastq
    fi
    """
}

process PORECHOP {
    tag { "adapter trimming" }
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    path reads

    output:
    path "*.trimmed.fastq.gz", emit: trimmed

    script:
    """
    porechop -i ${reads} -o ${reads.simpleName}.trimmed.fastq ${params.porechop_opts}
    gzip ${reads.simpleName}.trimmed.fastq
    """
}

process NANOFILT {
    tag { "quality filtering" }
    publishDir "${params.outdir}/filtered", mode: 'copy'

    input:
    path reads

    output:
    path "*.filtered.fastq.gz", emit: filtered

    script:
    """
    zcat ${reads} | NanoFilt ${params.nanofilt_opts} | gzip > ${reads.simpleName}.filtered.fastq.gz
    """
}

process FASTP {
    tag { "illumina preprocessing" }
    publishDir "${params.outdir}/illumina", mode: 'copy'

    input:
    tuple path(short_1), path(short_2)

    output:
    tuple path("*R1.filtered.fastq.gz"), path("*R2.filtered.fastq.gz"), emit: filtered_pe
    path "fastp.json", emit: json
    path "fastp.html", emit: html

    when:
    params.hybrid && params.input_short

    script:
    """
    fastp -i ${short_1} -I ${short_2} \
        -o ${short_1.simpleName}.filtered.fastq.gz \
        -O ${short_2.simpleName}.filtered.fastq.gz \
        ${params.fastp_opts} \
        --html fastp.html
    """
}

process FMLRC2 {
    tag { "hybrid correction" }
    publishDir "${params.outdir}/hybrid_corrected", mode: 'copy'

    input:
    path ont_reads
    tuple path(short_1), path(short_2)

    output:
    path "*.hybrid_corrected.fastq.gz", emit: corrected

    when:
    params.hybrid && params.input_short

    script:
    """
    fmlrc2 ${params.fmlrc2_opts} \
        --short_reads ${short_1} ${short_2} \
        --long_reads ${ont_reads} \
        --output ${ont_reads.simpleName}.hybrid_corrected.fastq
    gzip ${ont_reads.simpleName}.hybrid_corrected.fastq
    """
}

workflow {
    // Input validation
    if (!params.input_ont && !params.input_fast5) {
        error "Please provide either ONT FASTQ (--input_ont) or FAST5 (--input_fast5) files"
    }

    // Paired-end Illumina reads: expect input_short to be a directory with both R1 and R2 files
    if (params.input_short) {
        illumina_reads = Channel.fromFilePairs(
            "${params.input_short}/${params.pattern_short}",
            flat: true,
            checkIfExists: true
        )
    } else {
        illumina_reads = Channel.empty()
    }

    // Main workflow
    if (params.input_fast5) {
        ont_reads = DORADO(params.input_fast5)
    } else {
        ont_reads = Channel.fromPath(params.input_ont)
    }

    trimmed_reads = PORECHOP(ont_reads)
    filtered_reads = NANOFILT(trimmed_reads)

    if (params.hybrid && params.input_short) {
        illumina_processed = FASTP(illumina_reads)
        corrected_reads = FMLRC2(filtered_reads, illumina_processed.filtered_pe)
    } else {
        corrected_reads = filtered_reads
    }

    emit:
    filtered_fastq = filtered_reads
    corrected_fastq = corrected_reads
    short_reads = params.hybrid ? illumina_processed.filtered_pe : null
} 