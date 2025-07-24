#!/usr/bin/env nextflow

// Assembly and polishing module

process METAFLYE {
    tag { "metaFlye assembly" }
    publishDir "${params.outdir}/assembly/flye", mode: 'copy'

    when:
    !params.raven

    input:
    path reads

    output:
    path "*.flye.fasta", emit: assembly
    path "flye_info/*", emit: info

    script:
    """
    flye \
        --nano-raw ${reads} \
        --out-dir flye_info \
        --threads ${task.cpus} \
        --meta \
        ${params.flye_opts}
    
    cp flye_info/assembly.fasta ${reads.simpleName}.flye.fasta
    """
}

process RAVEN {
    tag { "Raven assembly" }
    publishDir "${params.outdir}/assembly/raven", mode: 'copy'

    when:
    params.raven

    input:
    path reads

    output:
    path "*.raven.fasta", emit: assembly

    script:
    """
    raven \
        --threads ${task.cpus} \
        ${params.raven_opts} \
        ${reads} > ${reads.simpleName}.raven.fasta
    """
}

process UNICYCLER {
    tag { "hybrid assembly" }
    publishDir "${params.outdir}/assembly/unicycler", mode: 'copy'

    when:
    params.hybrid

    input:
    path long_reads
    tuple path(short_1), path(short_2)

    output:
    path "*.unicycler.fasta", emit: assembly
    path "unicycler_info/*", emit: info

    script:
    """
    unicycler \
        -l ${long_reads} \
        -1 ${short_1} \
        -2 ${short_2} \
        -o unicycler_info \
        -t ${task.cpus} \
        ${params.unicycler_opts}
    
    cp unicycler_info/assembly.fasta ${long_reads.simpleName}.unicycler.fasta
    """
}

process MINIMAP2_ALIGN {
    tag { "read mapping" }
    publishDir "${params.outdir}/assembly/mapping", mode: 'copy'

    input:
    path assembly
    path reads

    output:
    path "*.bam", emit: bam
    path "*.bam.bai", emit: bai

    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} ${assembly} ${reads} | \
        samtools sort -@ ${task.cpus} -o ${assembly.simpleName}.bam
    samtools index ${assembly.simpleName}.bam
    """
}

process RACON {
    tag { "Racon polishing" }
    publishDir "${params.outdir}/assembly/racon", mode: 'copy'

    input:
    path assembly
    path reads
    path bam

    output:
    path "*.racon.fasta", emit: polished

    script:
    """
    racon \
        -t ${task.cpus} \
        ${reads} \
        ${bam} \
        ${assembly} > ${assembly.simpleName}.racon.fasta
    """
}

process MEDAKA {
    tag { "Medaka polishing" }
    container 'nanoporetech/medaka:latest'
    publishDir "${params.outdir}/assembly/medaka", mode: 'copy'

    input:
    path assembly
    path reads

    output:
    path "*.medaka.fasta", emit: polished

    script:
    """
    medaka_consensus \
        -i ${reads} \
        -d ${assembly} \
        -o medaka_output \
        -t ${task.cpus} \
        -m r941_min_sup_g507

    cp medaka_output/consensus.fasta ${assembly.simpleName}.medaka.fasta
    """
}

process PILON {
    tag { "Pilon polishing" }
    publishDir "${params.outdir}/assembly/pilon", mode: 'copy'

    when:
    params.hybrid

    input:
    path assembly
    tuple path(short_1), path(short_2)

    output:
    path "*.pilon.fasta", emit: polished

    script:
    """
    # Map short reads to assembly
    bwa index ${assembly}
    bwa mem -t ${task.cpus} ${assembly} ${short_1} ${short_2} | \
        samtools sort -@ ${task.cpus} -o illumina_mapped.bam
    samtools index illumina_mapped.bam

    # Run Pilon
    pilon \
        --genome ${assembly} \
        --frags illumina_mapped.bam \
        --output ${assembly.simpleName}.pilon \
        --changes \
        --threads ${task.cpus}

    mv ${assembly.simpleName}.pilon.fasta ${assembly.simpleName}.pilon.fasta
    """
}

process HOMOPOLISH {
    tag { "Homopolish polishing" }
    publishDir "${params.outdir}/assembly/homopolish", mode: 'copy'

    input:
    path assembly

    output:
    path "*.homopolished.fasta", emit: polished

    script:
    """
    homopolish polish \
        -a ${assembly} \
        -o ${assembly.simpleName}.homopolished.fasta \
        -m R941-sup \
        -d ${params.databases}/homopolish
    """
}

workflow {
    // Input validation
    if (!params.input_ont && !params.input_fast5) {
        error "Please provide either ONT FASTQ (--input_ont) or FAST5 (--input_fast5) files"
    }

    // Main workflow
    reads = Channel.fromPath(params.input_ont)

    // Assembly
    if (params.hybrid && params.short_1 && params.short_2) {
        short_reads = tuple(file(params.short_1), file(params.short_2))
        assembly = UNICYCLER(reads, short_reads)
    } else if (params.raven) {
        assembly = RAVEN(reads)
    } else {
        assembly = METAFLYE(reads)
    }

    // Polishing
    mapped_reads = MINIMAP2_ALIGN(assembly.assembly, reads)
    racon_polished = RACON(assembly.assembly, reads, mapped_reads.bam)
    medaka_polished = MEDAKA(racon_polished.polished, reads)

    if (params.hybrid) {
        short_reads = tuple(file(params.short_1), file(params.short_2))
        final_assembly = PILON(medaka_polished.polished, short_reads)
    } else {
        final_assembly = HOMOPOLISH(medaka_polished.polished)
    }

    emit:
    assembly = final_assembly.polished
    info = assembly.info
} 