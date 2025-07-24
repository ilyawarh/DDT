#!/usr/bin/env nextflow

// Read-based taxonomic classification module

process KRAKENUNIQ {
    tag { "taxonomic classification with KrakenUniq" }
    publishDir "${params.outdir}/taxonomy/krakenuniq", mode: 'copy'

    when:
    !params.read_expand

    input:
    path reads

    output:
    path "*.krakenuniq.tsv", emit: report
    path "*.kraken", emit: classifications

    script:
    """
    krakenuniq \
        --db ${params.databases}/krakenuniq \
        --threads ${task.cpus} \
        --report-file ${reads.simpleName}.krakenuniq.tsv \
        --output ${reads.simpleName}.kraken \
        ${reads}
    """
}

process KAIJU {
    tag { "taxonomic classification with Kaiju" }
    publishDir "${params.outdir}/taxonomy/kaiju", mode: 'copy'

    when:
    !params.read_expand

    input:
    path reads

    output:
    path "*.kaiju.tsv", emit: report
    path "*.kaiju.names", emit: names

    script:
    """
    kaiju \
        -t ${params.databases}/kaiju/nodes.dmp \
        -f ${params.databases}/kaiju/kaiju_db.fmi \
        -i ${reads} \
        -o ${reads.simpleName}.kaiju.tsv

    kaiju-addTaxonNames \
        -t ${params.databases}/kaiju/nodes.dmp \
        -n ${params.databases}/kaiju/names.dmp \
        -i ${reads.simpleName}.kaiju.tsv \
        -o ${reads.simpleName}.kaiju.names
    """
}

process CENTRIFUGE {
    tag { "taxonomic classification with Centrifuge" }
    publishDir "${params.outdir}/taxonomy/centrifuge", mode: 'copy'

    when:
    params.read_expand

    input:
    path reads

    output:
    path "*.centrifuge.tsv", emit: report
    path "*.centrifuge.kreport", emit: kreport

    script:
    """
    centrifuge \
        -x ${params.databases}/centrifuge/nt \
        -U ${reads} \
        --threads ${task.cpus} \
        -S ${reads.simpleName}.centrifuge.tsv \
        --report-file ${reads.simpleName}.centrifuge.kreport
    """
}

process DIAMOND {
    tag { "protein search with DIAMOND" }
    publishDir "${params.outdir}/taxonomy/diamond", mode: 'copy'

    when:
    params.read_expand

    input:
    path reads

    output:
    path "*.diamond.daa", emit: daa

    script:
    """
    diamond blastx \
        --db ${params.databases}/nr \
        --query ${reads} \
        --out ${reads.simpleName}.diamond.daa \
        --threads ${task.cpus} \
        --outfmt 100 \
        --sensitive
    """
}

process MEGAN_LR {
    tag { "taxonomic classification with MEGAN-LR" }
    container 'megan-lr:latest'  // Replace with actual Docker image
    publishDir "${params.outdir}/taxonomy/megan", mode: 'copy'

    when:
    params.read_expand

    input:
    path daa

    output:
    path "*.megan.txt", emit: classifications
    path "*.rma6", emit: rma

    script:
    """
    daa2rma \
        --in ${daa} \
        --out ${daa.simpleName}.rma6 \
        --reads ${daa.simpleName}.megan.txt \
        --mapDB ${params.databases}/megan/megan-map-Jan2021.db
    """
}

process MERGE_TAXONOMY {
    tag { "merging taxonomic classifications" }
    publishDir "${params.outdir}/taxonomy", mode: 'copy'

    input:
    path '*'

    output:
    path "merged_taxonomy.tsv", emit: merged

    script:
    """
    # Merge taxonomic classifications from different tools
    # This is a placeholder - implement actual merging logic
    echo "Merged taxonomy" > merged_taxonomy.tsv
    cat * >> merged_taxonomy.tsv
    """
}

workflow {
    // Input validation
    if (!params.databases) {
        error "Please provide the path to taxonomy databases (--databases)"
    }

    // Main workflow
    reads = Channel.fromPath(params.input_ont)

    if (params.read_expand) {
        // Extended taxonomic classification
        centrifuge_results = CENTRIFUGE(reads)
        diamond_results = DIAMOND(reads)
        megan_results = MEGAN_LR(diamond_results.daa)
        tax_results = centrifuge_results.mix(megan_results.classifications)
    } else {
        // Standard taxonomic classification
        krakenuniq_results = KRAKENUNIQ(reads)
        kaiju_results = KAIJU(reads)
        tax_results = krakenuniq_results.mix(kaiju_results.report)
    }

    merged_taxonomy = MERGE_TAXONOMY(tax_results.collect())

    emit:
    taxonomy = merged_taxonomy.merged
} 