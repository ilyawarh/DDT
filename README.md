# Metagenomics Nextflow Pipeline

## Overview
A modular, profile-driven Nextflow pipeline for metagenomic analysis of Nanopore (R10) and hybrid Nanopore+Illumina data. Supports raw signal (.fast5) or FASTQ input, with optional paired-end Illumina reads for hybrid assembly. Generates comprehensive reports on taxonomic composition, genome reconstruction, virulence factors, and more.

## Modules
- **Preprocess**: Adapter trimming, quality filtering, correction (ONT/Illumina)
- **Read_Based**: Taxonomic classification (KrakenUniq, Kaiju, etc.)
- **Assembly**: Assembly and polishing (metaFlye, Unicycler, etc.)
- **Binning**: MAG binning (SemiBin2, MetaBinner, etc.)
- **Annotation**: Taxonomic/functional annotation, virulence/plasmid/viral search
- **Evaluation**: Genome quality and abundance assessment

## Profiles
- `hybrid`: Enables hybrid assembly/correction
- `annotate_more`: Adds extra annotation tools
- `poly_binning`: Parallel binning and refinement
- `read_expand`: Alternative classifiers
- `basecall`: Dorado basecalling for raw signals
- `slurm`: Cluster execution

## Input Parameters
- `--input_ont`: Path to Nanopore FASTQ directory
- `--input_fast5`: Path to Nanopore .fast5 directory
- `--short_1`, `--short_2`: Illumina paired-end FASTQ directories
- `--databases`: Path to folder with all databases
- All other parameters set in `nextflow.config`

## QC Shortcut
Run QC on all reads before main pipeline:
```bash
nextflow run qc.nf --input_ont <ONT_FASTQ_DIR> --short_1 <ILLUMINA_R1_DIR> --short_2 <ILLUMINA_R2_DIR> --outdir reports/QC
```

## Example Pipeline Run
```bash
nextflow run main.nf --input_ont data/ONT --short_1 data/Illumina/R1 --short_2 data/Illumina/R2 --databases db/ --profile hybrid,slurm
```

## Output Structure
- `results/`: Intermediate and final results
- `reports/`: All reports (MultiQC, QUAST, Taxpasta, final summary, etc.)
- `raw_data/`: Symlinks and raw/processed sequencing files

## Final Report
- `final_report.tsv` and `final_report.xlsx`: Taxonomic composition, quantitative proportions, genome sizes/labels, virulence factors per genome
- Other detailed reports (functional annotation, plasmids, viral genomes) in `reports/`

## Help
To print help from the command line:
```bash
nextflow run main.nf --help
``` 