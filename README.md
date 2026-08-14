# BINF*6999 Summer Research

Bioinformatics research investigating copy number variation.

## Overall research question

Does variation in retrogene copy number contribute to increased cancer resistance? What are the connections across species, across different cancer-linked genes, investigate mechanisms underlying it.

Phase 1: focuses on the **proliferating cell nuclear antigen (PCNA)** gene across cetacean species
Phase 2: extends the comparison

## This repository brings together three related components of the project:

1. a comparative survey of PCNA loci and available genomic and RNA-sequencing resources across cetaceans
2. validation of amplicon-sequencing analyses 
3. characterization and comparison of PCNA copies in selected species.

> **Project status:** Analyses targetign phase 2 objectives to begin in September, 2026. File organization and reproducibility documentation are being refined.

## Goals

- Identify PCNA parent genes and putative retrogene copies in cetacean genomes.
- Compile publicly available cetacean genome and RNA-seq resources for downstream analysis.
- Validate experimental method of using long-range PCR to select retrogene copies followed by  Oxford Nanopore amplicon sequencing on PCR products
- Compare the number, sequence identity, and genomic organization of PCNA copies among species.

## Repository structure

```text
.
├── amp_seq_validation/
│   ├── README.md
│   ├── ampseq_NARwhale.botdolphin/
│   ├── reference sequences and genome assemblies
│   ├── sample sheets
│   ├── shell and R analysis scripts
│   └── generated workflow outputs
└── cetacean_PCNA_project/
    ├── annotations/
    ├── blast_db/
    ├── blast_results/
    ├── genomes/
    ├── pcna_refs/
    ├── proteins/
    ├── rna_seq/
    ├── scripts/
    ├── tables/
    └── tools/

```

### `cetacean_PCNA_project/`

Evaluation of publically available data, primary resource NCBI

Key contents include:

- genome assemblies, annotations, and protein sequences
- PCNA query sequences and locally generated BLAST databases
- BLAST search results for candidate PCNA loci
- inventories of available cetacean genomes and RNA-seq datasets
- scripts for assembling and annotating species-level metadata
- summary tables of candidate PCNA loci

### `amp_seq_validation/`

Amplicon-sequencing validation and sequence comparison of PCNA parent genes and retrogene copies.

The initial work focused on beluga (*Delphinapterus leucas*) and blue whale (*Balaenoptera musculus*), with subsequent analyses of North Atlantic right whale (*Eubalaena glacialis*) and bottlenose dolphin (*Tursiops truncatus*) samples. The directory contains reference sequences, sample sheets, workflow scripts, consensus sequences, alignments, variant calls, pairwise comparisons and figures generated.

See [`amp_seq_validation/README.md`](amp_seq_validation/README.md) for the detailed biological background and analysis log.

## Analysis overview

The project uses a combination of:

- NCBI genome assemblies and annotations
- nucleotide and protein BLAST searches
- publicly available RNA-seq metadata
- Oxford Nanopore amplicon sequencing
- the EPI2ME Labs `wf-amplicon` Nextflow workflow conducted on the application interface
- read alignment and consensus-sequence analysis
- R-based sequence comparison, visualization, and exploratory analysis.

Because several analyses are exploratory, individual scripts may currently contain run-specific paths or assumptions. Script parameters and input paths before are currently being reviewed and modified to support running them in a new environment.

## Software

Software used in different parts of the project includes:

- BLAST+
- Nextflow
- EPI2ME Labs `wf-amplicon`
- SAMtools
- R
- Python

## Data and generated files

This repository contains a mixture of source files and generated outputs. Examples of generated files include:

- BLAST database indexes (`.ndb`, `.nhr`, `.nin`, `.nsq`, and related files);
- FASTA indexes (`.fai`);
- BAM alignment files and indexes;
- Nextflow logs and execution reports; and
- intermediate workflow output directories.

These files may be removed from version control in a future cleanup once scripts and instructions for regenerating them have been thoroughly validated. Large raw sequencing datasets are stored in an appropriate data repository rather than committed directly to Git.

## Reproducibility

Reproducibility documentation is in progress. currently an environment.yml has been created to support the dependencies used in the ongoing work of this project

## Citation and data sources

Genome assemblies, annotations, and sequencing-study metadata were obtained from public NCBI resources. Samples sequenced are from samples obtained by our lab. Relevant assembly accessions are recorded in corresponding analysis notes. Academic papers consulted in the course of this analysis are available in the references folder


