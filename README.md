# BINF*6999 Summer Research

Bioinformatics research investigating copy number variation in cancer-associated genes and its potential relationship with cancer resistance.

## Overall research question

Does variation in retrogene copy number contribute to increased cancer resistance? This project explores patterns across species and cancer-associated genes, with the broader goal of investigating the mechanisms that may underlie those patterns.

- **Phase 1:** Investigate the **proliferating cell nuclear antigen (PCNA)** parent gene and retrogene copies across cetacean species.
- **Phase 2:** Extend the comparative approach: additional cancer-associated genes and species.

> **Project status:** Phase 1 analyses and reproducibility documentation are being refined. Work toward the Phase 2 objectives is planned to begin in September 2026.

## Project components

This repository contains three related analytical components:

1. a survey of PCNA loci and publicly available genomic and RNA-sequencing resources across cetaceans
2. validation of a long-range PCR and Oxford Nanopore amplicon-sequencing approach for characterizing selected PCNA copies
3. comparative analysis of PCNA parent-gene and retrogene sequences from selected cetacean species

## Goals

- Identify PCNA parent genes and putative retrogene copies in cetacean genomes.
- Compile publicly available cetacean genome and RNA-seq resources for downstream analysis.
- Evaluate the use of long-range PCR followed by Oxford Nanopore amplicon sequencing to isolate and characterize individual PCNA copies.
- Compare the number, sequence identity, and genomic organization of PCNA copies among species.
- Develop a reproducible framework that can be extended to additional cancer-associated genes.

## Repository structure

```text
.
├── README.md
├── environment.yml
├── amp_seq_validation/
│   ├── README.md
│   ├── ampseq_NARwhale.botdolphin/
│   │   ├── inputs/
│   │   ├── references/
│   │   ├── scripts/
│   │   └── workflow_runs/
│   ├── scripts/
│   ├── reference sequences and genome assemblies
│   └── sample sheets and generated indexes
├── cetacean_PCNA_project/
│   ├── annotations/
│   ├── blast_db/
│   ├── blast_results/
│   ├── genomes/
│   ├── pcna_refs/
│   ├── proteins/
│   ├── rna_seq/
│   ├── scripts/
│   ├── tables/
│   └── tools/
├── pcna_retrogene_analyses/
│   ├── comparison sequence sets
│   ├── R analysis scripts
│   └── results/
│       ├── alignments
│       ├── pairwise-comparison tables
│       ├── phylogenetic trees
│       ├── identity heatmaps
│       └── poster figures and source data
└── project_documents/
    ├── LLM_validation_log.md
    ├── poster_srp_fin.pdf
    └── references/
        └── binf6999_references.csv
```

## Component descriptions

### `cetacean_PCNA_project/`

Survey and evaluation of publicly available cetacean genomic and RNA-sequencing resources, primarily obtained through NCBI.

Key contents include:

- genome assemblies, annotations, and protein sequences
- PCNA query sequences and locally generated BLAST databases
- BLAST results for candidate PCNA loci
- inventories of available cetacean genomes and RNA-seq datasets
- scripts for collecting and annotating species-level metadata
- summary tables of candidate PCNA loci and public datasets

### `amp_seq_validation/`

Validation of an experimental approach that combines long-range PCR with Oxford Nanopore amplicon sequencing to characterize selected PCNA parent genes and retrogene copies.

Initial work focused on beluga (*Delphinapterus leucas*) and blue whale (*Balaenoptera musculus*), followed by analyses of North Atlantic right whale (*Eubalaena glacialis*) and bottlenose dolphin (*Tursiops truncatus*) samples.

Directory contains:

- genome-derived and amplicon-specific reference sequences
- sample sheets and primer information
- scripts used to launch or inspect amplicon workflows
- `wf-amplicon` run parameters and generated outputs
- consensus sequences, alignments, and variant calls
- read-level investigations of allele frequencies, mixed sites, and phasing.

See [`amp_seq_validation/README.md`](amp_seq_validation/README.md) for the detailed biological background and analysis log.

### `pcna_retrogene_analyses/`

Comparative sequence analyses performed after candidate PCNA sequences and amplicon consensus sequences were identified.

It contains:

- curated sequence sets used for PCNA comparisons
- R scripts for preparing and comparing homologous sequence regions
- pairwise sequence-identity tables
- multiple-sequence alignments
- neighbour-joining trees
- parent-gene and retrogene identity heatmaps
- final poster figures with their source data

### `project_documents/`

Project-level supporting material, including the research poster, references consulted, and a log documenting validation of material LLM-assisted work.

## Analysis overview

Several analyses are exploratory, and some scripts may still contain run-specific paths or assumptions. Script parameters and input paths are being reviewed to support execution in a clean environment.

## Software and environment

Software used across the project includes:

- BLAST+
- Nextflow
- EPI2ME Labs `wf-amplicon`
- Minimap2
- SAMtools
- R
- Python

The root [`environment.yml`](environment.yml) provides a starting Conda environment for the command-line, Python, and R dependencies used in the project.

Create and activate the environment with:

```bash
conda env create -f environment.yml
conda activate binf6999
```

## Data and generated files

The repository contains both source material and selected generated outputs. Generated files may include:

- BLAST database indexes (`.ndb`, `.nhr`, `.nin`, `.nsq`, and related files);
- FASTA indexes (`.fai`);
- BAM alignments and indexes;
- compressed variant files and indexes;
- Nextflow logs and execution reports; and
- intermediate workflow output directories.

Generated files are being reviewed as the repository is reorganized. Files that can be recreated reliably may be removed from version control once their generating commands have been validated and documented. Large raw sequencing datasets are stored outside this Git repository.

## Data sources and attribution

Genome assemblies, annotations, and sequencing-study metadata were obtained from public NCBI resources. sequenced samples were from samples obtained by our lab. Relevant assembly accessions and sample information are recorded in corresponding analysis files and notes. Academic sources consulted during the project are listed in `project_documents/references/binf6999_references.csv`.

## Reproducibility

Reproducibility documentation is in progress. currently an environment.yml has been created to support the dependencies used in the ongoing work of this project

