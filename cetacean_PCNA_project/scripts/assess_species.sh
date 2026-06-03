#!/usr/bin/env bash
# assess_species.sh
# Purpos: to check genome assemblies + RNA-seq data for each cetacean species on NCBI.
#
# Usage:
#   ./scripts/assess_species.sh "scientific name" commonname
#
# Outputs (in tables/):
#   beluga_genomes.tsv   — one row per genome assembly
#   beluga_rnaseq.csv    — one row per RNA-seq run

TAXON="$1"    # scientific name  e.g. "Delphinapterus leucas"
COMMON="$2"   # common name      e.g. beluga  (used for output filenames)

if [[ -z "$TAXON" || -z "$COMMON" ]]; then
    echo "Usage: $0 \"Scientific name\" commonname"
    echo "  e.g. $0 \"Delphinapterus leucas\" beluga"
    exit 1
fi

mkdir -p tables

echo ""
echo "=== $COMMON ($TAXON) ==="


# 1. GENOME ASSEMBLIES  (uses tools/datasets)

echo ""
echo "Fetching genome assemblies from NCBI..."

GENOME_OUT="tables/${COMMON}_genomes.tsv"

./tools/datasets summary genome taxon "$TAXON" --as-json-lines \
  | ./tools/dataformat tsv genome \
    --fields accession,organism-name,assminfo-level,assminfo-name,assminfo-refseq-category,assmstats-contig-n50,assmstats-scaffold-n50,assmstats-total-sequence-len,assmstats-number-of-scaffolds,annotinfo-name,assminfo-release-date \
  > "$GENOME_OUT"

N_GENOMES=$(tail -n +2 "$GENOME_OUT" | grep -c '.') # subtract header
echo "  Found $N_GENOMES assembly/assemblies → $GENOME_OUT"


# 2. RNA-SEQ RUNS  (uses esearch + efetch)

echo ""
echo "Fetching RNA-seq run info from NCBI SRA..."

RNASEQ_OUT="rna_seq/${COMMON}_rnaseq_runinfo.csv"

esearch -db sra -query "\"${TAXON}\"[Organism] AND \"RNA-Seq\"[Strategy]" \
  | efetch -format runinfo \
  > "$RNASEQ_OUT"

N_RUNS=$(( $(wc -l < "$RNASEQ_OUT") - 1 ))
echo "  Found $N_RUNS RNA-seq run(s) → $RNASEQ_OUT"


# 3. QUICK SUMMARY to terminal

echo ""
echo "--- Genome assembly levels ---"
awk -F'\t' 'NR>1 {print $3}' "$GENOME_OUT" | sort | uniq -c | sort -rn

echo ""
echo "--- RNA-seq platforms ---"

head -1 "$RNASEQ_OUT" | tr ',' '\n' | grep -n "Platform"   # shows you the real column number
awk -F',' 'NR==1{for(i=1;i<=NF;i++) if($i=="Platform") col=i} NR>1{print $col}' \
    "$RNASEQ_OUT" | sort | uniq -c | sort -rn

echo "--- RNA-seq BioProjects ---"
awk -F',' 'NR==1{for(i=1;i<=NF;i++) if($i=="BioProject") col=i} NR>1{print $col}' \
    "$RNASEQ_OUT" | sort -u

echo ""
echo "Done! can open these in Excel:"
echo "  $GENOME_OUT"
echo "  $RNASEQ_OUT"
