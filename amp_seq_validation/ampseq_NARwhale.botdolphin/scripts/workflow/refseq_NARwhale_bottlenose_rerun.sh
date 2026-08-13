#!/usr/bin/env bash
set -euo pipefail

RUN_ID="$(date +%Y%m%d_%H%M%S)"
OUTDIR="refseq_rightwhale_bottlenose_rerun_${RUN_ID}"
REFDIR="${OUTDIR}/reference_genomes"

mkdir -p "$REFDIR"

############################################
# North Atlantic right whale RefSeq: GCF_028564815.1
# Species: Eubalaena glacialis
# Assembly: mEubGla1.1.hap2.+ XY
# FTP path encodes "+ XY" as "._XY"
############################################

RIGHT_WHALE_REFSEQ_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/564/815/GCF_028564815.1_mEubGla1.1.hap2._XY/GCF_028564815.1_mEubGla1.1.hap2._XY_genomic.fna.gz"
RIGHT_WHALE_REFSEQ_GZ="${REFDIR}/GCF_028564815.1_mEubGla1.1.hap2._XY_genomic.fna.gz"
RIGHT_WHALE_REFSEQ="${REFDIR}/GCF_028564815.1_mEubGla1.1.hap2._XY_genomic.fna"

wget -O "$RIGHT_WHALE_REFSEQ_GZ" "$RIGHT_WHALE_REFSEQ_URL"
gunzip "$RIGHT_WHALE_REFSEQ_GZ"
samtools faidx "$RIGHT_WHALE_REFSEQ"

grep "NC_083728.1" "$RIGHT_WHALE_REFSEQ"   # parent PCNA gene
grep "NC_083716.1" "$RIGHT_WHALE_REFSEQ"   # PCNA-like candidate

# Parent PCNA gene from RefSeq GFF:
# NC_083728.1 Gnomon gene 38362236-38368011, strand -
# Strand: Plus/Minus, so reverse complement is needed.
samtools faidx --reverse-complement "$RIGHT_WHALE_REFSEQ" \
    NC_083728.1:38362236-38368011 > "${OUTDIR}/rightwhale_parent_reference.fna"

# PCNA-like candidate from RefSeq GFF:

# Strand: Plus/Minus, so reverse complement is needed.
samtools faidx --reverse-complement "$RIGHT_WHALE_REFSEQ" \
    NC_083716.1:239260552-239261337 > "${OUTDIR}/rightwhale_pcna_like_candidate.fna"

# PCNA-like candidate with 150 bp flanks for B primer-style reference
samtools faidx --reverse-complement "$RIGHT_WHALE_REFSEQ" \
    NC_083716.1:239260402-239261487 > "${OUTDIR}/rightwhale_pcna_like_candidateB_reference.fna"

############################################
# Bottlenose dolphin RefSeq: GCF_011762595.2
# Species: Tursiops truncatus
# Assembly: mTurTru1.mat.Y
############################################

BOTTLENOSE_REFSEQ_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/762/595/GCF_011762595.2_mTurTru1.mat.Y/GCF_011762595.2_mTurTru1.mat.Y_genomic.fna.gz"
BOTTLENOSE_REFSEQ_GZ="${REFDIR}/GCF_011762595.2_mTurTru1.mat.Y_genomic.fna.gz"
BOTTLENOSE_REFSEQ="${REFDIR}/GCF_011762595.2_mTurTru1.mat.Y_genomic.fna"

wget -O "$BOTTLENOSE_REFSEQ_GZ" "$BOTTLENOSE_REFSEQ_URL"
gunzip "$BOTTLENOSE_REFSEQ_GZ"
samtools faidx "$BOTTLENOSE_REFSEQ"

grep "NC_047048.1" "$BOTTLENOSE_REFSEQ"   # parent PCNA gene
grep "NC_047038.1" "$BOTTLENOSE_REFSEQ"   # PCNA pseudogene 1
grep "NC_047040.1" "$BOTTLENOSE_REFSEQ"   # PCNA pseudogene 2

# Parent PCNA gene from RefSeq GFF:
# Strand: Plus/Plus, so no reverse complement needed.
samtools faidx "$BOTTLENOSE_REFSEQ" \
    NC_047048.1:55092925-55098911 > "${OUTDIR}/bottlenose_parent_reference.fna"

# PCNA pseudogene 1 from RefSeq GFF:
# NC_047038.1 Gnomon pseudogene 95022203-95023028, strand -
# In the bottlenose RefSeq GFF, the two PCNA-derived loci are annotated as pseudogene spans without exon subfeatures, so the full pseudogene span is used as the retrogene-candidate body.
# Strand: Plus/Minus, so reverse complement is needed.
samtools faidx --reverse-complement "$BOTTLENOSE_REFSEQ" \
    NC_047038.1:95022203-95023028 > "${OUTDIR}/bottlenose_pcna_pseudogene1.fna"

# PCNA pseudogene 2 from RefSeq GFF:
# Strand: Plus/Plus, so no reverse complement needed.
samtools faidx "$BOTTLENOSE_REFSEQ" \
    NC_047040.1:1531464-1533838 > "${OUTDIR}/bottlenose_pcna_pseudogene2.fna"

# PCNA pseudogene 1 with 150 bp flanks for B primer-style reference
samtools faidx --reverse-complement "$BOTTLENOSE_REFSEQ" \
    NC_047038.1:95022053-95023178 > "${OUTDIR}/bottlenose_pcna_pseudogene1B_reference.fna"

# PCNA pseudogene 2 with 150 bp flanks for B primer-style reference
samtools faidx "$BOTTLENOSE_REFSEQ" \
    NC_047040.1:1531314-1533988 > "${OUTDIR}/bottlenose_pcna_pseudogene2B_reference.fna"

############################################
# Rename headers to match amplicon workflow
############################################

sed -i 's/>.*/>NorthAtlanticRightWhalePCNA_AParent/' "${OUTDIR}/rightwhale_parent_reference.fna"
sed -i 's/>.*/>NorthAtlanticRightWhalePCNA_like_candidate/' "${OUTDIR}/rightwhale_pcna_like_candidate.fna"
sed -i 's/>.*/>NorthAtlanticRightWhalePCNA_like_candidateB/' "${OUTDIR}/rightwhale_pcna_like_candidateB_reference.fna"

sed -i 's/>.*/>BottlenoseDolphinPCNA_AParent/' "${OUTDIR}/bottlenose_parent_reference.fna"
sed -i 's/>.*/>BottlenoseDolphinPCNA_pseudogene1/' "${OUTDIR}/bottlenose_pcna_pseudogene1.fna"
sed -i 's/>.*/>BottlenoseDolphinPCNA_pseudogene2/' "${OUTDIR}/bottlenose_pcna_pseudogene2.fna"
sed -i 's/>.*/>BottlenoseDolphinPCNA_pseudogene1B/' "${OUTDIR}/bottlenose_pcna_pseudogene1B_reference.fna"
sed -i 's/>.*/>BottlenoseDolphinPCNA_pseudogene2B/' "${OUTDIR}/bottlenose_pcna_pseudogene2B_reference.fna"

############################################
# Sanity checks
############################################

echo
echo "Output FASTA headers:"
grep ">" "${OUTDIR}"/*.fna

echo
echo "Output sequence lengths:"
printf "rightwhale_parent\t"
grep -v ">" "${OUTDIR}/rightwhale_parent_reference.fna" | tr -d '\n' | wc -c
printf "rightwhale_pcna_like_candidate\t"
grep -v ">" "${OUTDIR}/rightwhale_pcna_like_candidate.fna" | tr -d '\n' | wc -c
printf "rightwhale_pcna_like_candidateB\t"
grep -v ">" "${OUTDIR}/rightwhale_pcna_like_candidateB_reference.fna" | tr -d '\n' | wc -c
printf "bottlenose_parent\t"
grep -v ">" "${OUTDIR}/bottlenose_parent_reference.fna" | tr -d '\n' | wc -c
printf "bottlenose_pcna_pseudogene1\t"
grep -v ">" "${OUTDIR}/bottlenose_pcna_pseudogene1.fna" | tr -d '\n' | wc -c
printf "bottlenose_pcna_pseudogene2\t"
grep -v ">" "${OUTDIR}/bottlenose_pcna_pseudogene2.fna" | tr -d '\n' | wc -c
printf "bottlenose_pcna_pseudogene1B\t"
grep -v ">" "${OUTDIR}/bottlenose_pcna_pseudogene1B_reference.fna" | tr -d '\n' | wc -c
printf "bottlenose_pcna_pseudogene2B\t"
grep -v ">" "${OUTDIR}/bottlenose_pcna_pseudogene2B_reference.fna" | tr -d '\n' | wc -c

echo
echo "Done. New outputs are in: $OUTDIR"
