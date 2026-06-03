# download GCF RefSeq assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/288/925/GCF_002288925.2_ASM228892v3/GCF_002288925.2_ASM228892v3_ge
nomic.fna.gz

# Decompress
gunzip GCF_002288925.2_ASM228892v3_genomic.fna.gz

# Index
samtools faidx GCF_002288925.2_ASM228892v3_genomic.fna
grep ">" GCF_002288925.2_ASM228892v3_genomic.fna | head -50
# Ensuring that scaffolds of interest are present
grep "NW_022098045" GCF_002288925.2_ASM228892v3_genomic.fna
grep "NW_022098076" GCF_002288925.2_ASM228892v3_genomic.fna
grep "NW_022098051" GCF_002288925.2_ASM228892v3_genomic.fna

# Extraction of each sequence
# Parent gene — plus strand, no reverse complement needed
samtools faidx GCF_002288925.2_ASM228892v3_genomic.fna \
    NW_022098051.1:7066550-7072506 > parent_reference.fna

# Retrogene 1 — minus strand, reverse complement required
# Also add 150bp flanks for B primer set
samtools faidx --reverse-complement GCF_002288925.2_ASM228892v3_genomic.fna \
    NW_022098045.1:8167847-8168931 > retrogene1B_reference.fna

# Retrogene 2 — minus strand, reverse complement required
# Also add 150bp flanks for B primer set
samtools faidx --reverse-complement GCF_002288925.2_ASM228892v3_genomic.fna \
    NW_022098076.1:94205968-94207053 > retrogene2B_reference.fna

# for ARetrogenes mixed sample:
# Retrogene 1 exon body only 
samtools faidx --reverse-complement GCF_002288925.2_ASM228892v3_genomic.fna \
    NW_022098045.1:8167997-8168781 > retrogene1_exons.fna
# Retrogene 2 exon body only
samtools faidx --reverse-complement GCF_002288925.2_ASM228892v3_genomic.fna \
    NW_022098076.1:94206118-94206903 > retrogene2_exons.fna

# sanity-checking the output
# Check headers and sequence lengths
grep ">" *.fna
grep -v ">" parent_reference.fna | tr -d '\n' | wc -c      # expect ~4554
grep -v ">" retrogene1_exons.fna | tr -d '\n' | wc -c      # expect ~786
grep -v ">" retrogene2_exons.fna | tr -d '\n' | wc -c      # expect ~786
grep -v ">" retrogene1B_reference.fna | tr -d '\n' | wc -c # expect ~1086
grep -v ">" retrogene2B_reference.fna | tr -d '\n' | wc -c # expect ~1086

# renaming to match amplicon workflow
sed -i 's/>.*/>B12BelugaPCNA_AParent/' parent_reference.fna
sed -i 's/>.*/>B12BelugaPCNA_Retrogene1_exons/' retrogene1_exons.fna
sed -i 's/>.*/>B12BelugaPCNA_Retrogene2_exons/' retrogene2_exons.fna
sed -i 's/>.*/>B12BelugaRetrogene1B/' retrogene1B_reference.fna
sed -i 's/>.*/>B12BelugaRetrogene2B/' retrogene2B_reference.fna

# Repeating the process for blue whale
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/873/245/GCF_009873245.2_mBalMus1.pri.v3/GCF_009873245.2_mBalMus1.pri.v3_genomic.fna.gz

gunzip GCF_009873245.2_mBalMus1.pri.v3_genomic.fna.gz

samtools faidx GCF_009873245.2_mBalMus1.pri.v3_genomic.fna

# verify scaffold acessions from blast results are present
grep "NC_041915028" GCF_009873245.2_mBalMus1.pri.v3_genomic.fna
grep "NC_045799" GCF_009873245.2_mBalMus1.pri.v3_genomic.fna
# retrogene not placed - troubleshooting

# from the blast results, 1 retrogene: at NC_041915028.1:111914528–111915313 on the minus strand
# extracting exon body only for BlueWhaleRetrogenes, and another extraction with flanking regions for BlueWhaleRetrogene2B

# For BlueWhaleRetrogenes — exon body only
samtools faidx --reverse-complement GCF_009873245.2_mBalMus1.pri.v3_genomic.fna \
    NC_041915028.1:111914528-111915313 > bluewhale_retrogene_exons.fna

# For BlueWhaleRetrogene2B — same sequence plus flanks
samtools faidx --reverse-complement GCF_009873245.2_mBalMus1.pri.v3_genomic.fna \
    NC_041915028.1:111914228-111915613 > bluewhale_retrogene2B_reference.fna
# (111914528 - 300bp upstream; 111915313 + 300bp downstream)
# Using 300bp flanks since exact primer binding sites in blue whale unknown

# rename headers
sed -i 's/>.*/>BlueWhaleRetrogenes/' bluewhale_retrogene_exons.fna
sed -i 's/>.*/>BlueWhaleRetrogene2B/' bluewhale_retrogene2B_reference.fna