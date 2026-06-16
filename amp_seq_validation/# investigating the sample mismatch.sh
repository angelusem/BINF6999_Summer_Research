# investigating the sample mismatch
# blasting retrogene2b consensus sequences for both beluga and blue whale 
# create folder for this 
# make blast dbs for each reference genome 
makeblastdb \
  -in /home/angelusm/binf6999_summer_research/amp_seq_validation/GCF_002288925.2_ASM228892v3_genomic.fna \
  -dbtype nucl \
  -out beluga_genome_db

makeblastdb \
  -in /home/angelusm/binf6999_summer_research/amp_seq_validation/GCF_009873245.2_mBalMus1.pri.v3_genomic.fna \
  -dbtype nucl \
  -out bluewhale_genome_db
# copying in the 2 consensus sequence files
cp ../B12BelugaRetrogene2B/consensus/medaka.consensus.fasta B12BelugaRetrogene2B.consensus.fasta
cp ../BlueWhaleRetrogene2B/consensus/medaka.consensus.fasta BlueWhaleRetrogene2B.consensus.fasta
# 4 consensus blasts run
blastn \
  -query B12BelugaRetrogene2B.consensus.fasta \
  -db beluga_genome_db \
  -out B12BelugaRetrogene2B_vs_beluga_genome.tsv \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -max_target_seqs 10

blastn \
  -query B12BelugaRetrogene2B.consensus.fasta \
  -db bluewhale_genome_db \
  -out B12BelugaRetrogene2B_vs_bluewhale_genome.tsv \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -max_target_seqs 10

blastn \
  -query BlueWhaleRetrogene2B.consensus.fasta \
  -db beluga_genome_db \
  -out BlueWhaleRetrogene2B_vs_beluga_genome.tsv \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -max_target_seqs 10

blastn \
  -query BlueWhaleRetrogene2B.consensus.fasta \
  -db bluewhale_genome_db \
  -out BlueWhaleRetrogene2B_vs_bluewhale_genome.tsv \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -max_target_seqs 10

# view top hits
for f in *.tsv; do
  echo "### $f"
  sort -k12,12nr "$f" | head -5
done
#B12BelugaRetrogene2B.consensus best against blue whale genome?
# BlueWhaleRetrogene2B.consensus best against beluga genome?
# results:
