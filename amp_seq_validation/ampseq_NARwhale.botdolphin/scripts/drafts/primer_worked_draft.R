#to be confident that primer B amplified the intended retrocopy locus—not merely some PCNA-related sequence—you need three layers of evidence:
#The B consensus contains the expected PCNA-like core.
#the sequence on both sides matches the unique flanks of the intended genomic locus.
#Individual reads physically span each flank-to-core junction.

#load sequences
library(Biostrings)

fasta_file <- paste0(
  "C:/Users/angel/Downloads/Research project summer/",
  "ampseq_NARwhale.botdolphin/",
  "compilation_consensus_comparison_NRwhale"
)

seqs <- readDNAStringSet(fasta_file)

names(seqs)
width(seqs)
# extract relevant records
B_denovo <- seqs[["NARW01_B_denovo"]]
A_denovo <- seqs[["NARW01_A_denovo"]]

pcna_like <- seqs[["REF_PCNA_like"]]
pcna_cds  <- seqs[["REF_PCNA_CDS"]]

# check their lengths
c(
  B_denovo = length(B_denovo),
  A_denovo = length(A_denovo),
  PCNA_like = length(pcna_like),
  PCNA_CDS = length(pcna_cds)
)
# orient B sequences consistently
# searching the fowards orientation 
forward_match <- matchPattern(
  pcna_like,
  B_denovo
)
# and the reverse orientation
reverse_match <- matchPattern(
  pcna_like,
  reverseComplement(B_denovo)
)

reverse_match
forward_match

# generally speaking, a orientation rule could be 
if (length(forward_match) > 0) {
  B_oriented <- B_denovo
} else if (length(reverse_match) > 0) {
  B_oriented <- reverseComplement(B_denovo)
} else {
  stop("The complete PCNA-like reference is not an exact match in B.")
}

# locate the pcna-like core
core_match <- matchPattern(
  pcna_like,
  B_oriented
)

core_match
# extract observed b bases
core_start <- start(core_match)[1]
core_end <- end(core_match)[1]

B_core <- subseq(
  B_oriented,
  start = core_start,
  end = core_end
)

B_core
length(B_core)
# confirm that extracted bases are identical 
identical(
  as.character(B_core),
  as.character(pcna_like)
)


#extract the B flanks
B_left_flank <- subseq(
  B_oriented,
  start = 1,
  end = core_start - 1
)

B_right_flank <- subseq(
  B_oriented,
  start = core_end + 1,
  end = length(B_oriented)
)

c(
  left_flank = length(B_left_flank),
  core = length(B_core),
  right_flank = length(B_right_flank)
)
# not sure if useful, but just in case: save them separately
# important to determine which genomic locus was amplified
B_components <- DNAStringSet(c(
  NARW_B_left_flank = B_left_flank,
  NARW_B_PCNA_like_core = B_core,
  NARW_B_right_flank = B_right_flank
))

writeXStringSet(
  B_components,
  "NARW_B_components.fasta"
)

# verifying that A and B independently agree
# wkt that A sequence is in the opposite orientation and also contains the 20 bp primers
# orienting it
A_oriented <- reverseComplement(A_denovo)
# removing 20 bp from each end
A_primer_free <- subseq(
  A_oriented,
  start = 21,
  end = length(A_oriented) - 20
)

length(A_primer_free)
# can we find A inside B
A_in_B <- matchPattern(
  A_primer_free,
  B_oriented
)

A_in_B

A_in_B_core <- matchPattern(
  A_primer_free,
  B_core
)

A_in_B_core
# just serves as additional technical confirmation 
# confirming that the sequence differs from parent pCNA 
like_bases <- strsplit(
  as.character(B_core),
  ""
)[[1]]

parent_bases <- strsplit(
  as.character(pcna_cds),
  ""
)[[1]]

diagnostic_positions <- which(
  like_bases != parent_bases
)

length(diagnostic_positions)
# create a diagnostic table
diagnostic_table <- data.frame(
  position = diagnostic_positions,
  B_core_base = like_bases[diagnostic_positions],
  parent_PCNA_base = parent_bases[diagnostic_positions]
)

diagnostic_table

#compare the complete B sequence with the intended genomic locus
B_genomic_reference <- readDNAStringSet(
  "path/to/rightwhale_PCNA_like_with_flanks.fasta"
)[[1]]

# try both orientations
match_forward <- matchPattern(
  B_oriented,
  B_genomic_reference
)

match_reverse <- matchPattern(
  B_oriented,
  reverseComplement(B_genomic_reference)
)

match_forward
match_reverse
# 8. check whether the B flanks are unique in the genome 
# align the complete B consensus to right whale genome? blast or minimap2
minimap2 -x asm20 right_whale_genome.fna NARW_B_denovo.fasta > NARW_B_vs_genome.paf
# Look for:
#One near-full-length primary alignment
#Approximately 992 bp aligned
#High identity
#Both flanks included
#No equally good alignment at parental PCNA or another PCNA-like locus

#Align the flanks separately as well:
minimap2 -x asm20 right_whale_genome.fna NARW_B_components.fasta \
> NARW_B_components_vs_genome.paf
#The ideal result is:
#Left flank uniquely maps adjacent to the annotated PCNA-like locus.
#Right flank uniquely maps adjacent to the same locus.
#The core maps between them.
#The complete B sequence does not align equally well to parental PCNA.
#If the flanks map to multiple places, use longer flanks or inspect mapping quality and secondary alignments.

#Confirm the primer sites in the intended reference
# B primers are 
B_forward_primer <- DNAString(
  "GACCAGATTTGACTTTGGACTTTA"
)

B_reverse_primer <- DNAString(
  "GGGGTAAGAAGAGACTGCTTG"
)
# reverse-binding primer sites in the fwd-orientated ref is
B_reverse_binding <- reverseComplement(
  B_reverse_primer
)

B_reverse_binding

# search for them in the full reference used for variant calling mode:
find_forward <- matchPattern(
  B_forward_primer,
  B_genomic_reference
)

find_reverse <- matchPattern(
  B_reverse_binding,
  B_genomic_reference
)

find_forward
find_reverse
# They should:
#Occur in the correct inward-facing orientation
#Surround the PCNA-like core
#Produce the expected PCR product
#Not produce another similarly sized high-quality product elsewhere in the genome
#Use NCBI Primer-BLAST or another in-silico PCR tool against the right-whale genome for the final specificity check.

# 10. confirm the locus with individual readss: using alignment file generated to verfify that my individual reds connect the flanks to the core -- IGV
# for a junction spanning read count in R
library(Rsamtools)
library(GenomicAlignments)
library(GenomicRanges)

bam_file <- "path/to/B_reads_mapped_to_B_denovo.bam"

alignments <- readGAlignments(
  bam_file,
  use.names = TRUE,
  param = ScanBamParam(
    flag = scanBamFlag(
      isUnmappedQuery = FALSE,
      isSecondaryAlignment = FALSE,
      isSupplementaryAlignment = FALSE
    ),
    what = c("mapq", "cigar")
  )
)

alignment_ranges <- granges(alignments)

minimum_anchor <- 20

left_spanning <- (
  start(alignment_ranges) <= 115 - minimum_anchor
) & (
  end(alignment_ranges) >= 116 + minimum_anchor
)

right_spanning <- (
  start(alignment_ranges) <= 901 - minimum_anchor
) & (
  end(alignment_ranges) >= 902 + minimum_anchor
)

full_locus_spanning <- (
  start(alignment_ranges) <= 115 - minimum_anchor
) & (
  end(alignment_ranges) >= 902 + minimum_anchor
)

c(
  left_junction_spanning_reads = sum(left_spanning),
  right_junction_spanning_reads = sum(right_spanning),
  complete_core_plus_flanks = sum(full_locus_spanning)
) 

##
# saving by b primer core
library(Biostrings)

# B_core was created previously with subseq()
B_core_fasta <- DNAStringSet(B_core)

names(B_core_fasta) <- "NARW01_B_PCNA_like_core"

output_file <- paste0(
  "C:/Users/angel/Downloads/Research project summer/",
  "ampseq_NARwhale.botdolphin/",
  "NARW01_B_PCNA_like_core.fasta"
)

writeXStringSet(
  B_core_fasta,
  filepath = output_file,
  format = "fasta"
)

cat("FASTA written to:\n", output_file, "\n")
#
file.exists(output_file)

check <- readDNAStringSet(output_file)

names(check)
width(check)
check
# check if sequence is unchanged
identical(
  as.character(check[[1]]),
  as.character(B_core)
)
#
B_parts <- DNAStringSet(c(
  NARW01_B_full_denovo = B_oriented,
  NARW01_B_left_flank = B_left_flank,
  NARW01_B_PCNA_like_core = B_core,
  NARW01_B_right_flank = B_right_flank
))

writeXStringSet(
  B_parts,
  filepath = paste0(
    "C:/Users/angel/Downloads/Research project summer/",
    "ampseq_NARwhale.botdolphin/",
    "NARW01_B_full_and_components.fasta"
  ),
  format = "fasta"
)