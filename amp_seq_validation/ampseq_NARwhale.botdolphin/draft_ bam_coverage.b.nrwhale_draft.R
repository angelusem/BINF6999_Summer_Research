# Purpose:
  # Validate primer removal and expected amplicon boundaries for the NARW A and B consensus sequences.

# For the B amplicon, quantify primary-read coverage, count reads spanning both PCNA-like/flank junctions, examine alignment, endpoints, and verify that primer sequences are absent from the final de novo consensus.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "Biostrings",
  "Rsamtools",
  "GenomicAlignments",
  "GenomicRanges"
))

install.packages("tidyverse")

library(Biostrings)
library(tidyverse)
library(Rsamtools)
library(GenomicAlignments)
library(GenomicRanges)
library(ggplot2)

#------ 1. Confirm primer positions ----------
# loading in the sequences
fasta_file <- paste0(
  "C:/Users/angel/Downloads/Research project summer/",
  "ampseq_NARwhale.botdolphin/",
  "compilation_consensus_comparison_NRwhale"
)

seqs <- readDNAStringSet(fasta_file)

A_reference <- seqs[["NARW01_A_variant_excelref"]]
B_reference <- seqs[["NARW01_B_variant_excelref"]]
PCNA_like   <- seqs[["REF_PCNA_like"]]
# define my primers from known data
primers <- DNAStringSet(c(
  A_primer1 = "ATGTCTCCTTGGTGCAGCTC",
  A_primer2 = "CATCCTCGATCTTGGGAGCC",
  B_primer1 = "GACCAGATTTGACTTTGGACTTTA",
  B_primer2 = "GGGGTAAGAAGAGACTGCTTG"
))
# search sequences for primers in both orientations
find_both_orientations <- function(primer, reference) {
  forward <- matchPattern(primer, reference)
  reverse <- matchPattern(reverseComplement(primer), reference)
  
  bind_rows(
    tibble(
      orientation = "as_given",
      start = start(forward),
      end = end(forward)
    ),
    tibble(
      orientation = "reverse_complement",
      start = start(reverse),
      end = end(reverse)
    )
  )
}

find_both_orientations(primers[["A_primer1"]], A_reference)
find_both_orientations(primers[["A_primer2"]], A_reference)

find_both_orientations(primers[["B_primer1"]], B_reference)
find_both_orientations(primers[["B_primer2"]], B_reference)

# locate my pcna-like core
core_match <- matchPattern(
  reverseComplement(PCNA_like),
  B_reference
)

core_match
# as expected, 786 bp
# ------ 2. comparing observed consensus with expected intervals ----
A_denovo <- seqs[["NARW01_A_denovo"]]
B_denovo <- seqs[["NARW01_B_denovo"]]

A_expected <- subseq(
  A_reference,
  start = 15,
  end = 656
)

B_expected <- subseq(
  B_reference,
  start = 4,
  end = 1165
)
# checking
A_denovo == A_expected
# since  b denovo sequnece is in opposite orientation
B_denovo_oriented <- reverseComplement(B_denovo)

matchPattern(B_denovo_oriented, B_reference)
# building a summary table
primer_summary <- tibble(
  design = c("A", "B"),
  reference_length = c(
    length(A_reference),
    length(B_reference)
  ),
  expected_start = c(15, 4),
  expected_end = c(656, 1165),
  expected_amplicon_length = c(642, 1162),
  denovo_length = c(
    length(A_denovo),
    length(B_denovo)
  )
)

primer_summary

# ----- 3. calculate BAM coverage -----


bam_file <- "C:/Users/angel/Downloads/Research project summer/ampseq_NARwhale.botdolphin/originalrun_notref_amplicon_variants/amplicon_assembly_report/NARightWhaleRetrogeneB/alignments/NARightWhaleRetrogeneB.aligned.sorted.bam"

bam_header <- scanBamHeader(bam_file)
names(bam_header[[1]]$targets)
target_name <- "NARightWhaleRetrogeneB" 

# excluding secondary and supplementary alignments

bam_param <- ScanBamParam(
  which = GRanges(
    seqnames = target_name,
    ranges = IRanges(1, 1186)
  ),
  flag = scanBamFlag(
    isUnmappedQuery = FALSE,
    isSecondaryAlignment = FALSE,
    isSupplementaryAlignment = FALSE
  ),
  what = c("mapq", "cigar", "qname")
)

alignments <- readGAlignments(
  bam_file,
  param = bam_param,
  use.names = TRUE
)

alignments

# retain just the mappings that are reasonably confdent: alignments <- alignments[mcols(alignments)$mapq >= 20]

# coverage calculation
depth_rle <- coverage(alignments)[[target_name]]

depth <- as.integer(depth_rle)[1:1186]

coverage_table <- tibble(
  position = 1:1186,
  depth = depth,
  region = case_when(
    position < 4 ~ "external_left_context",
    position <= 27 ~ "B_primer1",
    position < 201 ~ "left_flank",
    position <= 986 ~ "PCNA_like_core",
    position < 1145 ~ "right_flank",
    position <= 1165 ~ "B_primer2_binding_site",
    TRUE ~ "external_right_context"
  )
)

# summarized by region
coverage_summary <- coverage_table |>
  group_by(region) |>
  summarise(
    start = min(position),
    end = max(position),
    number_positions = n(),
    positions_covered = sum(depth > 0),
    fraction_covered = mean(depth > 0),
    minimum_depth = min(depth),
    median_depth = median(depth),
    mean_depth = mean(depth),
    maximum_depth = max(depth),
    .groups = "drop"
  )

coverage_summary

# plotting it


ggplot(coverage_table, aes(position, depth)) +
  geom_area(fill = "grey45") +
  geom_vline(
    xintercept = c(4, 27, 201, 986, 1145, 1165),
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = 201,
    xmax = 986,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.08,
    fill = "blue"
  ) +
  labs(
    title = "NARW01 primer-B read depth",
    x = "B-reference position",
    y = "Read depth"
  ) +
  theme_classic()

# ------ 4. counting junction spanning reads -----
alignment_ranges <- granges(alignments)

left_junction <- 200
right_junction <- 986
minimum_anchor <- 20

left_spanning <- (
  start(alignment_ranges) <= left_junction - minimum_anchor
) & (
  end(alignment_ranges) >= left_junction + minimum_anchor
)

right_spanning <- (
  start(alignment_ranges) <= right_junction - minimum_anchor
) & (
  end(alignment_ranges) >= right_junction + minimum_anchor
)

junction_counts <- tibble(
  junction = c(
    "left_flank_to_PCNA_like",
    "PCNA_like_to_right_flank"
  ),
  coordinate = c(
    left_junction,
    right_junction
  ),
  reads_spanning_with_20bp_anchors = c(
    sum(left_spanning),
    sum(right_spanning)
  )
)

junction_counts
# coutning reads spanning the entire pcna core plus the "anchors"
full_core_spanning <- (
  start(alignment_ranges) <= 201 - 20
) & (
  end(alignment_ranges) >= 986 + 20
)

sum(full_core_spanning)
# ---5. plotting alignment start and end positions-----
read_boundaries <- tibble(
  read_name = names(alignments),
  alignment_start = start(alignment_ranges),
  alignment_end = end(alignment_ranges),
  mapq = mcols(alignments)$mapq
)

ggplot(read_boundaries, aes(alignment_start)) +
  geom_histogram(binwidth = 5) +
  geom_vline(xintercept = 4, color = "red") +
  labs(
    title = "B-read alignment start positions",
    x = "Reference position"
  ) +
  theme_classic()

ggplot(read_boundaries, aes(alignment_end)) +
  geom_histogram(binwidth = 5) +
  geom_vline(xintercept = 1165, color = "red") +
  labs(
    title = "B-read alignment end positions",
    x = "Reference position"
  ) +
  theme_classic()

# confirmation that my primers are not present in the consensus sequence.

matchPattern(DNAString("GACCAGATTTGACTTTGGACTTTA"), B_denovo)

matchPattern(
  reverseComplement(DNAString("GGGGTAAGAAGAGACTGCTTG")),
  B_denovo
)

# ----- end -----