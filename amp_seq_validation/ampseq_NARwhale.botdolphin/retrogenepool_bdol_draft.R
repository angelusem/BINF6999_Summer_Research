## Purpose: Detect candidate mixed allele sites in the dolphin retrogene BAM.These sites can indicate two co-amplified retrogene copies if the same alternate alleles recur at substantial frequency across many reads.

## Install once if needed:
## install.packages("BiocManager")
## BiocManager::install(c("Rsamtools", "GenomicRanges"))
## install.packages(c("dplyr", "ggplot2", "readr"))

library(Rsamtools)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(readr)

dolphinRpooled_bam <- "C:/Users/angel/Downloads/Research project summer/ampseq_NARwhale.botdolphin/wf-amplicon_nrw_bdol_ogrefs/output/BDolphinPCNAEx1Ex6RetrogenesB/alignments/BDolphinPCNAEx1Ex6RetrogenesB.aligned.sorted.bam"
dolphinRpooled_bai <- paste0(dolphinRpooled_bam, ".bai")
#
if (!file.exists(bam)) {
  stop("BAM file not found: ", bam)
}

if (!file.exists(bai)) {
  message("BAM index not found. Attempting to create it...")
  indexBam(bam)
}
#
# inspect reference names and lengths stored in the BAM.
dolphinRpooled_bam_header <- scanBamHeader(dolphinRpooled_bam)[[1]]
targets <- dolphinRpooled_bam_header$targets
print(targets)

#Pile up bases at every covered position.
pile <- pileup(
  dolphinRpooled_bam,
  pileupParam = PileupParam(
    distinguish_strands = FALSE,
    distinguish_nucleotides = TRUE,
    ignore_query_Ns = TRUE,
    min_base_quality = 10,
    min_mapq = 20
  )
)

allele_freqs <- pile %>%
  group_by(seqnames, pos) %>%
  mutate(
    depth = sum(count),
    freq = count / depth
  ) %>%
  ungroup() %>%
  arrange(seqnames, pos, desc(freq))

# Candidate mixed sites: second-most-common allele is substantial.
mixed_sites <- allele_freqs %>%
  group_by(seqnames, pos) %>%
  mutate(rank = row_number()) %>%
  filter(depth >= 50, rank == 2, freq >= 0.20) %>%
  ungroup() %>%
  arrange(seqnames, pos)

mixed_sites

write_csv(allele_freqs, "dolphin_allele_frequencies.csv")
write_csv(mixed_sites, "dolphin_candidate_mixed_sites.csv")

message("Candidate mixed sites: ", nrow(mixed_sites))
print(mixed_sites)

if (nrow(mixed_sites) > 0) {
  p <- ggplot(mixed_sites, aes(x = pos, y = freq)) +
    geom_point(size = 2) +
    facet_wrap(~ seqnames, scales = "free_x") +
    theme_bw() +
    labs(
      title = "Candidate mixed allele sites in dolphin retrogene sample",
      subtitle = "Second-most-common base frequency; depth >= 50, minor allele >= 20%",
      x = "Reference position",
      y = "Minor allele frequency"
    )
  
  ggsave("dolphin_candidate_mixed_sites.png", p, width = 8, height = 4.5, dpi = 300)
}

