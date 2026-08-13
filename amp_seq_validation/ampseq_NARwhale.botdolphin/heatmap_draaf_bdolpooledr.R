## Install once if needed:
## install.packages("BiocManager")
## BiocManager::install("Rsamtools")
## install.packages(c("dplyr", "tidyr", "ggplot2", "readr", "stringr"))

library(Rsamtools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)

bam <- "C:/Users/angel/Downloads/Research project summer/ampseq_NARwhale.botdolphin/wf-amplicon_nrw_bdol_ogrefs/output/BDolphinPCNAEx1Ex6RetrogenesB/alignments/BDolphinPCNAEx1Ex6RetrogenesB.aligned.sorted.bam"

mixed_sites_csv <- "C:/Users/angel/Downloads/Research project summer/ampseq_NARwhale.botdolphin/dolphin_candidate_mixed_sites.csv"

if (!file.exists(paste0(bam, ".bai"))) {
  indexBam(bam)
}

mixed_sites <- read_csv(mixed_sites_csv, show_col_types = FALSE) %>%
  mutate(
    pos = as.integer(pos),
    freq = as.numeric(freq),
    depth = as.integer(depth)
  ) %>%
  arrange(pos)

seqname <- mixed_sites$seqnames[1]

sites <- mixed_sites %>%
  filter(depth >= 50, freq >= 0.20) %>%
  pull(pos) %>%
  unique()

sites <- sites[seq_len(min(length(sites), 40))]

param <- ScanBamParam(
  which = GRanges(seqname, IRanges(min(sites), max(sites))),
  what = c("qname", "pos", "cigar", "seq", "mapq")
)

aln <- scanBam(bam, param = param)[[1]]

reads <- tibble(
  read_id = aln$qname,
  ref_start = aln$pos,
  cigar = aln$cigar,
  seq = as.character(aln$seq),
  mapq = aln$mapq
) %>%
  filter(!is.na(read_id), !is.na(ref_start), !is.na(cigar), !is.na(seq), mapq >= 20)

parse_cigar <- function(cigar) {
  lens <- as.integer(str_extract_all(cigar, "\\d+")[[1]])
  ops <- str_extract_all(cigar, "[MIDNSHP=X]")[[1]]
  tibble(len = lens, op = ops)
}

base_at_sites <- function(ref_start, cigar, seq, sites) {
  seq_chars <- strsplit(seq, "", fixed = TRUE)[[1]]
  out <- rep(NA_character_, length(sites))
  names(out) <- as.character(sites)
  
  ref_pos <- ref_start
  query_pos <- 1L
  
  ops <- parse_cigar(cigar)
  
  for (i in seq_len(nrow(ops))) {
    len <- ops$len[i]
    op <- ops$op[i]
    
    if (op %in% c("M", "=", "X")) {
      ref_range <- ref_pos:(ref_pos + len - 1L)
      query_range <- query_pos:(query_pos + len - 1L)
      
      hits <- match(sites, ref_range)
      has_hit <- which(!is.na(hits))
      
      if (length(has_hit) > 0) {
        out[has_hit] <- seq_chars[query_range[hits[has_hit]]]
      }
      
      ref_pos <- ref_pos + len
      query_pos <- query_pos + len
      
    } else if (op == "D") {
      ref_range <- ref_pos:(ref_pos + len - 1L)
      
      hits <- match(sites, ref_range)
      has_hit <- which(!is.na(hits))
      
      if (length(has_hit) > 0) {
        out[has_hit] <- "-"
      }
      
      ref_pos <- ref_pos + len
      
    } else if (op == "N") {
      ref_pos <- ref_pos + len
      
    } else if (op %in% c("I", "S")) {
      query_pos <- query_pos + len
      
    } else if (op %in% c("H", "P")) {
      next
    }
  }
  
  out
}

allele_matrix <- lapply(seq_len(nrow(reads)), function(i) {
  base_at_sites(
    ref_start = reads$ref_start[i],
    cigar = reads$cigar[i],
    seq = reads$seq[i],
    sites = sites
  )
})

allele_df <- as_tibble(do.call(rbind, allele_matrix), .name_repair = "minimal") %>%
  mutate(read_id = reads$read_id, .before = 1)

min_sites_per_read <- max(5, ceiling(length(sites) * 0.50))

allele_df <- allele_df %>%
  mutate(n_called = rowSums(!is.na(across(-read_id)))) %>%
  filter(n_called >= min_sites_per_read)

write_csv(allele_df, "dolphin_read_by_site_alleles.csv")

long <- allele_df %>%
  select(-n_called) %>%
  pivot_longer(-read_id, names_to = "pos", values_to = "allele") %>%
  mutate(pos = as.integer(pos))

site_labels <- long %>%
  filter(!is.na(allele)) %>%
  count(pos, allele, name = "n") %>%
  group_by(pos) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  summarise(
    major = allele[rank == 1][1],
    minor = allele[rank == 2][1],
    .groups = "drop"
  )

long_labeled <- long %>%
  left_join(site_labels, by = "pos") %>%
  mutate(
    allele_class = case_when(
      is.na(allele) ~ "not covered",
      allele == major ~ "major",
      allele == minor ~ "minor",
      TRUE ~ "other"
    )
  )

read_order <- long_labeled %>%
  group_by(read_id) %>%
  summarise(
    minor_fraction = mean(allele_class == "minor", na.rm = TRUE),
    major_fraction = mean(allele_class == "major", na.rm = TRUE),
    called = sum(allele_class != "not covered"),
    .groups = "drop"
  ) %>%
  arrange(desc(minor_fraction), desc(called)) %>%
  pull(read_id)

plot_df <- long_labeled %>%
  mutate(
    read_id = factor(read_id, levels = rev(read_order)),
    pos = factor(pos, levels = sort(unique(pos)))
  )

p <- ggplot(plot_df, aes(x = pos, y = read_id, fill = allele_class)) +
  geom_tile(color = NA) +
  scale_fill_manual(
    values = c(
      "major" = "#f2f2f2",
      "minor" = "#2b6cb0",
      "other" = "#d95f02",
      "not covered" = "#222222"
    )
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Read-level allele pattern across dolphin mixed sites",
    subtitle = "Rows are reads ordered by minor-allele fraction; blue blocks indicate reads carrying the site-specific minor allele",
    x = "Candidate mixed site position",
    y = "Reads",
    fill = "Allele class"
  )

ggsave("dolphin_read_phasing_heatmap.png", p, width = 11, height = 8, dpi = 300)

write_csv(site_labels, "dolphin_phasing_site_major_minor_labels.csv")