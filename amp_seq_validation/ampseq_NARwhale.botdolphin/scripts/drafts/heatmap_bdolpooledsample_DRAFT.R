## Install these packages once if they are not already installed.
# Bioconductor is used to install Rsamtools; the remaining packages come from CRAN.
# install.packages("BiocManager")
# BiocManager::install("Rsamtools")
# install.packages(c("dplyr", "tidyr", "ggplot2", "readr", "stringr"))

# Load packages for reading BAM files, manipulating data, parsing strings, and plotting.
library(Rsamtools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)

# Path to the coordinate-sorted BAM file containing aligned dolphin reads.
bam <- "C:/Users/angel/Downloads/Research project summer/ampseq_NARwhale.botdolphin/wf-amplicon_nrw_bdol_ogrefs/output/BDolphinPCNAEx1Ex6RetrogenesB/alignments/BDolphinPCNAEx1Ex6RetrogenesB.aligned.sorted.bam"

# Path to the CSV containing candidate mixed sites and their summary statistics.
mixed_sites_csv <- "C:/Users/angel/Downloads/Research project summer/ampseq_NARwhale.botdolphin/dolphin_candidate_mixed_sites.csv"

# Create a BAM index if one does not already exist, for efficient retrieval of reads from a specified genomic interval.
if (!file.exists(paste0(bam, ".bai"))) {
  indexBam(bam)
}

# Read candidate mixed sites, enforce the expected column types, and sort by position.
mixed_sites <- read_csv(mixed_sites_csv, show_col_types = FALSE) %>%
  mutate(
    pos = as.integer(pos),
    freq = as.numeric(freq),
    depth = as.integer(depth)
  ) %>%
  arrange(pos)

# Use the reference sequence/chromosome named in the first row of the site table (assumes all candidate sites belong to the same reference sequence)
seqname <- mixed_sites$seqnames[1]

# Retain well-supported candidate sites: at least 50 reads of depth and an alternate/mixed-allele frequency of at least 20%.
sites <- mixed_sites %>%
  filter(depth >= 50, freq >= 0.20) %>%
  pull(pos) %>%
  unique()

# Limit the analysis to the first 40 qualifying sites to keep the matrix and heatmap manageable. If fewer than 40 sites qualify, use all of them.
sites <- sites[seq_len(min(length(sites), 40))]

# Define the genomic interval and BAM fields to retrieve. The interval spans from the first selected site through the last selected site.
param <- ScanBamParam(
  which = GRanges(seqname, IRanges(min(sites), max(sites))),
  what = c("qname", "pos", "cigar", "seq", "mapq")
)

# Read alignments overlapping the selected interval from the BAM file.
aln <- scanBam(bam, param = param)[[1]]

# Convert the alignment data to a tibble and remove incomplete or low-quality alignments. ( mapq threshold of 20 to retain reasonably confident mappings)
reads <- tibble(
  read_id = aln$qname,
  ref_start = aln$pos,
  cigar = aln$cigar,
  seq = as.character(aln$seq),
  mapq = aln$mapq
) %>%
  filter(!is.na(read_id), !is.na(ref_start), !is.na(cigar), !is.na(seq), mapq >= 20)

# Split a CIGAR string into paired operation lengths and operation codes.
# For example, "10M2I5M" becomes lengths 10, 2, 5 and operations M, I, M.
parse_cigar <- function(cigar) {
  lens <- as.integer(str_extract_all(cigar, "\\d+")[[1]])
  ops <- str_extract_all(cigar, "[MIDNSHP=X]")[[1]]
  tibble(len = lens, op = ops)
}

# Determine which base a single aligned read carries at each candidate site.
# The function defined walks through the CIGAR operations while separately tracking positions in the reference and in the read sequence.
base_at_sites <- function(ref_start, cigar, seq, sites) {
  # Split the read sequence into individual bases.
  seq_chars <- strsplit(seq, "", fixed = TRUE)[[1]]
  
  # Initialize one missing value per candidate site. Names preserve site positions.
  out <- rep(NA_character_, length(sites))
  names(out) <- as.character(sites)
  
  # Reference coordinates are genomic. query coordinates index the read sequence.
  ref_pos <- ref_start
  query_pos <- 1L
  
  ops <- parse_cigar(cigar)
  
  # Walk through each CIGAR operation and update reference/query coordinates.
  for (i in seq_len(nrow(ops))) {
    len <- ops$len[i]
    op <- ops$op[i]
    
    if (op %in% c("M", "=", "X")) {
      # Alignment match/mismatch: both reference and query positions advance.
      ref_range <- ref_pos:(ref_pos + len - 1L)
      query_range <- query_pos:(query_pos + len - 1L)
      
      # Identify candidate sites within this aligned block and record their bases.
      hits <- match(sites, ref_range)
      has_hit <- which(!is.na(hits))
      
      if (length(has_hit) > 0) {
        out[has_hit] <- seq_chars[query_range[hits[has_hit]]]
      }
      
      ref_pos <- ref_pos + len
      query_pos <- query_pos + len
      
    } else if (op == "D") {
      # Deletion relative to the reference: only the reference position advances.
      ref_range <- ref_pos:(ref_pos + len - 1L)
      
      hits <- match(sites, ref_range)
      has_hit <- which(!is.na(hits))
      
      # Represent a deletion at a candidate site with a dash.
      if (length(has_hit) > 0) {
        out[has_hit] <- "-"
      }
      
      ref_pos <- ref_pos + len
      
    } else if (op == "N") {
      # Skipped reference region (often an intron): advance reference only.
      ref_pos <- ref_pos + len
      
    } else if (op %in% c("I", "S")) {
      # Insertion or soft clipping consumes query bases but not reference bases.
      query_pos <- query_pos + len
      
    } else if (op %in% c("H", "P")) {
      # Hard clipping and padding consume neither stored query nor reference bases.
      next
    }
  }
  
  return(out)
}

# Apply the base-extraction function to every retained read.
allele_matrix <- lapply(seq_len(nrow(reads)), function(i) {
  base_at_sites(
    ref_start = reads$ref_start[i],
    cigar = reads$cigar[i],
    seq = reads$seq[i],
    sites = sites
  )
})

# Combine results into a read-by-site table and add read identifiers as column 1.
allele_df <- as_tibble(do.call(rbind, allele_matrix), .name_repair = "minimal") %>%
  mutate(read_id = reads$read_id, .before = 1)

# Require each read to have calls at a minimum of five sites or 50% of all selected sites, whichever requirement is larger.
min_sites_per_read <- max(5, ceiling(length(sites) * 0.50))

# Count non-missing site calls per read and discard reads with insufficient coverage.
allele_df <- allele_df %>%
  mutate(n_called = rowSums(!is.na(across(-read_id)))) %>%
  filter(n_called >= min_sites_per_read)

# Save the filtered read-by-site allele matrix for downstream inspection.
write_csv(allele_df, "dolphin_read_by_site_alleles.csv")

# Convert the wide allele matrix to long format: one row per read/site combination.
long <- allele_df %>%
  select(-n_called) %>%
  pivot_longer(-read_id, names_to = "pos", values_to = "allele") %>%
  mutate(pos = as.integer(pos))

# At each site, rank observed alleles by read count. The most frequent allele is labeled major and the second-most frequent allele is labeled minor.
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

# Attach the major/minor labels to every observation and classify each call.
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

# Order reads primarily by the fraction of sites carrying minor alleles, then by the number of called sites, to visually group similar allele patterns.
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

# Convert read IDs and positions to factors to ensure that ggplot uses the intended order. Reversing read_order places reads with higher minor-allele fractions near the top
plot_df <- long_labeled %>%
  mutate(
    read_id = factor(read_id, levels = rev(read_order)),
    pos = factor(pos, levels = sort(unique(pos)))
  )

# Draw a heatmap:  each tile represents one read at one candidate site
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

# Save the heatmap as a high-resolution PNG.
ggsave("dolphin_read_phasing_heatmap.png", p, width = 11, height = 8, dpi = 300)

# Save the inferred major and minor allele at each candidate site.
write_csv(site_labels, "dolphin_phasing_site_major_minor_labels.csv")
