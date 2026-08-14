#!/usr/bin/env Rscript

# PCNA/PCNA-like consensus comparison
#
# Required packages:
# install.packages(c("ape", "ggplot2", "readr", "tidyr", "dplyr", "purrr"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("DECIPHER"))

required <- c(
  "Biostrings", "DECIPHER", "ape", "ggplot2",
  "readr", "tidyr", "dplyr", "purrr"
)

missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0) {
  stop(
    "Install the missing packages before running: ",
    paste(missing, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(Biostrings)
  library(DECIPHER)
  library(ape)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(dplyr)
  library(purrr)
})

input_fasta <- paste0(
  "C:/Users/angel/Downloads/Research project summer/",
  "ampseq_NARwhale.botdolphin/",
  "PCNA_comparison_nrwhale_bdol.fasta"
)

result_dir <- paste0(
  "C:/Users/angel/Downloads/Research project summer/",
  "ampseq_NARwhale.botdolphin/results"
)

stopifnot(file.exists(input_fasta))

dir.create(
  result_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

sequences <- readDNAStringSet(input_fasta)

names(sequences)
width(sequences)
# The full-core analysis excludes A because A covers only a 602-bp subregion.
full_core_names <- setdiff(names(sequences), "NARW_A_primer_free")
full_core <- sequences[full_core_names]
full_alignment <- AlignSeqs(full_core, processors = NULL, verbose = FALSE)

# The shared-region analysis includes A. Align all records and retain only
# alignment columns in which A contains a nucleotide.
shared_alignment_all <- AlignSeqs(
  sequences,
  processors = NULL,
  verbose = FALSE
)
shared_matrix_all <- as.matrix(shared_alignment_all)
a_row <- shared_matrix_all["NARW_A_primer_free", ]
a_columns <- a_row %in% c("A", "C", "G", "T")
shared_matrix <- shared_matrix_all[, a_columns, drop = FALSE]
shared_alignment <- DNAStringSet(
  apply(shared_matrix, 1, paste0, collapse = "")
)

writeXStringSet(
  full_alignment,
  file.path(result_dir, "PCNA_full_core_alignment.fasta")
)
writeXStringSet(
  shared_alignment,
  file.path(result_dir, "PCNA_A_shared_region_alignment.fasta")
)

pairwise_stats <- function(aligned_sequences, dataset_name) {
  x <- as.matrix(aligned_sequences)
  sequence_names <- rownames(x)
  pairs <- combn(sequence_names, 2, simplify = FALSE)
  
  map_dfr(pairs, function(pair) {
    a <- x[pair[1], ]
    b <- x[pair[2], ]
    valid_a <- a %in% c("A", "C", "G", "T")
    valid_b <- b %in% c("A", "C", "G", "T")
    comparable <- valid_a & valid_b
    gap_columns <- xor(a == "-", b == "-")
    matches <- sum(a[comparable] == b[comparable])
    substitutions <- sum(a[comparable] != b[comparable])
    compared <- sum(comparable)
    
    tibble(
      dataset = dataset_name,
      sequence_1 = pair[1],
      sequence_2 = pair[2],
      alignment_columns = length(a),
      compared_nucleotides = compared,
      identical_nucleotides = matches,
      substitutions = substitutions,
      gap_columns = sum(gap_columns),
      percent_identity = ifelse(
        compared > 0,
        100 * matches / compared,
        NA_real_
      )
    )
  })
}

pairwise_results <- bind_rows(
  pairwise_stats(full_alignment, "full_core"),
  pairwise_stats(shared_alignment, "A_shared_region")
)

pairwise_results

write_csv(
  pairwise_results,
  file.path(result_dir, "PCNA_pairwise_comparisons.csv")
)

# Create a symmetrical identity heatmap for the full-core dataset.
# Use one explicitly constructed record per matrix cell. This prevents
# overlapping diagonal labels if the script is sourced repeatedly or if
# duplicated pair records enter the plotting table.
full_pairs <- pairwise_results |>
  filter(dataset == "full_core") |>
  distinct(sequence_1, sequence_2, .keep_all = TRUE)

identity_long <- bind_rows(
  full_pairs |>
    select(sequence_1, sequence_2, percent_identity),
  full_pairs |>
    transmute(
      sequence_1 = sequence_2,
      sequence_2 = sequence_1,
      percent_identity
    ),
  tibble(
    sequence_1 = names(full_alignment),
    sequence_2 = names(full_alignment),
    percent_identity = 100
  )
) |>
  distinct(sequence_1, sequence_2, .keep_all = TRUE)

display_labels <- c(
  NARW_B_PCNA_like_core = "NARW B core",
  NARW_NCBI_PCNA_like = "NARW NCBI PCNA-like",
  NARW_parent_PCNA_CDS = "NARW parent PCNA CDS",
  Bottlenose_PCNA_pseudogene1_core = "Bottlenose pseudogene 1",
  Bottlenose_PCNA_pseudogene2_core = "Bottlenose pseudogene 2"
)

sequence_order <- names(full_alignment)

identity_long <- identity_long |>
  mutate(
    column_position = match(sequence_1, sequence_order),
    row_position    = match(sequence_2, sequence_order)
  ) |>
  # Keep the lower triangle, excluding the 100% diagonal
 filter(row_position > column_position) |>
  mutate(
    sequence_1 = factor(
      sequence_1,
      levels = sequence_order,
      labels = unname(display_labels[sequence_order])
    ),
    sequence_2 = factor(
      sequence_2,
      levels = rev(sequence_order),
      labels = unname(display_labels[rev(sequence_order)])
    )
  )
identity_plot <- ggplot(
  identity_long,
  aes(sequence_1, sequence_2, fill = percent_identity)
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = sprintf("%.1f", percent_identity)),
    size = 3
  ) +
  scale_fill_viridis_c() +
  coord_equal() +
  labs(
    title = "PCNA-derived sequence identity",
    subtitle = "Identity calculated across comparable A/C/G/T positions; gap columns excluded",
    x = NULL,
    y = NULL,
    fill = "Identity (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank(),
    plot.margin = margin(15, 30, 20, 15)
  )

ggsave(
  file.path(result_dir, "PCNA_identity_heatmap.png"),
  identity_plot,
  width = 9.5,
  height = 8,
  dpi = 300
)

# Exploratory unrooted neighbor-joining tree. The NARW B core and the NCBI
# PCNA-like sequence are identical, so collapse them into one tip rather than
# plotting two labels at the same zero-distance coordinate.
tree_alignment <- full_alignment[
  names(full_alignment) != "NARW_NCBI_PCNA_like"
]
names(tree_alignment)[
  names(tree_alignment) == "NARW_B_PCNA_like_core"
] <- "NARW_sample_and_NCBI_PCNA_like_identical"

tree_tip_labels <- c(
  NARW_sample_and_NCBI_PCNA_like_identical =
    "NARW sample + NCBI PCNA-like (identical)",
  NARW_parent_PCNA_CDS = "NARW parent PCNA CDS",
  Bottlenose_PCNA_pseudogene1_core = "Bottlenose pseudogene 1",
  Bottlenose_PCNA_pseudogene2_core = "Bottlenose pseudogene 2"
)

full_dnabin <- as.DNAbin(as.matrix(tree_alignment))
full_distance <- dist.dna(
  full_dnabin,
  model = "raw",
  pairwise.deletion = TRUE
)
tree <- nj(full_distance)
tree$tip.label <- unname(tree_tip_labels[tree$tip.label])

write.tree(
  tree,
  file = file.path(result_dir, "PCNA_full_core_NJ_tree.newick")
)

png(
  file.path(result_dir, "PCNA_full_core_NJ_tree.png"),
  width = 2700,
  height = 2100,
  res = 300
)
par(mar = c(2, 2, 4, 2))
plot(
  tree,
  type = "unrooted",
  cex = 0.9,
  label.offset = 0.0015,
  main = "Exploratory unrooted NJ tree of PCNA-derived sequences",
  no.margin = FALSE
)
add.scale.bar()
dev.off()

# Detailed differences between the NARW B core and each reference.
difference_table <- function(aligned_sequences, focal_name) {
  x <- as.matrix(aligned_sequences)
  focal <- x[focal_name, ]
  map_dfr(setdiff(rownames(x), focal_name), function(other_name) {
    other <- x[other_name, ]
    differing <- focal != other
    tibble(
      focal_sequence = focal_name,
      comparison_sequence = other_name,
      alignment_position = which(differing),
      focal_base = focal[differing],
      comparison_base = other[differing]
    )
  })
}

differences <- difference_table(
  full_alignment,
  "NARW_B_PCNA_like_core"
)

write_csv(
  differences,
  file.path(result_dir, "NARW_B_core_sequence_differences.csv")
)

message("Analysis complete. Results written to: ", result_dir)
