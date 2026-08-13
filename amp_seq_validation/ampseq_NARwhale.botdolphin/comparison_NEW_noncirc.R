#!/usr/bin/env Rscript

## Expanded PCNA / PCNA-like comparison -------
# Purpose: comparison of all PCNA retrogene sequences to date
# Inputs:
#   1. Existing NARW + bottlenose comparison FASTA
#   2. Manually compiled beluga + blue-whale FASTA
# Outputs are written to:
#   ampseq_NARwhale.botdolphin/results_expanded/

# Install packages once if needed:
# install.packages(c("ape", "ggplot2", "readr", "dplyr", "purrr", "tidyr"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("Biostrings", "DECIPHER"))

required_packages <- c(
  "Biostrings",
  "DECIPHER",
  "pwalign",
  "ape",
  "ggplot2",
  "readr",
  "dplyr",
  "purrr",
  "tidyr"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Install these packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(Biostrings)
  library(DECIPHER)
  library(ape)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(pwalign)
})

# ----------------------------------------------------
# ------- File locations -----------
# ---------------------------------------------------

existing_file <- paste0(
  "C:/Users/angel/Downloads/Research project summer/",
  "ampseq_NARwhale.botdolphin/",
  "PCNA_comparison_nrwhale_bdol.fasta"
)

narw_original_file <- paste0(
  "C:/Users/angel/Downloads/Research project summer/",
  "ampseq_NARwhale.botdolphin/",
  "compilation_consensus_comparison_NRwhale"
)

new_whale_file <- paste0(
  "C:/Users/angel/Downloads/Research project summer/",
  "amplicon_sequencing_rerun/",
  "consensus_seqs_belugablue"
)

project_dir <- paste0(
  "C:/Users/angel/Downloads/Research project summer/",
  "ampseq_NARwhale.botdolphin"
)

result_dir <- file.path(
  project_dir,
  "results_expanded"
)

expanded_fasta <- file.path(
  project_dir,
  "PCNA_comparison_expanded.fasta"
)

stopifnot(file.exists(existing_file))
stopifnot(file.exists(narw_original_file))
stopifnot(file.exists(new_whale_file))

dir.create(
  result_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# -------------------------------------------------
# -------- Read and validate inputs ----------
# ----------------------------------------------------

existing <- readDNAStringSet(existing_file)
new_whales <- readDNAStringSet(new_whale_file)

narw_original <- readDNAStringSet(
  narw_original_file
)

stopifnot(
  "NARW01_B_denovo" %in%
    names(narw_original),
  
  "REF_PCNA_CDS" %in%
    names(narw_original)
)

if (anyDuplicated(names(existing))) {
  stop("The existing comparison FASTA contains duplicate identifiers.")
}

if (anyDuplicated(names(new_whales))) {
  stop("The beluga/blue-whale FASTA contains duplicate identifiers.")
}

required_existing_names <- c(
  "NARW_A_primer_free",
  "NARW_B_PCNA_like_core",
  "NARW_NCBI_PCNA_like",
  "NARW_parent_PCNA_CDS",
  "Bottlenose_PCNA_pseudogene1_core",
  "Bottlenose_PCNA_pseudogene2_core"
)

required_new_names <- c(
  "Beluga_retrogene1B_vc_denovoref",
  "Beluga_retrogene2B_vc_denovoref",
  "BlueWhale_retrogene2B_vc_denovo",
  "Beluga_parent_PCNA_CDS",
  "BlueWhale_parent_PCNA_CDS"
)

missing_existing_names <- setdiff(
  required_existing_names,
  names(existing)
)

missing_new_names <- setdiff(
  required_new_names,
  names(new_whales)
)

if (length(missing_existing_names) > 0) {
  stop(
    "Missing identifiers from the existing comparison FASTA: ",
    paste(missing_existing_names, collapse = ", ")
  )
}

if (length(missing_new_names) > 0) {
  stop(
    "Missing identifiers from consensus_seqs_belugablue: ",
    paste(missing_new_names, collapse = ", ")
  )
}

# -----------------------------------------------------
# Independently define the NARW B PCNA-homologous core
# The full de novo B consensus is oriented and aligned to parental PCNA CDS instead of NCBI PCNA-like (felt like the latter was biased reasoning)


orient_to_parent <- function(
    query,
    parent_cds
) {
  forward_alignment <- pairwiseAlignment(
    query,
    parent_cds,
    type = "local"
  )
  
  reverse_alignment <- pairwiseAlignment(
    reverseComplement(query),
    parent_cds,
    type = "local"
  )
  
  if (
    score(reverse_alignment) >
    score(forward_alignment)
  ) {
    reverseComplement(query)
  } else {
    query
  }
}

extract_parent_homologous_core <- function(
    query,
    parent_cds
) {
  query_oriented <- orient_to_parent(
    query,
    parent_cds
  )
  
  core_alignment <- AlignSeqs(
    DNAStringSet(c(
      parent_CDS =
        as.character(parent_cds),
      
      de_novo_consensus =
        as.character(query_oriented)
    )),
    processors = NULL,
    verbose = FALSE
  )
  
  core_matrix <- as.matrix(
    core_alignment
  )
  
  parent_columns <- which(
    core_matrix["parent_CDS", ] %in%
      c("A", "C", "G", "T")
  )
  
  homologous_span <- min(parent_columns):
    max(parent_columns)
  
  query_core <- core_matrix[
    "de_novo_consensus",
    homologous_span
  ]
  
  DNAString(
    paste0(
      query_core[query_core != "-"],
      collapse = ""
    )
  )
}

narw_B_full_consensus <- narw_original[[
  "NARW01_B_denovo"
]]

narw_parent_CDS <- narw_original[[
  "REF_PCNA_CDS"
]]

narw_B_core_independent <-
  extract_parent_homologous_core(
    query = narw_B_full_consensus,
    parent_cds = narw_parent_CDS
  )

length(narw_B_core_independent)
#
narw_core_index <- match(
  "NARW_B_PCNA_like_core",
  names(existing)
)

stopifnot(!is.na(narw_core_index))

existing[narw_core_index] <- DNAStringSet(
  as.character(
    narw_B_core_independent
  )
)

names(existing)[narw_core_index] <-
  "NARW_B_PCNA_like_core"

narw_ncbi_like <- existing[[
  "NARW_NCBI_PCNA_like"
]]

narw_ncbi_validation <- tibble(
  independently_defined_core_length =
    length(narw_B_core_independent),
  
  ncbi_reference_length =
    length(narw_ncbi_like),
  
  identical_to_ncbi =
    identical(
      as.character(
        narw_B_core_independent
      ),
      as.character(
        narw_ncbi_like
      )
    )
)

narw_ncbi_validation
# -----------------------------------------------------
# Orient and extract the new PCNA-homologous cores
#
# Coordinates were established by local alignment against each species'
# 786-bp parent PCNA CDS:
#   Beluga retrogene 1: positions 151-935, forward orientation (785 bp)
#   Beluga retrogene 2: positions 52-837, forward orientation (786 bp)
#   Blue-whale retrogene 2: reverse complement, then positions 52-837 (786 bp)
# -----------------------------------------------------
beluga_retrogene1_core <- subseq(
  new_whales[["Beluga_retrogene1B_vc_denovoref"]],
  start = 151,
  end = 935
)

beluga_retrogene2_core <- subseq(
  new_whales[["Beluga_retrogene2B_vc_denovoref"]],
  start = 52,
  end = 837
)

blue_retrogene2_oriented <- reverseComplement(
  new_whales[["BlueWhale_retrogene2B_vc_denovo"]]
)

blue_retrogene2_core <- subseq(
  blue_retrogene2_oriented,
  start = 52,
  end = 837
)

new_comparison_sequences <- DNAStringSet(c(
  Beluga_retrogene1_core =
    as.character(beluga_retrogene1_core),
  
  Beluga_retrogene2_core =
    as.character(beluga_retrogene2_core),
  
  BlueWhale_retrogene2_core =
    as.character(blue_retrogene2_core),
  
  Beluga_parent_PCNA_CDS =
    as.character(
      new_whales[["Beluga_parent_PCNA_CDS"]]
    ),
  
  BlueWhale_parent_PCNA_CDS =
    as.character(
      new_whales[["BlueWhale_parent_PCNA_CDS"]]
    )
))

expected_new_lengths <- c(
  Beluga_retrogene1_core = 785,
  Beluga_retrogene2_core = 786,
  BlueWhale_retrogene2_core = 786,
  Beluga_parent_PCNA_CDS = 786,
  BlueWhale_parent_PCNA_CDS = 786
)

observed_new_lengths <- setNames(
  width(new_comparison_sequences),
  names(new_comparison_sequences)
)

if (!identical(
  as.integer(observed_new_lengths),
  as.integer(expected_new_lengths[names(observed_new_lengths)])
)) {
  stop(
    "Unexpected extracted sequence lengths. Review the input FASTA ",
    "and extraction coordinates."
  )
}

# The 5,957-bp B12BelugaPCNA_AParent genomic locus is intentionally excluded.
expanded <- c(
  existing,
  new_comparison_sequences
)

if (anyDuplicated(names(expanded))) {
  stop("Duplicate identifiers detected after combining the FASTAs.")
}

writeXStringSet(
  expanded,
  expanded_fasta
)

sequence_inventory <- tibble(
  sequence_id = names(expanded),
  length = width(expanded),
  ambiguous_N = vcountPattern(
    "N",
    expanded,
    fixed = TRUE
  )
)

write_csv(
  sequence_inventory,
  file.path(
    result_dir,
    "PCNA_expanded_sequence_inventory.csv"
  )
)

# -----------------------------------------------------
# ----- Multiple-sequence alignments ----
# -----------------------------------------------------

# The 602-bp NARW A sequence is excluded from the full-core alignment.
full_core_names <- setdiff(
  names(expanded),
  "NARW_A_primer_free"
)

full_core <- expanded[full_core_names]

full_alignment <- AlignSeqs(
  full_core,
  processors = NULL,
  verbose = FALSE
)

writeXStringSet(
  full_alignment,
  file.path(
    result_dir,
    "PCNA_expanded_full_core_alignment.fasta"
  )
)

# A second alignment keeping only the alignment span covered by primer-free A.
shared_alignment_all <- AlignSeqs(
  expanded,
  processors = NULL,
  verbose = FALSE
)

shared_matrix_all <- as.matrix(
  shared_alignment_all
)

a_row <- shared_matrix_all[
  "NARW_A_primer_free",
]

a_columns <- a_row %in% c(
  "A",
  "C",
  "G",
  "T"
)

shared_matrix <- shared_matrix_all[
  ,
  a_columns,
  drop = FALSE
]

shared_strings <- apply(
  shared_matrix,
  1,
  paste0,
  collapse = ""
)

shared_alignment <- DNAStringSet(
  shared_strings
)

names(shared_alignment) <- rownames(
  shared_matrix
)

writeXStringSet(
  shared_alignment,
  file.path(
    result_dir,
    "PCNA_expanded_A_shared_region_alignment.fasta"
  )
)

# -----------------------------------------------------
# ---- Pairwise sequence statistics ------
# -----------------------------------------------------

pairwise_stats <- function(
    aligned_sequences,
    dataset_name
) {
  x <- as.matrix(aligned_sequences)
  sequence_names <- rownames(x)
  
  pairs <- combn(
    sequence_names,
    2,
    simplify = FALSE
  )
  
  map_dfr(pairs, function(pair) {
    a <- x[pair[1], ]
    b <- x[pair[2], ]
    
    valid_a <- a %in% c(
      "A",
      "C",
      "G",
      "T"
    )
    
    valid_b <- b %in% c(
      "A",
      "C",
      "G",
      "T"
    )
    
    comparable <- valid_a & valid_b
    gap_columns <- xor(
      a == "-",
      b == "-"
    )
    
    matches <- sum(
      a[comparable] == b[comparable]
    )
    
    substitutions <- sum(
      a[comparable] != b[comparable]
    )
    
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
  pairwise_stats(
    full_alignment,
    "full_core"
  ),
  pairwise_stats(
    shared_alignment,
    "A_shared_region"
  )
)

write_csv(
  pairwise_results,
  file.path(
    result_dir,
    "PCNA_expanded_pairwise_comparisons.csv"
  )
)

# -----------------------------------------------------
# Expanded full-core identity heatmap
# -----------------------------------------------------

display_labels <- c(
  NARW_B_PCNA_like_core =
    "NARW B core",
  
  NARW_NCBI_PCNA_like =
    "NARW NCBI PCNA-like",
  
  NARW_parent_PCNA_CDS =
    "NARW parent PCNA CDS",
  
  Bottlenose_PCNA_pseudogene1_core =
    "Bottlenose pseudogene 1",
  
  Bottlenose_PCNA_pseudogene2_core =
    "Bottlenose pseudogene 2",
  
  Beluga_retrogene1_core =
    "Beluga retrogene 1",
  
  Beluga_retrogene2_core =
    "Beluga retrogene 2",
  
  BlueWhale_retrogene2_core =
    "Blue-whale retrogene 2",
  
  Beluga_parent_PCNA_CDS =
    "Beluga parent PCNA CDS",
  
  BlueWhale_parent_PCNA_CDS =
    "Blue-whale parent PCNA CDS"
)

unlabelled_sequences <- setdiff(
  names(full_alignment),
  names(display_labels)
)

if (length(unlabelled_sequences) > 0) {
  stop(
    "Missing display labels for: ",
    paste(unlabelled_sequences, collapse = ", ")
  )
}

full_pairs <- pairwise_results |>
  filter(dataset == "full_core") |>
  distinct(
    sequence_1,
    sequence_2,
    .keep_all = TRUE
  )

identity_long <- bind_rows(
  full_pairs |>
    select(
      sequence_1,
      sequence_2,
      percent_identity
    ),
  
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
  distinct(
    sequence_1,
    sequence_2,
    .keep_all = TRUE
  )

sequence_order <- names(full_alignment)

identity_long <- identity_long |>
  mutate(
    sequence_1 = factor(
      sequence_1,
      levels = sequence_order,
      labels = unname(
        display_labels[sequence_order]
      )
    ),
    
    sequence_2 = factor(
      sequence_2,
      levels = rev(sequence_order),
      labels = unname(
        display_labels[rev(sequence_order)]
      )
    )
  )

identity_plot <- ggplot(
  identity_long,
  aes(
    sequence_1,
    sequence_2,
    fill = percent_identity
  )
) +
  geom_tile(
    color = "white"
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.1f",
        percent_identity
      )
    ),
    size = 2.7
  ) +
  scale_fill_viridis_c() +
  coord_equal() +
  labs(
    title = "PCNA-derived sequence identity",
    subtitle = paste(
      "Identity across comparable A/C/G/T positions;",
      "gap columns excluded"
    ),
    x = NULL,
    y = NULL,
    fill = "Identity (%)"
  ) +
  theme_minimal(
    base_size = 11
  ) +
  theme(
    axis.text.x = element_text(
      angle = 40,
      hjust = 1
    ),
    panel.grid = element_blank(),
    plot.margin = margin(
      15,
      30,
      20,
      15
    )
  )

ggsave(
  file.path(
    result_dir,
    "PCNA_expanded_identity_heatmap.png"
  ),
  identity_plot,
  width = 12,
  height = 11,
  dpi = 300
)

# -----------------------------------------------------
# Retrogene-only identity heatmap
# -----------------------------------------------------

retrogene_names <- c(
  "NARW_B_PCNA_like_core",
  "Bottlenose_PCNA_pseudogene1_core",
  "Bottlenose_PCNA_pseudogene2_core",
  "Beluga_retrogene1_core",
  "Beluga_retrogene2_core",
  "BlueWhale_retrogene2_core"
)

retrogene_alignment <- full_alignment[
  retrogene_names
]

retrogene_pairs <- pairwise_stats(
  retrogene_alignment,
  "retrogenes_only"
)

write_csv(
  retrogene_pairs,
  file.path(
    result_dir,
    "PCNA_retrogene_only_pairwise_comparisons.csv"
  )
)

retrogene_identity_long <- bind_rows(
  retrogene_pairs |>
    select(
      sequence_1,
      sequence_2,
      percent_identity
    ),
  
  retrogene_pairs |>
    transmute(
      sequence_1 = sequence_2,
      sequence_2 = sequence_1,
      percent_identity
    ),
  
  tibble(
    sequence_1 = retrogene_names,
    sequence_2 = retrogene_names,
    percent_identity = 100
  )
) |>
  distinct(
    sequence_1,
    sequence_2,
    .keep_all = TRUE
  ) |>
  mutate(
    sequence_1 = factor(
      sequence_1,
      levels = retrogene_names,
      labels = unname(
        display_labels[retrogene_names]
      )
    ),
    
    sequence_2 = factor(
      sequence_2,
      levels = rev(retrogene_names),
      labels = unname(
        display_labels[rev(retrogene_names)]
      )
    )
  )

retrogene_identity_plot <- ggplot(
  retrogene_identity_long,
  aes(
    sequence_1,
    sequence_2,
    fill = percent_identity
  )
) +
  geom_tile(
    color = "white"
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.1f",
        percent_identity
      )
    ),
    size = 3
  ) +
  scale_fill_viridis_c() +
  coord_equal() +
  labs(
    title = "PCNA-derived retrogene identity",
    subtitle = paste(
      "Identity across comparable A/C/G/T positions;",
      "gap columns excluded"
    ),
    x = NULL,
    y = NULL,
    fill = "Identity (%)"
  ) +
  theme_minimal(
    base_size = 11
  ) +
  theme(
    axis.text.x = element_text(
      angle = 40,
      hjust = 1
    ),
    panel.grid = element_blank()
  )

ggsave(
  file.path(
    result_dir,
    "PCNA_retrogene_only_identity_heatmap.png"
  ),
  retrogene_identity_plot,
  width = 9,
  height = 8,
  dpi = 300
)

# -----------------------------------------------------
# --------- Exploratory neighbor-joining tree --------
# -----------------------------------------------------

tree_alignment <- full_alignment

narw_sequences_identical <- identical(
  as.character(
    full_alignment[[
      "NARW_B_PCNA_like_core"
    ]]
  ),
  as.character(
    full_alignment[[
      "NARW_NCBI_PCNA_like"
    ]]
  )
)

if (narw_sequences_identical) {
  tree_alignment <- tree_alignment[
    names(tree_alignment) !=
      "NARW_NCBI_PCNA_like"
  ]
  
  names(tree_alignment)[
    names(tree_alignment) ==
      "NARW_B_PCNA_like_core"
  ] <- "NARW_sample_and_NCBI_PCNA_like_identical"
}

tree_display_labels <- c(
  NARW_sample_and_NCBI_PCNA_like_identical =
    "NARW sample + NCBI PCNA-like",
  
  NARW_B_PCNA_like_core =
    "NARW B core",
  
  NARW_NCBI_PCNA_like =
    "NARW NCBI PCNA-like",
  
  NARW_parent_PCNA_CDS =
    "NARW parent PCNA CDS",
  
  Bottlenose_PCNA_pseudogene1_core =
    "Bottlenose pseudogene 1",
  
  Bottlenose_PCNA_pseudogene2_core =
    "Bottlenose pseudogene 2",
  
  Beluga_retrogene1_core =
    "Beluga retrogene 1",
  
  Beluga_retrogene2_core =
    "Beluga retrogene 2",
  
  BlueWhale_retrogene2_core =
    "Blue-whale retrogene 2",
  
  Beluga_parent_PCNA_CDS =
    "Beluga parent PCNA CDS",
  
  BlueWhale_parent_PCNA_CDS =
    "Blue-whale parent PCNA CDS"
)

unlabelled_tree_tips <- setdiff(
  names(tree_alignment),
  names(tree_display_labels)
)

if (length(unlabelled_tree_tips) > 0) {
  stop(
    "Missing tree labels for: ",
    paste(unlabelled_tree_tips, collapse = ", ")
  )
}

tree_dnabin <- as.DNAbin(
  as.matrix(tree_alignment)
)

tree_distance <- dist.dna(
  tree_dnabin,
  model = "raw",
  pairwise.deletion = TRUE
)

tree <- nj(
  tree_distance
)

tree$tip.label <- unname(
  tree_display_labels[
    tree$tip.label
  ]
)

if (anyNA(tree$tip.label)) {
  stop("At least one tree label became NA.")
}

write.tree(
  tree,
  file = file.path(
    result_dir,
    "PCNA_expanded_NJ_tree.newick"
  )
)

png(
  file.path(
    result_dir,
    "PCNA_expanded_NJ_tree.png"
  ),
  width = 3300,
  height = 2500,
  res = 300
)

par(
  mar = c(
    4,
    3,
    4,
    12
  )
)

plot(
  tree,
  type = "phylogram",
  direction = "rightwards",
  cex = 0.85,
  label.offset = 0.001,
  main = paste(
    "Exploratory NJ tree of",
    "PCNA-derived sequences"
  )
)

add.scale.bar()

dev.off()

# -----------------------------------------------------
# ------ Detailed aligned differences relative to the NARW B core ----------
# -----------------------------------------------------

difference_table <- function(
    aligned_sequences,
    focal_name
) {
  x <- as.matrix(aligned_sequences)
  focal <- x[focal_name, ]
  
  map_dfr(
    setdiff(
      rownames(x),
      focal_name
    ),
    function(other_name) {
      other <- x[other_name, ]
      differing <- focal != other
      
      tibble(
        focal_sequence = focal_name,
        comparison_sequence = other_name,
        alignment_position = which(differing),
        focal_base = focal[differing],
        comparison_base = other[differing]
      )
    }
  )
}

differences <- difference_table(
  full_alignment,
  "NARW_B_PCNA_like_core"
)

write_csv(
  differences,
  file.path(
    result_dir,
    "NARW_B_core_expanded_sequence_differences.csv"
  )
)

message(
  "Expanded PCNA analysis complete."
)

message(
  "Expanded FASTA: ",
  expanded_fasta
)

message(
  "Results directory: ",
  result_dir
)
# --- end ----


