#!/usr/bin/env Rscript

## Retrogene-only PCNA identity heatmap ------
# Purpose:
#   1. reads the existing retrogene-only pairwise comparison CSV
#   2. confirms that all 6 expected sequences and all 15 unique pairs exist
#   3. constructs the mirrored half without overwriting sequence names
#   4. adds one 100% self-comparison for every sequence
#   5. saves a full symmetric heatmap and its plotting table.
#
# Required packages:
# install.packages(c("ggplot2", "dplyr", "readr", "tibble"))

library(ggplot2)
library(dplyr)
library(readr)
library(tibble)

# ---- File locations ------

input_file <- paste0(
  "C:/Users/angel/Downloads/Research project summer/",
  "ampseq_NARwhale.botdolphin/results_expanded/",
  "PCNA_retrogene_only_pairwise_comparisons.csv"
)

output_directory <- paste0(
  "C:/Users/angel/Downloads/Research project summer/",
  "ampseq_NARwhale.botdolphin/results_expanded/",
  "corrected_retrogene_heatmap"
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

# ---- Sequence order and display label -----------

sequence_order <- c(
  "NARW_B_PCNA_like_core",
  "Bottlenose_PCNA_pseudogene1_core",
  "Bottlenose_PCNA_pseudogene2_core",
  "Beluga_retrogene1_core",
  "Beluga_retrogene2_core",
  "BlueWhale_retrogene2_core"
)

display_labels <- c(
  NARW_B_PCNA_like_core =
    "NA right whale RG",

  Bottlenose_PCNA_pseudogene1_core =
    "Bottlenose RG 1 (NCBI)",

  Bottlenose_PCNA_pseudogene2_core =
    "Bottlenose RG 2 (NCBI)",

  Beluga_retrogene1_core =
    "Beluga RG 1",

  Beluga_retrogene2_core =
    "Beluga RG 2",

  BlueWhale_retrogene2_core =
    "Blue whale RG"
)

# ---- Read and validate pairwise comparisons --------

pairwise <- read_csv(
  input_file,
  show_col_types = FALSE
) |>
  filter(dataset == "retrogenes_only") |>
  distinct(
    sequence_1,
    sequence_2,
    .keep_all = TRUE
  )

observed_sequences <- union(
  pairwise$sequence_1,
  pairwise$sequence_2
)

missing_sequences <- setdiff(
  sequence_order,
  observed_sequences
)

unexpected_sequences <- setdiff(
  observed_sequences,
  sequence_order
)

if (length(missing_sequences) > 0) {
  stop(
    "The following expected sequences are missing from the CSV: ",
    paste(missing_sequences, collapse = ", ")
  )
}

if (length(unexpected_sequences) > 0) {
  stop(
    "The CSV contains unexpected sequences: ",
    paste(unexpected_sequences, collapse = ", ")
  )
}

expected_pair_count <- choose(
  length(sequence_order),
  2
)

if (nrow(pairwise) != expected_pair_count) {
  stop(
    "Expected ",
    expected_pair_count,
    " unique pairwise comparisons for 6 sequences, but found ",
    nrow(pairwise),
    "."
  )
}

# ---- Build a correct full symmetric matrix ---------

# Temporary names are used when swapping the sequence columns: to avoid dplyr's sequential evaluation changing both columns to the same sequence, which caused incorrect diagonal cells in older draft of script for plot.

mirrored_pairs <- pairwise |>
  transmute(
    sequence_1_new = sequence_2,
    sequence_2_new = sequence_1,
    percent_identity
  ) |>
  rename(
    sequence_1 = sequence_1_new,
    sequence_2 = sequence_2_new
  )

self_comparisons <- tibble(
  sequence_1 = sequence_order,
  sequence_2 = sequence_order,
  percent_identity = 100
)

heatmap_data <- bind_rows(
  pairwise |>
    select(
      sequence_1,
      sequence_2,
      percent_identity
    ),
  mirrored_pairs,
  self_comparisons
) |>
  distinct(
    sequence_1,
    sequence_2,
    .keep_all = TRUE
  )

expected_cell_count <- length(sequence_order)^2

if (nrow(heatmap_data) != expected_cell_count) {
  stop(
    "The completed heatmap should contain ",
    expected_cell_count,
    " cells, but contains ",
    nrow(heatmap_data),
    "."
  )
}

# Verifying that every diagonal value is exactly 100%.
diagonal_check <- heatmap_data |>
  filter(sequence_1 == sequence_2)

if (
  nrow(diagonal_check) != length(sequence_order) ||
  any(diagonal_check$percent_identity != 100)
) {
  stop("The self-comparison diagonal was not constructed correctly.")
}

# Set the display order: Using the same underlying sequence order on both axes to produce a conventional top-left to bottom-right 100% diagonal.

heatmap_data <- heatmap_data |>
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

# Saving the exact 36-cell table to be used for plotting.
write_csv(
  heatmap_data |>
    mutate(
      sequence_1 = as.character(sequence_1),
      sequence_2 = as.character(sequence_2)
    ),
  file.path(
    output_directory,
    "corrected_retrogene_heatmap_values.csv"
  )
)

# ---- Draw the heatmap -------------------

identity_minimum <- floor(
  min(heatmap_data$percent_identity)
)

heatmap_plot <- ggplot(
  heatmap_data,
  aes(
    x = sequence_1,
    y = sequence_2,
    fill = percent_identity
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.1f",
        percent_identity
      )
    ),
    colour = "black",
    size = 4
  ) +
  scale_fill_viridis_c(
    name = "Identity (%)",
    limits = c(
      identity_minimum,
      100
    ),
    breaks = seq(
      identity_minimum,
      100,
      by = 1
    ),
    option = "D"
  ) +
  coord_equal() +
  labs(
    title = "PCNA-derived retrogene sequence identity",
    subtitle = paste(
      "Identity across comparable A/C/G/T positions;",
      "gap columns excluded"
    ),
    x = NULL,
    y = NULL,
    caption = paste(
      "Each sequence is compared with every other sequence.",
      "Diagonal cells are self-comparisons."
    )
  ) +
  theme_minimal(
    base_size = 13
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 40,
      hjust = 1,
      vjust = 1,
      colour = "#003B68",
      size = 11
    ),
    axis.text.y = element_text(
      colour = "#003B68",
      size = 11
    ),
    plot.title = element_text(
      colour = "#003B68",
      face = "bold",
      size = 19
    ),
    plot.subtitle = element_text(
      colour = "#4D5B64",
      size = 11
    ),
    plot.caption = element_text(
      colour = "#4D5B64",
      hjust = 0,
      size = 9
    ),
    legend.title = element_text(
      face = "bold"
    ),
    plot.margin = margin(
      15,
      25,
      15,
      15
    )
  )

# ---- Save high-resolution and vector versions ------------

png_file <- file.path(
  output_directory,
  "PCNA_retrogene_only_heatmap_corrected.png"
)

pdf_file <- file.path(
  output_directory,
  "PCNA_retrogene_only_heatmap_corrected.pdf"
)

ggsave(
  filename = png_file,
  plot = heatmap_plot,
  width = 10,
  height = 8.5,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = heatmap_plot,
  width = 10,
  height = 8.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

print(heatmap_plot)

message(
  "Corrected heatmap written to: ",
  output_directory
)


# ---- Create nonredundant triangular heatmap data ---------

# Blue whale is excluded from the x-axis because its column would be empty.
x_sequence_order <- sequence_order[
  -length(sequence_order)
]

# NARWhale is excluded from the y-axis because its row would be empty.
y_sequence_order <- sequence_order[-1]

# The pairwise object already contains each unique comparison once.
heatmap_data_nonredundant <- pairwise |>
  select(
    sequence_1,
    sequence_2,
    percent_identity
  ) |>
  mutate(
    sequence_1 = factor(
      sequence_1,
      levels = x_sequence_order,
      labels = unname(
        display_labels[x_sequence_order]
      )
    ),
    
    sequence_2 = factor(
      sequence_2,
      levels = rev(y_sequence_order),
      labels = unname(
        display_labels[rev(y_sequence_order)]
      )
    )
  )
#plotting it 
nonredundant_heatmap <- ggplot(
  heatmap_data_nonredundant,
  aes(
    x = sequence_1,
    y = sequence_2,
    fill = percent_identity
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(
      label = sprintf("%.1f", percent_identity),
      colour = percent_identity >= 97.5
    ),
    size = 7,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c(
      `FALSE` = "white",
      `TRUE` = "#003B68"
    )
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_viridis_c(
    name = "Identity (%)",
    limits = c(
      floor(min(
        heatmap_data_nonredundant$percent_identity
      )),
      100
    )
  ) +
  coord_equal() +
  labs(
    title = "PCNA-derived retrogene identity",
    subtitle = paste(
      "Each unique sequence pair is shown once;",
      "self-comparisons omitted"
    ),
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    
    plot.title = element_text(
      size = 27,
      face = "bold",
      colour = "#003B68"
    ),
    
    plot.subtitle = element_text(
      size = 18,
      colour = "#4D5B64"
    ),
    
    axis.text.x = element_text(
      size = 18,
      face = "bold",
      angle = 40,
      hjust = 1,
      colour = "#003B68"
    ),
    
    axis.text.y = element_text(
      size = 18,
      face = "bold",
      colour = "#003B68"
    ),
    
    legend.title = element_text(
      size = 18,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 16
    ),
    
    plot.margin = margin(
      20, 25, 20, 20
    )
  )


print(nonredundant_heatmap)

ggsave(
  filename = file.path(
    output_directory,
    "PCNA_retrogene_heatmap_poster.png"
  ),
  plot = nonredundant_heatmap,
  width = 12.77,
  height = 9.42,
  units = "in",
  dpi = 400,
  bg = "white"
)
