#!/usr/bin/env Rscript

# Poster figure: identity of each sequenced PCNA-derived core to the
# parent PCNA CDS from the same species.
#
# Required packages:
# install.packages(c("ggplot2", "dplyr", "readr", "tidyr"))

library(ggplot2)
library(dplyr)
library(readr)

# ---- 1. File locations -------------------------------------------------------

comparison_file <- paste0(
  "C:/Users/angel/Downloads/Research project summer/",
  "ampseq_NARwhale.botdolphin/results_expanded/",
  "PCNA_expanded_pairwise_comparisons.csv"
)

output_directory <- paste0(
  "C:/Users/angel/Downloads/Research project summer/",
  "ampseq_NARwhale.botdolphin/results_expanded/poster_figures"
)

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

# ---- 2. Define the biologically appropriate comparisons ---------------------

# Each row identifies one sequenced PCNA-derived core and the parent PCNA CDS
# from the SAME species. Bottlenose sequences are not included because the
# present comparison contains NCBI pseudogene sequences, not sequenced
# bottlenose amplicon consensuses with a species-parent comparison.

wanted_pairs <- tibble::tribble(
  ~species,     ~locus,         ~retrogene_id,                 ~parent_id,
  "NARW",       "PCNA-like",    "NARW_B_PCNA_like_core",       "NARW_parent_PCNA_CDS",
  "Beluga",     "Retrogene 1",  "Beluga_retrogene1_core",      "Beluga_parent_PCNA_CDS",
  "Beluga",     "Retrogene 2",  "Beluga_retrogene2_core",      "Beluga_parent_PCNA_CDS",
  "Blue whale", "Retrogene 2",  "BlueWhale_retrogene2_core",   "BlueWhale_parent_PCNA_CDS"
)

# ---- 3. Read the previously calculated pairwise results ---------------------

pairwise <- read_csv(comparison_file, show_col_types = FALSE) |>
  filter(dataset == "full_core")

# The pairwise table contains each comparison once, but either sequence can
# occur in sequence_1 or sequence_2. Match in both orientations.

plot_data <- wanted_pairs |>
  rowwise() |>
  mutate(
    matching_row = list(
      pairwise |>
        filter(
          (sequence_1 == retrogene_id & sequence_2 == parent_id) |
            (sequence_1 == parent_id & sequence_2 == retrogene_id)
        )
    )
  ) |>
  ungroup() |>
  tidyr::unnest(matching_row)

if (nrow(plot_data) != nrow(wanted_pairs)) {
  stop(
    "Not all four species-matched comparisons were found. ",
    "Check the sequence names and the input CSV."
  )
}

# ---- 4. Prepare labels and plotting order -----------------------------------

plot_data <- plot_data |>
  mutate(
    display_name = paste(species, locus, sep = "\n"),
    display_name = factor(
      display_name,
      levels = rev(c(
        "NARW\nPCNA-like",
        "Beluga\nRetrogene 1",
        "Beluga\nRetrogene 2",
        "Blue whale\nRetrogene 2"
      ))
    ),
    species = factor(
      species,
      levels = c("NARW", "Beluga", "Blue whale")
    ),
    difference_label = if_else(
      gap_columns > 0,
      paste0(
        sprintf("%.1f%%", percent_identity),
        "  (", substitutions, " substitutions, ",
        gap_columns, " gap column)"
      ),
      paste0(
        sprintf("%.1f%%", percent_identity),
        "  (", substitutions, " substitutions)"
      )
    ),
    n_label = paste0("n = ", compared_nucleotides, " compared nt"),
    # Put labels for values near 99% to the left so they remain inside the plot.
    label_to_left = percent_identity >= 98.6,
    label_x = if_else(
      label_to_left,
      percent_identity - 0.08,
      percent_identity + 0.08
    ),
    label_hjust = if_else(label_to_left, 1, 0)
  )

# Save the exact source values used in the figure.
write_csv(
  plot_data |>
    select(
      species, locus, percent_identity, substitutions,
      gap_columns, compared_nucleotides
    ),
  file.path(output_directory, "retrogene_parent_identity_source_data.csv")
)

# ---- 5. Construct the plot ---------------------------------------------------

species_colours <- c(
  "NARW" = "#079BB5",
  "Beluga" = "#4E69B1",
  "Blue whale" = "#8064A2"
)

p <- ggplot(
  plot_data,
  aes(y = display_name, colour = species)
) +
  # A line from the lower plot boundary to the observed identity.
  geom_segment(
    aes(
      x = 95.5,
      xend = percent_identity,
      yend = display_name
    ),
    linewidth = 2.2,
    lineend = "round"
  ) +
  geom_point(
    aes(x = percent_identity),
    size = 5.5
  ) +
  geom_text(
    aes(
      x = label_x,
      label = difference_label,
      hjust = label_hjust
    ),
    colour = "#003B68",
    fontface = "bold",
    size = 4.2,
    nudge_y = 0.12,
    show.legend = FALSE
  ) +
  geom_text(
    aes(
      x = label_x,
      label = n_label,
      hjust = label_hjust
    ),
    colour = "#4D5B64",
    size = 3.4,
    nudge_y = -0.13,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = species_colours) +
  scale_x_continuous(
    limits = c(95.5, 100.15),
    breaks = 96:100,
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = paste(
      "PCNA-derived cores remain highly similar to",
      "species-matched parent PCNA"
    ),
    subtitle = "Descriptive nucleotide comparison",
    x = "Nucleotide identity to species-matched parent PCNA CDS (%)",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(
      colour = "#003B68",
      face = "bold",
      size = 22,
      hjust = 0.5,
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      colour = "#003B68",
      face = "bold",
      size = 12,
      hjust = 0,
      margin = margin(b = 14)
    ),
    axis.title.x = element_text(
      colour = "#003B68",
      face = "bold",
      margin = margin(t = 12)
    ),
    axis.text.x = element_text(
      colour = "#4D5B64",
      size = 12
    ),
    axis.text.y = element_text(
      colour = "#003B68",
      size = 13,
      lineheight = 0.95,
      margin = margin(r = 14)
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(
      colour = "#E8EEF1",
      linewidth = 0.8
    ),
    legend.position = "none",
    plot.margin = margin(18, 24, 18, 18)
  )

# ---- 6. Save high-resolution and vector versions ----------------------------

ggsave(
  filename = file.path(
    output_directory,
    "PCNA_retrogene_parent_identity_poster.png"
  ),
  plot = p,
  width = 11,
  height = 6,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_directory,
    "PCNA_retrogene_parent_identity_poster.pdf"
  ),
  plot = p,
  width = 11,
  height = 6,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

print(p)

message("Figure and source-data CSV written to: ", output_directory)
