library(ape)
library(circlize)

# 1. Read the alignment
alignment <- read.dna(
  "circos_parent_only/aligned_PCNA_parent_only_sep30.fasta",
  format = "fasta"
)

# 2. Select parent CDS sequences using their exact FASTA labels
# Replace these examples with your actual sequence names
parent_names <- c(
  "BelugaWhale_parent",
  "BlueWhale_parent",
  "BottlenosDolphin_parent",
  "NARightWhale_parent",
  "Bowhead_parent",
  "Narwhal_parent",
  "PygmySpermWhale_parent",
  "SpermWhale_parent",
  "HumpbackWhale_parent",
  "GrayWhale_parent",
  "MinkeWhale_parent",
  "Orca_parent"
)

stopifnot(
  length(parent_names) >= 2,
  !anyDuplicated(parent_names),
  all(parent_names %in% rownames(alignment))
)

parents <- alignment[parent_names, , drop = FALSE]

# 3. Calculate percentage nucleotide differences
divergence_pct <- 100 * as.matrix(
  dist.dna(
    parents,
    model = "raw",
    pairwise.deletion = TRUE
  )
)
diag(divergence_pct) <- 0

# 4. Count comparable A/C/G/T sites for each pair
bases <- toupper(as.character(parents))

valid <- matrix(
  bases %in% c("A", "C", "G", "T"),
  nrow = nrow(bases),
  ncol = ncol(bases),
  dimnames = dimnames(bases)
)

comparable_sites <- tcrossprod(valid * 1)

# Save complete results before applying display filters
write.csv(divergence_pct, "parent_CDS_divergence_pct.csv")
write.csv(comparable_sites, "parent_CDS_comparable_sites.csv")

# 5. Prepare a plotting copy
# Choose a minimum overlap appropriate for your CDS length
expected_cds_length <- 786  # Replace with your reference CDS length
min_sites <- ceiling(0.90 * expected_cds_length)

plot_matrix <- divergence_pct
plot_matrix[
  !is.finite(divergence_pct) |
    comparable_sites < min_sites
] <- 0
diag(plot_matrix) <- 0

# Optional divergence filter; leave disabled initially
# plot_matrix[plot_matrix < 1 | plot_matrix > 15] <- 0

# Remove sequences without drawable links
active <- rowSums(plot_matrix) > 0

if (any(!active)) {
  message(
    "No drawable links for: ",
    paste(rownames(plot_matrix)[!active], collapse = ", ")
  )
}

plot_matrix <- plot_matrix[active, active, drop = FALSE]

if (nrow(plot_matrix) < 2) {
  stop("No positive differences meet the overlap requirement.")
}

# 6. Assign colors and preserve your chosen sequence order
sequence_names <- rownames(plot_matrix)

grid_colors <- setNames(
  hcl.colors(length(sequence_names), palette = "Dark 3"),
  sequence_names
)

# 7. Draw the chord diagram
circos.clear()
circos.par(
  start.degree = 90,
  gap.after = 3,
  canvas.xlim = c(-1.3, 1.3),
  canvas.ylim = c(-1.3, 1.3)
)

chordDiagram(
  x = plot_matrix,
  symmetric = TRUE,
  directional = 0,
  order = sequence_names,
  grid.col = grid_colors,
  col = "#666666",  # Neutral link color; width encodes divergence
  transparency = 0.5,
  link.sort = TRUE,
  link.decreasing = TRUE,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.18)
)

circos.track(
  track.index = 1,
  panel.fun = function(x, y) {
    circos.text(
      x = CELL_META$xcenter,
      y = CELL_META$ylim[1] + 0.05,
      labels = CELL_META$sector.index,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.6
    )
  },
  bg.border = NA
)

title(
  main = "Parent-gene CDS nucleotide differences",
  sub = paste0(
    "Wider chords = greater difference; minimum ",
    min_sites, " comparable bases"
  )
)

circos.clear()