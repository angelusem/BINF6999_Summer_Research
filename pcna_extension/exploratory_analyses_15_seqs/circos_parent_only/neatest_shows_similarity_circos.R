# ============================================================
# PCNA pairwise similarity Circos plot
# ============================================================

library(ape)
library(circlize)

# ------------------------------------------------------------
# 1. Read alignment
# ------------------------------------------------------------

alignment <- read.dna(
  "circos_parent_only/aligned_PCNA_parent_only_sep30.fasta",
  format = "fasta"
)

# ------------------------------------------------------------
# 2. Calculate pairwise nucleotide identity
# ------------------------------------------------------------

# raw = proportion of sites that differ
raw_dist <- as.matrix(
  dist.dna(
    alignment,
    model = "raw",
    pairwise.deletion = TRUE
  )
)

# Convert proportion different -> percent identical
sim_matrix <- (1 - raw_dist) * 100

# We do not want self-comparisons
diag(sim_matrix) <- NA


# ------------------------------------------------------------
# 3. Convert matrix to pairwise table
# ------------------------------------------------------------

species <- rownames(sim_matrix)
# ------------------------------------------------------------
# CLEAN DISPLAY LABELS
# ------------------------------------------------------------

# Keep original FASTA names for calculations
# but create prettier names for the figure.

pretty_labels <- c(
  "BlueWhale_parent"         = "Blue Whale",
  "NARightWhale_parent"      = "N.A. right Whale",
  "Bowhead_parent"           = "Bowhead Whale",
  "HumpbackWhale_parent"     = "Humpback Whale",
  "GrayWhale_parent"         = "Gray Whale",
  "MinkeWhale_parent"        = "Minke Whale",
  
  "BelugaWhale_parent"       = "Beluga Whale",
  "BottlenosDolphin_parent"  = "Bottlenose Dolphin",
  "Narwhal_parent"           = "Narwhal",
  "Orca_parent"              = "Orca",
  "SpermWhale_parent"        = "Sperm Whale",
  "PygmySpermWhale_parent"   = "Pygmy Sperm Whale"
)

pairwise <- data.frame(
  species1 = character(),
  species2 = character(),
  identity = numeric()
)

for (i in 1:(length(species) - 1)) {
  for (j in (i + 1):length(species)) {
    
    pairwise <- rbind(
      pairwise,
      data.frame(
        species1 = species[i],
        species2 = species[j],
        identity = sim_matrix[i, j]
      )
    )
  }
}

# Remove comparisons that could not be calculated
pairwise <- pairwise[!is.na(pairwise$identity), ]


# ------------------------------------------------------------
# 4. Inspect actual identity values
# ------------------------------------------------------------

range(pairwise$identity)

summary(pairwise$identity)

# Optional:
# See every comparison from least -> most similar
pairwise[order(pairwise$identity), ]


# ------------------------------------------------------------
# 5. Define whale groups
# ------------------------------------------------------------
#
# IMPORTANT:
# Replace the names below with the EXACT sequence names
# from your FASTA file.
#
# Check them with:
rownames(alignment)
baleen <- c(
  "BlueWhale_parent",
  "NARightWhale_parent",
  "Bowhead_parent",
  "HumpbackWhale_parent",
  "GrayWhale_parent",
  "MinkeWhale_parent"
)

toothed <- c(
  "BelugaWhale_parent",
  "BottlenosDolphin_parent",
  "Narwhal_parent",
  "SpermWhale_parent",
  "Orca_parent",
  "PygmySpermWhale_parent"
)

group <- ifelse(
  species %in% baleen,
  "Baleen",
  ifelse(
    species %in% toothed,
    "Toothed",
    "Other"
  )
)

names(group) <- species


# ------------------------------------------------------------
# 6. Sector colors = whale group
# ------------------------------------------------------------

# Change these colors however you want
sector_colors <- ifelse(
  group == "Baleen", "#bdbdbd",
  ifelse(
    group == "Toothed", "#636363",
    "grey70"
  )
)

names(sector_colors) <- species


# ------------------------------------------------------------
# 7. Similarity color scale
# ------------------------------------------------------------
#
# THIS is what communicates sequence similarity.
#
# Your alignment is approximately 98.73–100% identical,
# so use the actual observed range rather than 0–100%.

identity_min <- min(pairwise$identity)
identity_max <- max(pairwise$identity)

similarity_col_fun <- colorRamp2(
  c(
    identity_min,
    (identity_min + identity_max) / 2,
    identity_max
  ),
  c(
    "#FFF7BC",   # relatively lower identity
    "#41B6C4",
    "#253494"    # highest identity
  )
)

# Give each pair its color based on % identity
pairwise$link_color <- similarity_col_fun(pairwise$identity)


# ------------------------------------------------------------
# 8. 
# Give every ribbon the SAME weight
# ------------------------------------------------------------
#
# We don't want 98.8 vs 99.8 controlling ribbon width because
# visually those values are almost identical.
#
# Instead:
#     ribbon width = constant
#     ribbon color = identity
#

links <- data.frame(
  from  = pairwise$species1,
  to    = pairwise$species2,
  value = 1
)


# ------------------------------------------------------------
# 9. Choose species order
# ------------------------------------------------------------
#
# This places baleen together, followed by toothed whales.

species_order <- c(
  species[group == "Baleen"],
  species[group == "Toothed"],
  species[group == "Other"]
)


# ------------------------------------------------------------
# 10. Draw Circos plot
# ------------------------------------------------------------

circos.clear()

circos.par(
  start.degree = 90,
  gap.degree = 2.5,
  track.margin = c(0.01, 0.01),
  cell.padding = c(0, 0, 0, 0)
)

chordDiagram(
  x = links,
  
  order = species_order,
  
  grid.col = sector_colors,
  
  # Ribbon color = pairwise identity
  col = pairwise$link_color,
  
  # Transparency helps with crossings
  transparency = 0.30,
  
  annotationTrack = "grid",
  
  preAllocateTracks = list(
    track.height = 0.14
  ),
  
  
  link.sort = TRUE,
  
  link.decreasing = FALSE
)


# ------------------------------------------------------------
# 11. Add species labels
# ------------------------------------------------------------

# ============================================================
# CURVED SPECIES LABELS
# ============================================================

circos.trackPlotRegion(
  track.index = 1,

  panel.fun = function(x, y) {

    original_name <- CELL_META$sector.index

    # Get clean display name
    display_name <- pretty_labels[original_name]

    # Fall back to original name if no pretty label exists
    if (is.na(display_name)) {
      display_name <- original_name
    }

    circos.text(
      x = CELL_META$xcenter,

      # Position near outer part of label track
      y = CELL_META$ycenter,

      labels = display_name,

      # Make text follow circular curvature
      facing = "bending.outside",

      niceFacing = TRUE,

      cex = 0.7,
      font = 
    )
  },

  bg.border = NA
)

# ------------------------------------------------------------
# 12. Add similarity legend
# ------------------------------------------------------------

# ------------------------------------------------------------
# 12. Add similarity legend 
# ------------------------------------------------------------

baleen_color <- "grey75"
toothed_color <- "grey35"

sector_colors <- ifelse(
  group == "Baleen",
  baleen_color,
  toothed_color
)

names(sector_colors) <- species

legend(
  "topleft",
  
  legend = c(
    "Baleen whale",
    "Toothed whale"
  ),
  
  fill = c(
    baleen_color,
    toothed_color
  ),
  
  border = NA,
  
  title = "Cetacean group",
  
  bty = "n",
  
  cex = 0.9
)

legend_values <- seq(
  identity_min,
  identity_max,
  length.out = 5
)

legend_colors <- similarity_col_fun(legend_values)

legend(
  "topright",
  legend = sprintf("%.2f%%", legend_values),
  fill = legend_colors,
  border = NA,
  title = "Pairwise identity",
  bty = "n",
  cex = 0.8
)


circos.clear()

