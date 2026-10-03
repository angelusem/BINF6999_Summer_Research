
library(ape)
library(circlize)

folder <- paste0(
  "C:/Users/angel/Downloads/Research project summer/",
  "exploratory_analyses_15_seqs"
)

metadata <- read.csv(
  file.path(folder, "pcna_circos_subset_metadata.csv"),
  stringsAsFactors = FALSE
)

dna <- read.dna(
  file.path(
    folder,
    "pcna_rg_circos_subset_aligned_fullbwhead_potncbiNarwal.fasta"
  ),
  format = "fasta",
  as.matrix = TRUE
)

# Convert the alignment to a matrix of nucleotide characters.
aln <- as.character(dna)

stopifnot(
  !anyDuplicated(metadata$Identifier),
  !anyDuplicated(rownames(aln)),
  all(metadata$Identifier %in% rownames(aln))
)

# Group sequences by species, placing parents first.
metadata <- metadata[
  order(metadata$Species, metadata$Type != "Parent", metadata$Copy),
]

aln <- aln[metadata$Identifier, , drop = FALSE]

# Short, readable labels for the plot.
metadata$Label <- paste(metadata$Species, metadata$Copy)

out_dir <- file.path(folder, "circos_outputs")
dir.create(out_dir, showWarnings = FALSE)

## defining all possible comparison pairs
pair_index <- t(combn(seq_len(nrow(metadata)), 2))

pairs <- data.frame(
  id1 = metadata$Identifier[pair_index[, 1]],
  id2 = metadata$Identifier[pair_index[, 2]],
  species1 = metadata$Species[pair_index[, 1]],
  species2 = metadata$Species[pair_index[, 2]],
  type1 = metadata$Type[pair_index[, 1]],
  type2 = metadata$Type[pair_index[, 2]],
  stringsAsFactors = FALSE
)

pairs$within_species <- pairs$species1 == pairs$species2

table(pairs$within_species)
## calculate identiy in alignment windows
# identity defined as: matching A/C/G/T positions ÷ positions where both sequences have A/C/G/T
window_size <- 60

window_starts <- seq(1, ncol(aln), by = window_size)
valid_bases <- c("a", "c", "g", "t")

calculate_windows <- function(pair_table) {
  
  results <- lapply(seq_len(nrow(pair_table)), function(i) {
    
    a <- tolower(aln[pair_table$id1[i], ])
    b <- tolower(aln[pair_table$id2[i], ])
    
    do.call(rbind, lapply(window_starts, function(start) {
      
      end <- min(start + window_size - 1, ncol(aln))
      positions <- start:end
      
      comparable <- (
        a[positions] %in% valid_bases &
          b[positions] %in% valid_bases
      )
      
      n_comparable <- sum(comparable)
      
      identity <- if (n_comparable > 0) {
        mean(a[positions][comparable] == b[positions][comparable])
      } else {
        NA_real_
      }
      
      data.frame(
        id1 = pair_table$id1[i],
        id2 = pair_table$id2[i],
        start = start,
        end = end,
        n_comparable = n_comparable,
        coverage = n_comparable / length(positions),
        identity = identity
      )
    }))
  })
  
  do.call(rbind, results)
}

windows <- calculate_windows(pairs)

write.csv(
  windows,
  file.path(out_dir, "pairwise_window_identity.csv"),
  row.names = FALSE
)
# all nucleotide comparisons
## defining plotting function
#

plot_comparisons <- function(pair_table, filename, plot_title,
                             min_coverage = 0.5) {
  
  stopifnot(nrow(pair_table) > 0)
  
  ids <- unique(c(pair_table$id1, pair_table$id2))
  m <- metadata[metadata$Identifier %in% ids, , drop = FALSE]
  ids <- m$Identifier
  
  # Select windows belonging to the requested pairs.
  pair_key <- function(a, b) paste(a, b, sep = "|")
  
  wanted <- pair_key(pair_table$id1, pair_table$id2)
  
  w <- windows[
    pair_key(windows$id1, windows$id2) %in% wanted,
    ,
    drop = FALSE
  ]
  
  w <- w[
    !is.na(w$identity) & w$coverage >= min_coverage,
    ,
    drop = FALSE
  ]
  
  # Fixed identity scale, shared by every plot.
  identity_color <- colorRamp2(
    c(0, 0.8, 1),
    c("#2166AC", "#F7F7F7", "#B2182B")
  )
  
  species_levels <- unique(metadata$Species)
  
  species_colors <- setNames(
    hcl.colors(length(species_levels), palette = "Dark 3"),
    species_levels
  )
  
  # Larger gaps separate species.
  n <- nrow(m)
  last_in_species <- c(
    m$Species[-n] != m$Species[-1],
    TRUE
  )
  gaps <- ifelse(last_in_species, 8, 2)
  
  pdf(filename, width = 10, height = 10)
  
  on.exit({
    circos.clear()
    dev.off()
  }, add = TRUE)
  
  par(mar = c(2, 2, 4, 2))
  
  circos.clear()
  circos.par(
    start.degree = 90,
    gap.after = gaps,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.01, 0.01)
  )
  
  circos.initialize(
    factors = factor(ids, levels = ids),
    xlim = cbind(rep(0, n), rep(ncol(aln), n))
  )
  
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = 0.14,
    bg.border = NA,
    panel.fun = function(x, y) {
      
      id <- CELL_META$sector.index
      j <- match(id, m$Identifier)
      
      circos.rect(
        0, 0, ncol(aln), 0.35,
        col = species_colors[m$Species[j]],
        border = NA
      )
      
      circos.text(
        CELL_META$xcenter, 0.72,
        m$Label[j],
        facing = "bending.inside",
        niceFacing = TRUE,
        cex = 0.65
      )
      
      circos.axis(
        h = "bottom",
        major.at = seq(0, ncol(aln), by = 300),
        labels.cex = 0.45,
        minor.ticks = 0
      )
    }
  )
  
  for (i in seq_len(nrow(w))) {
    
    # Alignment column 1 is represented by interval [0, 1].
    interval <- c(w$start[i] - 1, w$end[i])
    
    circos.link(
      w$id1[i], interval,
      w$id2[i], interval,
      col = identity_color(w$identity[i]),
      border = NA
    )
  }
  
  title(main = plot_title)
  
  legend(
    "topleft",
    legend = c("0%", "80%", "100%"),
    fill = identity_color(c(0, 0.8, 1)),
    title = "Nucleotide identity",
    bty = "n",
    cex = 0.8
  )
  
  mtext(
    "Coordinates: alignment columns; windows with <50% comparable sites omitted",
    side = 1,
    cex = 0.7
  )
}
# within species plot
within_pairs <- subset(pairs, within_species)

for (species in unique(metadata$Species)) {
  
  selected <- subset(
    within_pairs,
    species1 == species
  )
  
  safe_name <- gsub("[^A-Za-z0-9]+", "_", species)
  
  plot_comparisons(
    selected,
    file.path(out_dir, paste0("within_", safe_name, ".pdf")),
    paste("Within species:", species)
  )
}
# comparison 1 for test 
species_a <- "Beluga"
species_b <- "Narwhal"

selected <- subset(
  pairs,
  (species1 == species_a & species2 == species_b) |
    (species1 == species_b & species2 == species_a)
)

plot_comparisons(
  selected,
  file.path(out_dir, "between_Beluga_Narwhal.pdf"),
  "Between species: Beluga and Narwhal"
)