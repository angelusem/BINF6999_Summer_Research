# Run once if needed:
# install.packages(c("ape", "circlize"))

library(ape)
library(circlize)

folder <- paste0(
  "C:/Users/angel/Downloads/Research project summer/",
  "exploratory_analyses_15_seqs"
)

meta <- read.csv(
  file.path(folder, "pcna_circos_subset_metadata.csv"),
  stringsAsFactors = FALSE
)

dna <- read.dna(
  file.path(
    folder,
    "Aligned_pcna_circos_subset_ncbinarwhal.fasta"
  ),
  format = "fasta",
  as.matrix = TRUE
)

aln <- as.character(dna)
 
stopifnot(
  !anyDuplicated(meta$Identifier),
  all(meta$Identifier %in% rownames(aln))
)

# Group by species, with parents before retrogenes.
meta <- meta[
  order(meta$Species, meta$Type != "Parent", meta$Copy),
]

aln <- aln[meta$Identifier, , drop = FALSE]
aln[] <- tolower(aln)

meta$Label <- paste(meta$Species, meta$Copy)


## classify each alignment position
site_colors <- c(
  "Match" = "#EEEEEE",
  "Substitution" = "#D73027",
  "Gap in B" = "#4575B4",
  "Gap in A" = "#FDAE61",
  "Both gaps" = "#FFFFFF",
  "Ambiguous" = "#969696"
)

classify_sites <- function(id_a, id_b) {
  
  a <- aln[id_a, ]
  b <- aln[id_b, ]
  
  # Treat either gap notation as a gap.
  a[a == "."] <- "-"
  b[b == "."] <- "-"
  
  a_base <- a %in% c("a", "c", "g", "t")
  b_base <- b %in% c("a", "c", "g", "t")
  
  state <- rep("Ambiguous", length(a))
  
  state[a_base & b_base & a == b] <- "Match"
  state[a_base & b_base & a != b] <- "Substitution"
  state[a_base & b == "-"] <- "Gap in B"
  state[a == "-" & b_base] <- "Gap in A"
  state[a == "-" & b == "-"] <- "Both gaps"
  
  state
}
#create my comparison table
idx <- t(combn(seq_len(nrow(meta)), 2))

pairs <- data.frame(
  A = meta$Identifier[idx[, 1]],
  B = meta$Identifier[idx[, 2]],
  species_A = meta$Species[idx[, 1]],
  species_B = meta$Species[idx[, 2]],
  label_A = meta$Label[idx[, 1]],
  label_B = meta$Label[idx[, 2]],
  stringsAsFactors = FALSE
)

pairs$within <- pairs$species_A == pairs$species_B
## draw circos
#
plot_differences <- function(
    p,
    title = "Aligned sequence differences"
) {
  stopifnot(nrow(p) > 0)
  
  if (nrow(p) > 15) {
    stop("Select 15 or fewer comparisons.")
  }
  
  L <- ncol(aln)
  n <- nrow(p)
  
  old_par <- par(no.readonly = TRUE)
  
  circos.clear()
  on.exit({
    circos.clear()
    layout(1)
    par(old_par)
  }, add = TRUE)
  
  # Left: circular plot. Right: comparison key and legend.
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.25, 1))
  par(
    mar = c(2, 1, 2, 1),
    oma = c(1, 0, 3, 0)
  )
  
  circos.par(
    start.degree = 90,
    gap.after = 70,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.003, 0.003),
    canvas.xlim = c(-1.12, 1.12),
    canvas.ylim = c(-1.12, 1.12),
    points.overflow.warning = FALSE
  )
  
  circos.initialize("Alignment", xlim = c(0, L))
  
  # Outer alignment-coordinate axis.
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = 0.06,
    bg.border = NA,
    panel.fun = function(x, y) {
      columns <- sort(unique(c(
        1,
        if (L >= 200) seq(200, L, by = 200),
        L
      )))
      
      circos.axis(
        h = "top",
        major.at = columns - 0.5,
        labels = columns,
        labels.cex = 0.65,
        minor.ticks = 0
      )
    }
  )
  
  ring_height <- min(0.12, 0.60 / n)
  
  for (i in seq_len(n)) {
    states <- classify_sites(p$A[i], p$B[i])
    
    blocks <- rle(states)
    ends <- cumsum(blocks$lengths)
    starts <- ends - blocks$lengths
    
    circos.trackPlotRegion(
      ylim = c(0, 1),
      track.height = ring_height,
      bg.border = NA,
      panel.fun = function(x, y) {
        circos.rect(
          xleft = starts,
          ybottom = 0,
          xright = ends,
          ytop = 1,
          col = unname(site_colors[blocks$values]),
          border = NA
        )
        
        # Outline each ring without outlining individual sites.
        circos.rect(
          xleft = 0,
          ybottom = 0,
          xright = L,
          ytop = 1,
          col = NA,
          border = "#CCCCCC",
          lwd = 0.4
        )
        
        # Short, horizontal ring number just inside the opening.
        circos.text(
          x = L + L * 8 / 290,
          y = 0.5,
          labels = i,
          facing = "downward",
          adj = c(0.5, 0.5),
          cex = 0.7,
          font = 2
        )
      }
    )
  }
  
  mtext(
    "Coordinates are alignment columns",
    side = 1,
    line = 0.4,
    cex = 0.8
  )
  
  # Separate panel for comparison labels and legend.
  par(
    mar = c(1, 0, 2, 1),
    cex = 1,
    mex = 1,
    pty = "m",
    xpd = FALSE
  )
  
  plot.new()
  plot.window(
    xlim = c(0, 1),
    ylim = c(0, 1),
    xaxs = "i",
    yaxs = "i"
  )
  
  text(
    0, 0.98,
    "Comparisons: outermost to innermost",
    adj = c(0, 1),
    font = 2,
    cex = 0.9
  )
  
  # Column headings.
  text(
    c(0, 0.07, 0.54),
    rep(0.91, 3),
    c("Ring", "A", "B"),
    adj = c(0, 0.5),
    font = 2,
    cex = 0.75
  )
  
  # Wrap longer sequence labels onto two lines.
  wrap_label <- function(z) {
    vapply(
      z,
      function(s) paste(strwrap(s, width = 25), collapse = "\n"),
      character(1)
    )
  }
  
  row_y <- 0.86 - (seq_len(n) - 1) * 0.038
  
  text(
    rep(0, n), row_y,
    labels = seq_len(n),
    adj = c(0, 0.5),
    font = 2,
    cex = 0.7
  )
  
  text(
    rep(0.07, n), row_y,
    labels = wrap_label(p$label_A),
    adj = c(0, 0.5),
    cex = 0.7
  )
  
  text(
    rep(0.54, n), row_y,
    labels = wrap_label(p$label_B),
    adj = c(0, 0.5),
    cex = 0.7
  )
  
  # Fixed-position legend below the comparison table.
  text(
    0, 0.28,
    "Site classification (A / B)",
    adj = c(0, 0.5),
    font = 2,
    cex = 0.85
  )
  
  legend_y <- 0.235 - (seq_along(site_colors) - 1) * 0.038
  
  rect(
    xleft = 0,
    ybottom = legend_y - 0.012,
    xright = 0.035,
    ytop = legend_y + 0.012,
    col = unname(site_colors),
    border = "#777777"
  )
  
  text(
    rep(0.055, length(site_colors)),
    legend_y,
    labels = names(site_colors),
    adj = c(0, 0.5),
    cex = 0.8
  )
  
  mtext(
    title,
    side = 3,
    outer = TRUE,
    line = 1,
    font = 2,
    cex = 1.2
  )
}

#plot within species plots
# narwhal
narwhal_pairs <- subset(
  pairs,
  within & species_A == "Narwhal"
)

plot_differences(
  narwhal_pairs,
  "Narwhal: parent and retrogene differences"
)
# beluga
plot_differences(
  subset(pairs, within & species_A == "Beluga"),
  "Beluga: parent and retrogene differences"
)

# all 12 within species comparisons together
plot_differences(
  subset(pairs, within),
  "Within-species comparisons"
)
# plot between species, starting with test: beluga and narwhal
parent_ids <- meta$Identifier[meta$Type == "Parent"]

selected <- subset(
  pairs,
  A %in% parent_ids &
    B %in% parent_ids &
    (
      (species_A == "Beluga" & species_B == "Narwhal") |
        (species_A == "Narwhal" & species_B == "Beluga")
    )
)

plot_differences(
  selected,
  "Beluga and narwhal: parent comparison"
)
# all nine cross species pairs for test
selected <- subset(
  pairs,
  (species_A == "Beluga" & species_B == "Narwhal") |
    (species_A == "Narwhal" & species_B == "Beluga")
)

plot_differences(
  selected,
  "Beluga versus narwhal: all sequence pairs"
)