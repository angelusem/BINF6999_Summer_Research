# install.packages(c("circlize")) 
# BiocManager::install(c("Biostrings", "DECIPHER"))
library(Biostrings)
library(DECIPHER)
library(circlize)

# 1. Input paths and settings -----------------------------------------------
input_dir <- "C:/Users/angel/Downloads/Research project summer"
seqs <- readDNAStringSet(file.path(input_dir, "pcnaRG_plus_parents_all.fasta"))
meta <- read.csv(file.path(input_dir, "retrogene_parent_pcna_metadata.csv"),
                 fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
reference <- "Beluga_parent_PCNA_CDS"  # Shared coordinates/orientation only; NOT the identity baseline
window_size <- 60                     # Smaller = more positional detail
output_dir <- file.path(getwd(), "circos_own_parent_output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Validate, orient, and align -------------------------------------------
stopifnot(!anyDuplicated(names(seqs)), !anyDuplicated(meta$Identifier),
          setequal(names(seqs), meta$Identifier), reference %in% names(seqs))
stopifnot(all(table(meta$Species[meta$Type == "Parent"]) == 1),
          setequal(meta$Species, meta$Species[meta$Type == "Parent"]))
oriented <- OrientNucleotides(seqs, reference = match(reference, names(seqs)),
                              processors = 1)
alignment <- AlignSeqs(oriented, processors = 1)
writeXStringSet(alignment, file.path(output_dir, "full_alignment.fasta"))
m <- as.matrix(alignment)
ref <- m[reference, ]
# Keep reference nucleotide positions only: excludes flanks and insertions
# relative to this reference. Full alignment above retains those columns.
m <- m[, ref != "-", drop = FALSE]
ref <- m[reference, ]
n <- ncol(m)

# 3. Order species, then parent / RG1 / RG2 ---------------------------------
species_order <- c("Bowhead whale", "North Atlantic right whale", "Blue whale",
                   "Bottlenose dolphin", "Beluga", "Narwhal")
meta <- meta[order(match(meta$Species, species_order),
                   match(meta$Copy, c("Parent", "RG1", "RG2"))), ]
ids <- meta$Identifier
species_colors <- setNames(c("#8C6BB1", "#A86432", "#596BB3",
                             "#C65067", "#218C74", "#247BA0"), species_order)
short <- c("Bowhead whale" = "Bowhead", "North Atlantic right whale" = "NA right",
           "Blue whale" = "Blue whale", "Bottlenose dolphin" = "Dolphin",
           "Beluga" = "Beluga", "Narwhal" = "Narwhal")
labels <- setNames(paste(short[meta$Species], meta$Copy), ids)
labels <- setNames(
  paste(short[meta$Species],
        ifelse(meta$Type == "Parent", "P", meta$Copy)),
  ids
)
identity_color <- colorRamp2(c(80, 90, 100), c("#B34D73", "#F0CF84", "#247BA0"))

# 4. Compare each sequence with its OWN species' parent --------------------
# Beluga remains only the shared coordinate/orientation reference.
# Parent gaps or ambiguous bases are excluded, not scored as differences.
# Internal retrogene gaps count as differences; uncovered ends are excluded.
parents <- meta[meta$Type == "Parent", ]
parent_by_species <- setNames(parents$Identifier, parents$Species)

stats <- do.call(rbind, lapply(ids, function(id) {
  species <- meta$Species[match(id, meta$Identifier)]
  parent_id <- unname(parent_by_species[species])
  s <- m[id, ]
  parent_seq <- m[parent_id, ]
  covered <- which(s != "-")
  stopifnot(length(covered) > 0)
  
  do.call(rbind, lapply(seq(1, n, by = window_size), function(start) {
    end <- min(start + window_size - 1, n)
    p <- start:end
    usable <- p >= min(covered) & p <= max(covered) &
      s[p] %in% c("A", "C", "G", "T", "-") &
      parent_seq[p] %in% c("A", "C", "G", "T")
    
    identity <- if (mean(usable) >= .8) {
      100 * mean(s[p[usable]] == parent_seq[p[usable]])
    } else NA_real_
    
    data.frame(id, parent_id, start, end,
               usable_sites = sum(usable),
               coverage = mean(usable), identity)
  }))
}))
write.csv(stats, file.path(output_dir, "window_identity.csv"), row.names = FALSE)

# 5. Draw the circular plot ------------------------------------------------
draw_circos <- function() {
  circos.clear()
  on.exit(circos.clear())
  par(mar = c(5, 2, 4, 2), xpd = NA)
  
  last_in_species <- !duplicated(meta$Species, fromLast = TRUE)
  
  circos.par(
    start.degree = 90,
    gap.after = ifelse(last_in_species, 8, 2),
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(.008, .008),
    canvas.xlim = c(-1.25, 1.25),
    canvas.ylim = c(-1.25, 1.25)
  )
  circos.initialize(factor(ids, levels = ids), xlim = c(0, n))
  
  # Outer species strip with outward-facing labels
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = .06,
    bg.border = NA,
    panel.fun = function(x, y) {
      
      id <- CELL_META$sector.index
      species <- meta$Species[match(id, ids)]
      
      circos.rect(
        0, 0, n, 1,
        col = species_colors[species],
        border = NA
      )
      
      circos.text(
        CELL_META$xcenter,
        CELL_META$ylim[2] + mm_y(2),
        labels[id],
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = .65
      )
    }
  )
  
  # Identity heatmap, with an axis for EVERY sector
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = .10,
    bg.border = NA,
    panel.fun = function(x, y) {
      
      id <- CELL_META$sector.index
      z <- stats[stats$id == id, ]
      
      colors <- rep("#CCD2DA", nrow(z))
      ok <- !is.na(z$identity)
      colors[ok] <- identity_color(z$identity[ok])
      
      circos.rect(
        z$start - 1, 0, z$end, 1,
        col = colors,
        border = NA
      )
      
      circos.axis(
        h = "bottom",
        major.at = c(0, 400, n),
        labels = c("1", "400", as.character(n)),
        labels.cex = .45,
        labels.facing = "clockwise",
        major.tick.length = .025
      )
    }
  )
  # Color = identity for the surrounding window, matching the heatmap.
  # Split ribbons at gaps/ambiguous bases so they do not bridge missing data.
  
  for (i in which(meta$Type == "Retrogene")) {
    
    retro_id <- meta$Identifier[i]
    parent_id <- unname(parent_by_species[meta$Species[i]])
    
    retro_seq <- m[retro_id, ]
    parent_seq <- m[parent_id, ]
    
    # Use only windows that passed the existing coverage threshold
    windows <- stats[
      stats$id == retro_id & !is.na(stats$identity),
      ,
      drop = FALSE
    ]
    
    for (j in seq_len(nrow(windows))) {
      
      positions <- windows$start[j]:windows$end[j]
      
      # Draw only where both sequences contain an actual nucleotide.
      # Mismatching nucleotides are still included.
      valid <- positions[
        retro_seq[positions] %in% c("A", "C", "G", "T") &
          parent_seq[positions] %in% c("A", "C", "G", "T")
      ]
      
      if (length(valid) == 0L) next
      
      # Separate contiguous stretches at gaps or ambiguous positions
      stretches <- split(
        valid,
        cumsum(c(TRUE, diff(valid) != 1L))
      )
      
      ribbon_color <- adjustcolor(
        identity_color(windows$identity[j]),
        alpha.f = 0.4
      )
      
      for (stretch in stretches) {
        
        # Base 1 occupies plotting interval 0–1, etc.
        interval <- c(min(stretch) - 1, max(stretch))
        
        circos.link(
          sector.index1 = parent_id,
          point1 = interval,
          sector.index2 = retro_id,
          point2 = interval,
          col = ribbon_color,
          border = NA,
          directional = 0
        )
      }
    }
  }
  
  title("PCNA parents and retrogenes", cex.main = 1.1)
  
  legend(
    x = 0,
    y = -1.13,
    xjust = 0.5,
    yjust = 1,
    bty = "n",
    horiz = TRUE,
    cex = 0.65,
    legend = c("<=80%", "90%", "100%", "<80% usable sites"),
    fill = c(identity_color(c(80, 90, 100)), "#CCD2DA"),
    border = NA,
    title = "Window identity to own species parent"
  )
  
} # End draw_circos()
draw_circos() # RStudio Plots pane
png(file.path(output_dir, "exploratory_circos.png"), width = 2400, height = 2400, res = 240)
tryCatch(draw_circos(), finally = dev.off())
message("Saved files in: ", normalizePath(output_dir))