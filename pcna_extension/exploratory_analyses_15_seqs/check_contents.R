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
reference <- "Beluga_parent_PCNA_CDS"  # Display reference as just a placeholder, NOT an ancestor
window_size <- 30                     # Smaller = more positional detail
output_dir <- file.path(getwd(), "circos_output")
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
identity_color <- colorRamp2(c(80, 90, 100), c("#B34D73", "#F0CF84", "#247BA0"))

# 4. Identity in windows, relative to the same Beluga parent for every row ---
# Internal gaps count as differences. Uncovered ends are excluded.
# Windows with <80% usable coverage are gray; insertions are not scored.
stats <- do.call(rbind, lapply(ids, function(id) {
  s <- m[id, ]
  covered <- which(s != "-")
  do.call(rbind, lapply(seq(1, n, by = window_size), function(start) {
    end <- min(start + window_size - 1, n)
    p <- start:end
    usable <- p >= min(covered) & p <= max(covered) &
      s[p] %in% c("A", "C", "G", "T", "-") &
      ref[p] %in% c("A", "C", "G", "T")
    identity <- if (mean(usable) >= .8) 100 * mean(s[p[usable]] == ref[p[usable]]) else NA_real_
    data.frame(id, start, end, coverage = mean(usable), identity)
  }))
}))
write.csv(stats, file.path(output_dir, "window_identity.csv"), row.names = FALSE)

# 5. Draw the circular plot ------------------------------------------------
draw_circos <- function() {
  circos.clear()
  on.exit(circos.clear())
  par(mar = c(3, 1, 4, 1))
  last_in_species <- !duplicated(meta$Species, fromLast = TRUE)
  circos.par(start.degree = 90, gap.after = ifelse(last_in_species, 8, 2),
             cell.padding = c(0, 0, 0, 0), track.margin = c(.008, .008))
  circos.initialize(factor(ids, levels = ids), xlim = c(0, n))
  
  # Outer strip identifies species and labels each sequence.
  circos.trackPlotRegion(ylim = c(0, 1), track.height = .12, bg.border = NA,
                         panel.fun = function(x, y) {
                           id <- CELL_META$sector.index
                           species <- meta$Species[match(id, ids)]
                           circos.rect(0, 0, n, .2, col = species_colors[species], border = NA)
                           circos.text(n / 2, .65, labels[id], facing = "bending.inside",
                                       niceFacing = TRUE, cex = .57)
                         })
  
  # Heatmap ring: local identity to the common reference.
  circos.trackPlotRegion(ylim = c(0, 1), track.height = .10, bg.border = NA,
                         panel.fun = function(x, y) {
                           z <- stats[stats$id == CELL_META$sector.index, ]
                           colors <- rep("#CCD2DA", nrow(z))
                           ok <- !is.na(z$identity)
                           colors[ok] <- identity_color(z$identity[ok])
                           circos.rect(z$start - 1, 0, z$end, 1, col = colors, border = NA)
                           circos.axis(h = "bottom", major.at = c(0, 200, 400, 600, n),
                                       labels = c(1, 200, 400, 600, n), labels.cex = .35,
                                       major.tick.length = .025)
                         })
  
  # Thin links show supplied within-species parent/retrogene associations.
  # They do NOT encode inferred evolutionary relationships or direction.
  for (i in which(meta$Type == "Retrogene")) {
    parent <- meta$Identifier[meta$Species == meta$Species[i] & meta$Type == "Parent"]
    circos.link(parent, n / 2, meta$Identifier[i], n / 2,
                col = adjustcolor(species_colors[meta$Species[i]], alpha.f = .5), lwd = 1.5)
  }
  title("PCNA parents and retrogenes", cex.main = 1.1)
  legend("bottom", inset = -.05, bty = "n", horiz = TRUE, cex = .65,
         legend = c("<=80%", "90%", "100%", "Low coverage"),
         fill = c(identity_color(c(80, 90, 100)), "#CCD2DA"), border = NA,
         title = "Window identity to Beluga parental CDS")
}

draw_circos() # RStudio Plots pane
png(file.path(output_dir, "exploratory_circos.png"), width = 2400, height = 2400, res = 240)
tryCatch(draw_circos(), finally = dev.off())
message("Saved files in: ", normalizePath(output_dir))
# 6. Compare each RG1 to corresponding RG2
# Full alignment, including insertion columns
full <- as.matrix(alignment)
valid <- c("A", "C", "G", "T")

# Preserve a mapping to Beluga CDS positions
ref_full <- full["Beluga_parent_PCNA_CDS", ]
ref_position <- cumsum(ref_full != "-")

two_copy_species <- c("Bottlenose dolphin", "Beluga", "Narwhal")

pairwise_results <- lapply(two_copy_species, function(sp) {
  
  id1 <- meta$Identifier[
    meta$Species == sp & meta$Copy == "RG1"
  ]
  id2 <- meta$Identifier[
    meta$Species == sp & meta$Copy == "RG2"
  ]
  stopifnot(length(id1) == 1, length(id2) == 1)
  
  a <- full[id1, ]
  b <- full[id2, ]
  
  # Compare only the region covered by both sequences
  span_a <- range(which(a != "-"))
  span_b <- range(which(b != "-"))
  
  pos <- seq_along(a)
  overlap <- pos >= max(span_a[1], span_b[1]) &
    pos <= min(span_a[2], span_b[2])
  
  status <- rep("Unavailable", length(a))
  
  bases <- overlap & a %in% valid & b %in% valid
  status[bases & a == b] <- "Match"
  status[bases & a != b] <- "Substitution"
  
  gap <- overlap & (
    (a == "-" & b %in% valid) |
      (b == "-" & a %in% valid)
  )
  status[gap] <- "Gap difference"
  
  data.frame(
    Species = sp,
    alignment_column = pos,
    reference_position = ifelse(ref_full == "-", NA, ref_position),
    RG1 = a,
    RG2 = b,
    status = status
  )
})

pairwise <- do.call(rbind, pairwise_results)

write.csv(
  pairwise,
  file.path(output_dir, "RG1_vs_RG2_positions.csv"),
  row.names = FALSE
)

# Counts of matching bases, substitutions, and gap differences
with(pairwise, table(Species, status))
