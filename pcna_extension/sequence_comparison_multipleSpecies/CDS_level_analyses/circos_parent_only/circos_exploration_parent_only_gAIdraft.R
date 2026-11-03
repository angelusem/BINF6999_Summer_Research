library(ape)
library(circlize)

alignment <- read.dna(
  "circos_parent_only/aligned_PCNA_parent_only_sep30.fasta",
  format = "fasta"
)

divergence_pct <- 100 * as.matrix(
  dist.dna(alignment, model = "raw", pairwise.deletion = TRUE)
)
diag(divergence_pct) <- 0

# Count alignment positions containing A, C, G, or T in both sequences
bases <- toupper(as.character(alignment))
valid <- matrix(
  bases %in% c("A", "C", "G", "T"),
  nrow = nrow(bases),
  dimnames = dimnames(bases)
)

comparable_sites <- tcrossprod(valid * 1)

# 1. Read your MAFFT FASTA alignment file
# 'format = "fasta"' reads sequential/interleaved FASTA alignment data
alignment <- read.dna("circos_parent_only/aligned_PCNA_parent_only_sep30.fasta", format = "fasta")

# 2. Compute the pairwise distance matrix
# method = "raw" computes the proportion of sites that differ (0 to 1)
raw_dist <- dist.dna(alignment, model = "raw", pairwise.deletion = TRUE)

# 3. Convert distance to a Percentage Similarity Matrix
# Similarity % = (1 - distance) * 100
sim_matrix <- as.matrix(1 - raw_dist) * 100

# 4. Clean data for the Chord Diagram
# Set the diagonal (self-similarity of 100%) to 0 so species don't link to themselves
diag(sim_matrix) <- 0

# Optional: Replace any missing alignment values (NaNs from heavy gapping) with 0
sim_matrix[is.nan(sim_matrix)] <- 0

# Filter: Set a threshold so only strong alignments show up, greater than 70%
#sim_matrix_filtered <- sim_matrix
#sim_matrix_filtered[sim_matrix_filtered < 70] <- 0

#DIFFERENCES  --------------------------------------


# 1. Force R to completely wipe its active graphic device memory
if(!is.null(dev.list())) dev.off()

# 2. establish filtered divergence data structure
divergence_matrix <- 100 - sim_matrix
diag(divergence_matrix) <- 0

# Set a strict divergence range (e.g., only show mutations between 1% and 15%)
divergence_matrix[divergence_matrix < 1 | divergence_matrix > 15] <- 0

# 3. Configure a freshly cleared circular layout grid
circos.clear()
circos.par(
  gap.after = 2, 
  start.degree = 90, 
  canvas.xlim = c(-1.3, 1.3), # Create a wider buffer zone to contain labels
  canvas.ylim = c(-1.3, 1.3)
)

# 4. Map colors to your sectors
species_names <- rownames(divergence_matrix)
grid_colors <- setNames(rainbow(length(species_names), v = 0.75, s = 0.6), species_names)

# 5. Render the Symmetric Diagram with strict link styling
chordDiagram(
  x = divergence_matrix,
  symmetric = TRUE,               # Keep links as flat uniform bands
  link.sort = TRUE,               
  link.decreasing = FALSE,        
  grid.col = grid_colors,
  transparency = 0.4,
  annotationTrack = "grid",       # ONLY draw the grid rings here
  preAllocateTracks = list(track.height = 0.15) # Leave a pristine track empty for labels
)

# 6. Explicitly draw text labels into the allocated empty outer track
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(
    x = CELL_META$xcenter, 
    y = CELL_META$ylim[1] + 0.05,  # Force text just outside the color arc
    labels = CELL_META$sector.index, 
    facing = "clockwise", 
    niceFacing = TRUE, 
    adj = c(0, 0.5),                
    cex = 0.55                     # Scales down font size to cleanly fit
  )
}, bg.border = NA)

# 7. Close out active circlize background processes
circos.clear()
###
# Wipe active graphics
if(!is.null(dev.list())) dev.off()

# Use similarity directly
similarity_matrix <- sim_matrix
diag(similarity_matrix) <- 0

# FILTER: Only show links between 85% and 99% similarity 
# (Equivalent to your 1% to 15% divergence goal)
similarity_matrix[similarity_matrix < 85 | similarity_matrix > 99] <- 0

circos.clear()
circos.par(gap.after = 2, start.degree = 90, canvas.xlim = c(-1.3, 1.3), canvas.ylim = c(-1.3, 1.3))
grid_colors <- setNames(rainbow(length(rownames(similarity_matrix)), v = 0.75, s = 0.6), rownames(similarity_matrix))

chordDiagram(
  x = similarity_matrix,
  symmetric = TRUE,               
  grid.col = grid_colors,
  transparency = 0.4,
  annotationTrack = "grid",       
  preAllocateTracks = list(track.height = 0.15) 
)
# ... [Keep your circos.track text labeling code here] ...
####
# 1. Clear previous graphic devices and active configurations
if(!is.null(dev.list())) dev.off()
circos.clear()

# 2. Prepare pure similarity data matrix
similarity_matrix <- sim_matrix
diag(similarity_matrix) <- 0  # Remove self-links

# 3. Dynamic Filtering: Automatically find a threshold that retains variation
# This filters out the bottom 25% of background noise, keeping the strongest relative patterns
min_cutoff <- quantile(similarity_matrix[similarity_matrix > 0], 0.25)
similarity_matrix[similarity_matrix < min_cutoff] <- 0

# 4. Set Canvas Buffers for Clear Labels
circos.par(
  gap.after = 3,            # Slight gap increase to isolate sectors
  start.degree = 90, 
  canvas.xlim = c(-1.4, 1.4), 
  canvas.ylim = c(-1.4, 1.4)
)

# 5. Define Consistent Colors
species_names <- rownames(similarity_matrix)
grid_colors <- setNames(rainbow(length(species_names), v = 0.75, s = 0.6), species_names)

# 6. Render the undirected chord graph
chordDiagram(
  x = similarity_matrix,
  symmetric = TRUE,
  directional = 0,               # 0 avoids the directional bias distortion
  link.sort = TRUE,               
  link.decreasing = TRUE,        # Plots strongest relationships on top
  grid.col = grid_colors,
  transparency = 0.4,
  annotationTrack = "grid",       
  preAllocateTracks = list(track.height = 0.18) # Generous text spacing
)

# 7. Add readable text labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(
    x = CELL_META$xcenter, 
    y = CELL_META$ylim[1] + 0.08,  
    labels = CELL_META$sector.index, 
    facing = "clockwise", 
    niceFacing = TRUE, 
    adj = c(0, 0.5),                
    cex = 0.55                     
  )
}, bg.border = NA)

circos.clear()
