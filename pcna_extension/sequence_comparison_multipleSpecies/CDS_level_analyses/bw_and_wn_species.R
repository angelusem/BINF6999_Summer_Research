library(ape)
library(ggplot2)

# Read your existing alignment.
file <- paste0(
  "C:/Users/angel/Downloads/Research project summer/",
  "pcnaRG_plus_parents_ncbinarwhal_all_hippo_cow_aligned.fasta"
)

aln <- read.dna(
  file, format = "fasta",
  as.character = TRUE, as.matrix = TRUE
)

# Standardize bases and remove descriptions after FASTA identifiers.
aln[] <- toupper(aln)
rownames(aln) <- sub("\\s.*$", "", rownames(aln))
stopifnot(!anyDuplicated(rownames(aln)))

# Specify the nine comparisons explicitly.
pairs <- data.frame(
  Species = c(
    "Bowhead whale", "North Atlantic right whale", "Blue whale",
    "Beluga", "Beluga", "Narwhal", "Narwhal",
    "Bottlenose dolphin", "Bottlenose dolphin"
  ),
  Group = c(rep("Baleen", 3), rep("Toothed", 6)),
  Copy = c("RG", "RG", "RG", "RG1", "RG2",
           "Chr5", "Chr7", "RG1", "RG2"),
  Retrogene = c(
    "BowheadPCNARetrogene",
    "NArightW_B_PCNA_like_core",
    "BlueWhale_retrogene_core",
    "Beluga_retrogene1_core",
    "Beluga_retrogene2_core",
    "CM143180.1:c95429311-95428529",
    "CM143182.1:1488428-1489213",
    "Bottlenose_PCNA_pseudogene1_core",
    "Bottlenose_PCNA_pseudogene2_core"
  ),
  Parent = c(
    "Bowhead_parent_PCNA_CDS_bmy_16007",
    "NARightW_parent_PCNA_CDS",
    "BlueWhale_parent_PCNA_CDS",
    "Beluga_parent_PCNA_CDS",
    "Beluga_parent_PCNA_CDS",
    "Narwhal_parent_PCNA_CDS",
    "Narwhal_parent_PCNA_CDS",
    "dolphin_parent_PCNA_CDS",
    "dolphin_parent_PCNA_CDS"
  )
)

# Select only the whale sequences needed for these comparisons.
ids <- unique(c(pairs$Retrogene, pairs$Parent))
missing <- setdiff(ids, rownames(aln))
if (length(missing) > 0) {
  stop("Missing sequence IDs: ", paste(missing, collapse = ", "))
}
whales <- aln[ids, , drop = FALSE]

# Keep the SAME A/C/G/T-only columns for all comparisons.
# Gaps and ambiguous bases are excluded.
keep <- apply(whales, 2, function(x) all(x %in% c("A", "C", "G", "T")))
n_sites <- sum(keep)
stopifnot(n_sites > 0)
message("Positions compared for every pair: ", n_sites)

# Calculate matching bases and percent identity.
pairs$Matches <- vapply(seq_len(nrow(pairs)), function(i) {
  sum(
    whales[pairs$Retrogene[i], keep] ==
      whales[pairs$Parent[i], keep]
  )
}, integer(1))

pairs$Sites_compared <- n_sites
pairs$Mismatches <- n_sites - pairs$Matches
pairs$Identity <- 100 * pairs$Matches / n_sites

print(pairs[, c("Species", "Copy", "Matches",
                "Sites_compared", "Identity")], row.names = FALSE)

# Order the species and separate overlapping copies slightly.
species_order <- unique(pairs$Species)
pairs$Species <- factor(pairs$Species, levels = rev(species_order))
pairs$Group <- factor(pairs$Group, levels = c("Baleen", "Toothed"))

dodge <- position_dodge(width = 0.45)

p <- ggplot(pairs, aes(
  x = Species, y = Identity,
  color = Group, shape = Copy, group = Copy
)) +
  geom_point(size = 3.5, position = dodge) +
  geom_text(
    aes(label = sprintf("%.2f%%", Identity)),
    position = dodge, hjust = -0.2,
    size = 3.3, show.legend = FALSE
  ) +
  coord_flip() +
  scale_shape_manual(values = c(
    RG1 = 16, RG2 = 17, Chr5 = 15, Chr7 = 18
  )) +
  scale_color_manual(values = c(
    Baleen = "#0072B2", Toothed = "#D55E00"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.30))) +
  labs(
    title = "PCNA retrogene similarity to its own parent",
    subtitle = paste(n_sites, "shared ungapped alignment positions"),
    x = NULL,
    y = "DNA identity to same-species parent (%)",
    caption = "Gaps excluded. Narwhal copies labeled by chromosome."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom"
  )

print(p)

# Save the table and figure in a new folder beside your FASTA.
out <- file.path(dirname(file), "PCNA_identity_results")
dir.create(out, showWarnings = FALSE)

write.csv(pairs, file.path(out, "parent_retrogene_identity.csv"),
          row.names = FALSE)
ggsave(file.path(out, "parent_retrogene_identity.png"),
       plot = p, width = 9, height = 5, dpi = 300)

message("Results saved in: ", out)
# Change baleen labels to RG.
pairs$Copy <- as.character(pairs$Copy)
pairs$Copy[pairs$Group == "Baleen"] <- "RG"

# Build labels, keeping each species' parent and retrogenes together.
meta <- do.call(rbind, lapply(species_order, function(sp) {
  rows <- pairs[as.character(pairs$Species) == sp, ]
  
  data.frame(
    ID = c(rows$Parent[1], rows$Retrogene),
    Species = sp,
    Copy = c("Parent", rows$Copy),
    stringsAsFactors = FALSE
  )
}))

meta$Label <- paste(meta$Species, meta$Copy, sep = " | ")

# Select the same shared alignment positions as before.
seqs <- whales[meta$ID, keep, drop = FALSE]

# Calculate percent differences for all sequence pairs.
n <- nrow(seqs)
difference <- matrix(
  0, nrow = n, ncol = n,
  dimnames = list(meta$Label, meta$Label)
)

for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    difference[i, j] <- 100 * mean(seqs[i, ] != seqs[j, ])
  }
}

# Convert the matrix to a plotting table.
heat <- as.data.frame(as.table(difference))
names(heat) <- c("Row", "Column", "Difference")

heat$Column <- factor(heat$Column, levels = meta$Label)
heat$Row <- factor(heat$Row, levels = rev(meta$Label))

# Outline cells comparing sequences from the same species.
heat$Within_species <-
  meta$Species[match(as.character(heat$Row), meta$Label)] ==
  meta$Species[match(as.character(heat$Column), meta$Label)]

p_heat <- ggplot(heat, aes(Column, Row, fill = Difference)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_tile(
    data = heat[heat$Within_species, ],
    fill = NA, color = "grey30", linewidth = 0.5
  ) +
  geom_text(aes(label = sprintf("%.1f", Difference)), size = 2.5) +
  scale_fill_gradient(
    low = "#FFF7EC", high = "#B30000",
    limits = c(0, max(heat$Difference)),
    name = "Different\nbases (%)"
  ) +
  coord_fixed() +
  labs(
    title = "PCNA sequence differences within and between species",
    subtitle = paste(sum(keep), "shared ungapped positions"),
    x = NULL, y = NULL,
    caption = paste(
      "Outlined cells: within-species comparisons.",
      "Baleen copies: RG. Narwhal copies: Chr5 and Chr7."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1),
    axis.text.y = element_text(size = 9)
  )

print(p_heat)

ggsave(
  file.path(out, "within_between_species_differences.png"),
  p_heat, width = 12, height = 10, dpi = 300
)

write.csv(
  difference,
  file.path(out, "pairwise_percent_differences.csv")
)
#####
# Compare every retrogene with its own species' parent.
pairs$Copy <- as.character(pairs$Copy)
pairs$Copy[pairs$Group == "Baleen"] <- "RG"

positions <- which(keep)

changes <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  parent <- as.character(whales[pairs$Parent[i], keep])
  rg <- as.character(whales[pairs$Retrogene[i], keep])
  
  data.frame(
    Species = as.character(pairs$Species[i]),
    Copy = pairs$Copy[i],
    Position = positions,
    Parent_base = parent,
    RG_base = rg,
    Changed = parent != rg,
    Change = ifelse(parent == rg, "Match",
                    paste0(parent, ">", rg)),
    stringsAsFactors = FALSE
  )
}))

changes$Label <- paste(changes$Species, changes$Copy, sep = " | ")

# Count DISTINCT SPECIES sharing an exact change at the same position.
# Two copies in one species therefore don't count as two species.
events <- unique(
  changes[changes$Changed, c("Position", "Change", "Species")]
)

if (nrow(events) == 0) {
  stop("No substitutions found in the selected shared positions.")
}

shared <- aggregate(
  Species ~ Position + Change,
  data = events,
  FUN = function(x) length(unique(x))
)
names(shared)[3] <- "Species_count"

shared <- shared[shared$Species_count >= 2, ]
shared <- shared[order(shared$Position, shared$Change), ]

print(shared, row.names = FALSE)

write.csv(
  shared,
  file.path(out, "changes_shared_across_species.csv"),
  row.names = FALSE
)

# Plot positions containing a change shared across species.
if (nrow(shared) == 0) {
  message("No exact substitutions are shared across two or more species.")
} else {
  selected <- changes[
    changes$Position %in% shared$Position, ]
  
  selected$Position <- factor(
    selected$Position,
    levels = sort(unique(selected$Position))
  )
  selected$Label <- factor(
    selected$Label,
    levels = rev(unique(changes$Label))
  )
  
  # Color by retrogene base when different; matching cells stay grey.
  selected$Display <- ifelse(
    selected$Changed, selected$RG_base, "Match"
  )
  
  p_shared <- ggplot(
    selected, aes(Position, Label, fill = Display)
  ) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(
      aes(label = ifelse(Changed, Change, "\u00b7")),
      size = 2.5
    ) +
    scale_fill_manual(
      values = c(
        Match = "#EEEEEE",
        A = "#8DD3C7", C = "#80B1D3",
        G = "#FDB462", T = "#FB8072"
      ),
      name = "RG base\nwhen different"
    ) +
    labs(
      title = "Shared retrogene–parent differences",
      subtitle = "Positions with an exact change shared by at least two species",
      x = "Position in the original alignment",
      y = NULL,
      caption = "Letters show parent > RG; dots indicate a match to the own-species parent."
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
  
  print(p_shared)
  p_shared
  
  ggsave(
    file.path(out, "shared_retrogene_changes.png"),
    p_shared,
    width = max(10, 0.4 * length(unique(selected$Position))),
    height = 5.5,
    dpi = 300,
    limitsize = FALSE
  )
}