library(ape)
library(ggplot2)

input_file <- "C:/Users/angel/Downloads/Research project summer/exploratory_analyses_15_seqs/Aligned_extended_pcnaAndRG_pNArightwhale_gray_minke_orca.fasta"

mapping <- 'Species|Group|Copy|Retrogene|Parent
Bowhead whale|Baleen|RG|BowheadPCNARetrogene|Bowhead_parent_full_bmy16007
North Atlantic right whale|Baleen|RG|NArightW_B_PCNA_like_core|NARightW_parent_PCNA_CDS
Blue whale|Baleen|RG|BlueWhale_retrogene_core|BlueWhale_parent_PCNA_CDS
Humpback whale|Baleen|RG|Humpback_PCNA_RG|Humpback_PCNA_parent
Gray whale|Baleen|RG|Graywhale_PCNA_RG_chr5|Graywhale_PCNA_parent_GeneID_137750119
Minke whale|Baleen|RG|Minke_PCNA_RG_chr8|Minke_PCNA_parent_GeneID_102999134
Beluga|Toothed|RG1|Beluga_retrogene1_core|Beluga_parent_PCNA_CDS
Beluga|Toothed|RG2|Beluga_retrogene2_core|Beluga_parent_PCNA_CDS
Narwhal|Toothed|RG chr5|Monodon_monoceros_chr5_isolate_NGI_ecotype_Baffin_Bay_Population|Narwhal_parent_PCNA_CDS
Narwhal|Toothed|RG chr7|Monodon_monoceros_chr7_isolate_NGI_ecotype_Baffin_Bay_Population|Narwhal_parent_PCNA_CDS
Bottlenose dolphin|Toothed|RG1|Bottlenose_PCNA_pseudogene1_core|dolphin_parent_PCNA_CDS
Bottlenose dolphin|Toothed|RG2|Bottlenose_PCNA_pseudogene2_core|dolphin_parent_PCNA_CDS
Sperm whale|Toothed|RG|Spermwhale_PCNA_like_GeneID_102980166|Spermwhale_PCNA_parent_GeneID_102992318
Pygmy sperm whale|Toothed|RG chr2|PygmySpwhale_PCNA_RG_chr2|PygmySpwhale_PCNA_parent_GeneID_131741109
Orca|Toothed|RG chr4|Orca_PCNA_RG_chr4|Orca_PCNA_parent_GeneID_101287992
Orca|Toothed|RG chr7|Orca_PCNA_RG_ch7|Orca_PCNA_parent_GeneID_101287992'
# Read the alignment and your sequence-pairing table.
aln <- read.dna(input_file, format = "fasta",
                as.character = TRUE, as.matrix = TRUE)
aln[] <- toupper(aln)
rownames(aln) <- sub("\\s.*$", "", rownames(aln))

pairs <- read.delim(text = mapping, sep = "|",
                    stringsAsFactors = FALSE)

# Keep only the sequences listed in the mapping.
ids <- unique(c(pairs$Parent, pairs$Retrogene))
stopifnot(all(ids %in% rownames(aln)))
aln <- aln[ids, , drop = FALSE]

# Keep positions with A/C/G/T in every included sequence.
keep <- apply(aln, 2, function(x) all(x %in% c("A", "C", "G", "T")))
stopifnot(any(keep))

# Compare each RG with its own parent.
d <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  parent <- aln[pairs$Parent[i], keep]
  rg <- aln[pairs$Retrogene[i], keep]
  
  data.frame(
    Species = pairs$Species[i],
    Group = pairs$Group[i],
    Label = paste(pairs$Species[i], pairs$Copy[i], sep = " | "),
    Position = which(keep),
    Change = ifelse(parent == rg, ".", paste0(parent, ">", rg)),
    Base = ifelse(parent == rg, "Match", rg)
  )
}))

row_order <- unique(d$Label)

# reordering within the toothed-whale panel.

toothed_order <- c(
  "Sperm whale | RG candidate",
  "Pygmy sperm whale | RG chr2",
  
  # Remaining RG2 copies
  "Beluga | RG2",
  "Narwhal | RG chr7",
  "Bottlenose dolphin | RG2",
  "Orca | RG chr7",
  
  # RG1 copies
  "Beluga | RG1",
  "Narwhal | RG chr5",
  "Bottlenose dolphin | RG1",
  "Orca | RG chr4"
)

# Keep the existing baleen-whale order.
baleen_order <- unique(d$Label[d$Group == "Baleen"])
row_order <- c(baleen_order, toothed_order)

# Check that every row is included exactly once.
stopifnot(
  !anyDuplicated(row_order),
  setequal(row_order, unique(d$Label))
)
# Find exact differences shared by at least two distinct species.
events <- unique(d[d$Change != ".", c("Position", "Change", "Species")])
if (!nrow(events)) stop("No differences found.")

counts <- aggregate(Species ~ Position + Change, events,
                    function(x) length(unique(x)))

sites <- sort(unique(counts$Position[counts$Species >= 2]))
if (!length(sites)) stop("No differences shared across two species.")

# Show every RG at those selected positions.
d <- d[d$Position %in% sites, ]
d$Position <- factor(d$Position, levels = sites)
d$Label <- factor(d$Label, levels = rev(row_order))
d$Group <- factor(d$Group, levels = c("Baleen", "Toothed"))

p <- ggplot(d, aes(Position, Label, fill = Base)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Change), size = 2.5) +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(
    Match = "#EEEEEE", A = "#8DD3C7", C = "#80B1D3",
    G = "#FDB462", T = "#FB8072"
  )) +
  labs(
    title = "Shared parent–retrogene differences",
    subtitle = paste(sum(keep), "shared positions screened"),
    x = "Alignment position", y = NULL,
    caption = "Text: parent > RG. Dot: matches own parent."
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "bottom"
  )

print(p)

ggsave(
  file.path(dirname(input_file), "shared_differences_single_plot.png"),
  p, width = max(12, length(sites) * 0.35 + 4),
  height = 8, dpi = 300
)

### within species comparisons
# Find species with more than one retrogene.
by_species <- split(pairs, pairs$Species)
multiple <- by_species[sapply(by_species, nrow) > 1]

# Compare all within-species RG pairs.
within <- do.call(rbind, lapply(multiple, function(sp) {
  combinations <- combn(seq_len(nrow(sp)), 2)
  
  do.call(rbind, lapply(seq_len(ncol(combinations)), function(k) {
    i <- combinations[1, k]
    j <- combinations[2, k]
    
    rg_a <- aln[sp$Retrogene[i], keep]
    rg_b <- aln[sp$Retrogene[j], keep]
    
    differences <- sum(rg_a != rg_b)
    
    data.frame(
      Species = sp$Species[i],
      Comparison = paste(sp$Copy[i], "vs", sp$Copy[j]),
      Positions_compared = length(rg_a),
      Differences = differences,
      Percent_different = 100 * differences / length(rg_a)
    )
  }))
}))

print(within, row.names = FALSE)

# Simple dot plot.
p_within <- ggplot(
  within,
  aes(Percent_different, Species)
) +
  geom_point(size = 3, color = "#0072B2") +
  geom_text(
    aes(label = sprintf("%.2f%%", Percent_different)),
    hjust = -0.2
  ) +
  scale_x_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  labs(
    title = "Differences between retrogenes within each species",
    subtitle = paste(sum(keep), "shared ungapped positions compared"),
    x = "Different DNA bases (%)",
    y = NULL
  ) +
  theme_minimal()

print(p_within)

