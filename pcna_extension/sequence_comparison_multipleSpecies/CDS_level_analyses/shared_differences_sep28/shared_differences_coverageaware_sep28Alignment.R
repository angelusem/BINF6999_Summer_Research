# Standalone shared parent/RG differences, including alignment gaps.
library(ape)
library(ggplot2)

input_file <- "C:/Users/angel/Downloads/Research project summer/exploratory_analyses_15_seqs/shared_differences_sep28/aligned_cetPCNA_RGandParent_sep28_bowheadFull.fasta"
output_dir <- file.path(dirname(input_file), "shared_differences_fullbowhead_sep28")
min_species <- 2L

# Only listed sequences are included. Baleen copies are labeled RG.
# Sperm-whale PCNA-like is a provisional RG candidate.
# All 28 sequences in the September 28 alignment are mapped below.
# Narwhal RG1/RG2 follow the new FASTA labels; chromosome identity is not inferred.
mapping <- 'Species|Group|Copy|Retrogene|Parent
Bowhead whale|Baleen|RG1|bmy_21945T0|Bowhead_parent
North Atlantic right whale|Baleen|RG1|NArightWhale_retrogene|NARightWhale_parent
Blue whale|Baleen|RG1|BlueWhale_retrogene|BlueWhale_parent
Humpback whale|Baleen|RG1|HumpbackWhale_retrogene|HumpbackWhale_parent
Gray whale|Baleen|RG1|GrayWhale_retrogene_chr5|GrayWhale_parent
Minke whale|Baleen|RG1|MinkeWhale_retrogene_chr8|MinkeWhale_parent
Beluga|Toothed|RG2|BelugaWhale_retrogene1|BelugaWhale_parent
Beluga|Toothed|RG1|BelugaWhale_retrogene2|BelugaWhale_parent
Narwhal|Toothed|RG2|Narwhal_retrogene1|Narwhal_parent
Narwhal|Toothed|RG1|Narwhal_retrogene2|Narwhal_parent
Bottlenose dolphin|Toothed|RG2|BottlenoseDolphin_retrogene1|BottlenosDolphin_parent
Bottlenose dolphin|Toothed|RG1|BottlenoseDolphin_retrogene2|BottlenosDolphin_parent
Sperm whale|Toothed|RG1|SpermWhale_retrogene|SpermWhale_parent
Pygmy sperm whale|Toothed|RG1|PygmySpermWhale_retrogene_chr2|PygmySpermWhale_parent
Orca|Toothed|RG2|Orca_retrogene_chr4|Orca_parent
Orca|Toothed|RG1|Orca_retrogene_ch7|Orca_parent'

pairs <- read.delim(text = mapping, sep = "|", stringsAsFactors = FALSE)
stopifnot(file.exists(input_file), min_species >= 2L)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

aln <- read.dna(input_file, format = "fasta",
                as.character = TRUE, as.matrix = TRUE)
aln[] <- toupper(aln)
rownames(aln) <- sub("\\s.*$", "", rownames(aln))
stopifnot(!anyDuplicated(rownames(aln)))
aln[aln == "."] <- "-"
ids <- unique(c(pairs$Parent, pairs$Retrogene))
missing <- setdiff(ids, rownames(aln))
if (length(missing)) stop("Missing sequences: ", paste(missing, collapse = ", "))
aln <- aln[ids, , drop = FALSE]
# Positions outside each sequence's first/last non-gap character are
# conservatively labeled No coverage, not interpreted as deletions.
# This does not establish whether the biological copy is truncated.
positions <- seq_len(ncol(aln))
covered <- function(x) {
  occupied <- which(x != "-")
  if (!length(occupied)) return(rep(FALSE, length(x)))
  positions >= min(occupied) & positions <= max(occupied)
}

d <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  parent <- aln[pairs$Parent[i], ]
  rg <- aln[pairs$Retrogene[i], ]
  coverage <- covered(parent) & covered(rg)
  known <- parent %in% c("A", "C", "G", "T", "-") &
    rg %in% c("A", "C", "G", "T", "-")
  both_gaps <- parent == "-" & rg == "-"
  comparable <- coverage & known & !both_gaps
  change <- ifelse(!coverage, "NC",
                   ifelse(!known, "?", ifelse(both_gaps, "--",
                                              ifelse(parent == rg, ".", paste0(parent, ">", rg)))))
  display <- ifelse(!coverage, "No coverage",
                    ifelse(!known, "Ambiguous", ifelse(both_gaps, "Both gaps",
                                                       ifelse(parent == rg, "Match", ifelse(rg == "-", "RG gap", rg)))))
  data.frame(
    Species = pairs$Species[i], Group = pairs$Group[i],
    Label = paste(pairs$Species[i], pairs$Copy[i], sep = " | "),
    Position = positions, Parent_base = parent, RG_base = rg,
    Comparable = comparable,
    Different = comparable & parent != rg,
    Change = change, Base = display
  )
}))
# A column remains usable if at least one parent/RG pair is comparable.
screened_positions <- unique(d$Position[d$Comparable])
d <- d[d$Position %in% screened_positions, ]
# reordering within the toothed-whale panel.

toothed_order <- c(
  "Sperm whale | RG1",
  "Pygmy sperm whale | RG1",
  
  # Remaining RG1 copies
  "Beluga | RG1",
  "Narwhal | RG1",
  "Bottlenose dolphin | RG1",
  "Orca | RG1",
  
  # RG2 copies
  "Beluga | RG2",
  "Narwhal | RG2",
  "Bottlenose dolphin | RG2",
  "Orca | RG2"
)

# Keep the existing baleen-whale order.
baleen_order <- unique(d$Label[d$Group == "Baleen"])
row_order <- c(baleen_order, toothed_order)

# Check that every row is included exactly once.
stopifnot(
  !anyDuplicated(row_order),
  setequal(row_order, unique(d$Label))
)

write.csv(d, file.path(output_dir, "all_screened_comparisons.csv"), row.names = FALSE)

# Count species, not copies, sharing each exact parent/RG combination.
events <- unique(d[d$Different, c("Position", "Change", "Species")])
if (!nrow(events)) stop("No differences found. screened comparisons saved.")
counts <- aggregate(Species ~ Position + Change, events,
                    function(x) length(unique(x)))
names(counts)[3] <- "Species_count"
write.csv(counts, file.path(output_dir, "difference_species_counts.csv"), row.names = FALSE)
sites <- sort(unique(counts$Position[counts$Species_count >= min_species]))
if (!length(sites)) stop("No differences meet the sharing threshold. tables saved.")

# Denominators are position-specific and include only comparable pairs.
# For species counts, a species qualifies if at least one RG is comparable;
# sharing means at least one copy has that exact difference.
selected_events <- counts[counts$Species_count >= min_species, ]
summary <- do.call(rbind, lapply(seq_len(nrow(selected_events)), function(i) {
  pos <- selected_events$Position[i]
  change <- selected_events$Change[i]
  do.call(rbind, lapply(c("Baleen", "Toothed"), function(g) {
    z <- d[d$Position == pos & d$Group == g, ]
    hit <- z$Comparable & z$Change == change
    nc <- sum(z$Comparable)
    ns <- length(unique(z$Species[z$Comparable]))
    hc <- sum(hit)
    hs <- length(unique(z$Species[hit]))
    data.frame(Position = pos, Change = change, Group = g,
               RG_copies_sharing = hc, RG_copies_comparable = nc,
               Species_sharing = hs, Species_comparable = ns,
               Copy_fraction = paste0(hc, "/", nc),
               Species_fraction = paste0(hs, "/", ns))
  }))
}))
write.csv(summary, file.path(output_dir, "shared_differences_by_group.csv"),
          row.names = FALSE)
writeLines(c(
  "NC: outside the overlapping spans of the parent and RG; treated as unknown.",
  "Terminal absence could be incomplete extraction or a true truncation; not resolved here.",
  "Internal gaps are retained; ambiguous bases and both-gap cells are not comparable.",
  "Each difference must occur in at least two distinct species globally to select a site.",
  "Denominators vary by position. Species counted once if any copy is comparable.",
  "A 0/0 fraction means no comparable data, not zero divergence.",
  "Gap columns are not separate biological indel events. Alignment is assumed correct."
), file.path(output_dir, "method_notes.txt"))

# All RGs are displayed at selected sites, including non-shared differences.
d <- d[d$Position %in% sites, ]
write.csv(d, file.path(output_dir, "plotted_comparisons.csv"), row.names = FALSE)
d$Position <- factor(d$Position, levels = sites)
d$Label <- factor(d$Label, levels = rev(row_order))
d$Group <- factor(d$Group, levels = c("Baleen", "Toothed"))

p <- ggplot(d, aes(Position, Label, fill = Base)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Change), size = 2.5) +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(
    Match = "#EEEEEE", "Both gaps" = "#FFFFFF", "RG gap" = "#999999",
    "No coverage" = "#DDD7EE", "Ambiguous" = "#FFF2B2",
    A = "#8DD3C7", C = "#80B1D3", G = "#FDB462", T = "#FB8072"
  ), name = "RG base / gap status") +
  labs(
    title = "Shared parent-retrogene differences, accounting for coverage",
    subtitle = paste(length(screened_positions), "alignment columns screened,",
                     length(sites), "shown. Depicts exact difference in at least",
                     min_species, "species"),
    x = "Original alignment column (selected positions only)", y = NULL,
    caption = paste("Text: parent > RG. '-' = gap; '.' = matching bases; '--' = both gaps.",
                    "\nNC = no overlapping coverage (excluded from counts); ? = ambiguous base. Internal gaps retained.")
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom")
print(p)
ggsave(file.path(output_dir, "shared_differences_colors_sep28.png"),
       p, width = max(12, length(sites) * 0.35 + 4), height = 8,
       dpi = 300, bg = "white", limitsize = FALSE)
message("Results saved in: ", normalizePath(output_dir, winslash = "/"))

#diff shades of blue
scale_fill_manual(
  values = c(
    "A"           = "#BDD7E7",  # Light blue
    "C"           = "#6BAED6",  # Medium blue
    "G"           = "#BCBDDC",  # Light purple
    "T"           = "#807DBA",  # Medium purple
    "Match"       = "#F0F0F0",  # Light grey
    "RG gap"      = "#969696",  # Darker grey
    "Both gaps"   = "#FFFFFF",  # White
    "No coverage" = "#E2D5ED",  # Pale lavender
    "Ambiguous"   = "#C7C7C7"   # Medium grey
  ),
  name = "RG base / gap status"
) +