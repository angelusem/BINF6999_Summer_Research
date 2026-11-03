# Standalone shared AMINO-ACID parent/RG differences, including alignment gaps.
# Run once if needed: install.packages("ggplot2")
# Open this complete file in RStudio and click Source.
library(ggplot2)

input_file <- "C:/Users/angel/Downloads/Research project summer/exploratory_analyses_15_seqs/protein_seq_analysis/MSA_baleen_protein.fa"
output_dir <- file.path(dirname(input_file), "baleen_protein_differences")
min_species <- 1L
# Set to 1 to include differences unique to one species.
# Clean only the final _1 suffix; leave RG1/RG2 intact.

# Six baleen species; each RG compared with its own parent.
# Trailing _1 is removed; plot copy labels are RG.
mapping <- 'Species|Group|Copy|Retrogene|Parent
Bowhead whale|Baleen|RG|Bowhead_RG1|Bowhead_parent
North Atlantic right whale|Baleen|RG|NArightWhale_RG1|NARightWhale_parent
Blue whale|Baleen|RG|BlueWhale_RG1|BlueWhale_parent
Humpback whale|Baleen|RG|HumpbackWhale_RG1|HumpbackWhale_parent
Gray whale|Baleen|RG|GrayWhale_RG1|GrayWhale_parent
Minke whale|Baleen|RG|MinkeWhale_RG1|MinkeWhale_parent'

pairs <- read.delim(text = mapping, sep = "|", stringsAsFactors = FALSE)
stopifnot(file.exists(input_file), min_species >= 1L)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read protein FASTA directly: preserve amino acids, X, *, and gaps.
lines <- trimws(readLines(input_file, warn = FALSE))
lines <- lines[nzchar(lines)]
starts <- which(startsWith(lines, ">"))
if (!length(starts) || starts[1] != 1L) stop("Invalid FASTA.")
ends <- c(starts[-1] - 1L, length(lines))
headers <- substring(lines[starts], 2L)
original_ids <- sub("\\s.*$", "", headers)
clean_ids <- sub("_1$", "", original_ids)
if (anyDuplicated(clean_ids)) stop("Removing _1 produces duplicate identifiers.")
sequences <- vapply(seq_along(starts), function(i) {
  if (ends[i] <= starts[i]) stop("Empty FASTA sequence.")
  toupper(gsub("\\s", "", paste0(lines[(starts[i]+1L):ends[i]], collapse = "")))
}, character(1))
if (length(unique(nchar(sequences))) != 1L) stop("Sequences are not equally aligned.")
aln <- do.call(rbind, strsplit(sequences, "", fixed = TRUE))
rownames(aln) <- clean_ids
aln[aln == "."] <- "-"
# Save a cleaned copy; the original FASTA is unchanged.
clean_headers <- sub("_1(?=\\s|$)", "", headers, perl = TRUE)
clean_fasta <- unlist(lapply(seq_along(sequences), function(i) {
  c(paste0(">", clean_headers[i]), sequences[i])
}), use.names = FALSE)
writeLines(clean_fasta, file.path(output_dir, "protein_alignment_clean_headers.fasta"))
write.csv(data.frame(Original_ID = original_ids, Clean_ID = clean_ids),
          file.path(output_dir, "header_lookup.csv"), row.names = FALSE)
# Standard amino acids and explicit stop symbols are comparable.
# X and other ambiguous/unknown residues are excluded pair by pair.
valid <- c(strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]], "*", "-")
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
  known <- parent %in% valid & rg %in% valid
  both_gaps <- parent == "-" & rg == "-"
  comparable <- coverage & known & !both_gaps
  change <- ifelse(!coverage, "NC",
                   ifelse(!known, "?", ifelse(both_gaps, "--",
                                              ifelse(parent == rg, ".", paste0(parent, ">", rg)))))
  display <- ifelse(!coverage, "No coverage",
                    ifelse(!known, "Unknown AA", ifelse(both_gaps, "Both gaps",
                                                        ifelse(parent == rg, "Match", ifelse(rg == "-", "RG gap", rg)))))
  data.frame(
    Species = pairs$Species[i], Group = pairs$Group[i],
    Label = paste(pairs$Species[i], pairs$Copy[i], sep = " | "),
    Position = positions, Parent_AA = parent, RG_AA = rg,
    Comparable = comparable,
    Different = comparable & parent != rg,
    Change = change, Base = display
  )
}))
# A column remains usable if at least one parent/RG pair is comparable.
screened_positions <- unique(d$Position[d$Comparable])
d <- d[d$Position %in% screened_positions, ]
row_order <- unique(d$Label)
write.csv(d, file.path(output_dir, "all_screened_comparisons.csv"), row.names = FALSE)

# Count species, not copies, sharing each exact parent/RG combination.
events <- unique(d[d$Different, c("Position", "Change", "Species")])
if (!nrow(events)) stop("No differences found; screened comparisons saved.")
counts <- aggregate(Species ~ Position + Change, events,
                    function(x) length(unique(x)))
names(counts)[3] <- "Species_count"
write.csv(counts, file.path(output_dir, "difference_species_counts.csv"), row.names = FALSE)
sites <- sort(unique(counts$Position[counts$Species_count >= min_species]))
if (!length(sites)) stop("No differences meet the sharing threshold; tables saved.")

# Denominators are position-specific and include only comparable pairs.
# For species counts, a species qualifies if at least one RG is comparable;
# sharing means at least one copy has that exact difference.
selected_events <- counts[counts$Species_count >= min_species, ]
summary <- do.call(rbind, lapply(seq_len(nrow(selected_events)), function(i) {
  pos <- selected_events$Position[i]
  change <- selected_events$Change[i]
  do.call(rbind, lapply(c("Baleen"), function(g) {
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
  "Internal gaps are retained; unknown amino acids and both-gap cells are not comparable.",
  "Each difference must occur in at least min_species distinct species globally to select a site.",
  "Denominators vary by position. Species counted once if any copy is comparable.",
  "A 0/0 fraction means no comparable data, not zero divergence.",
  "Gap columns are not separate biological indel events. Alignment is assumed correct.",
  "This input already contains explicit * symbols. X is treated as unknown.",
  "Protein coordinates do not directly identify nucleotide coordinates or frameshifts.",
  "A>B denotes present-day parent/RG residues, not evolutionary mutation direction."
), file.path(output_dir, "method_notes.txt"))

# All RGs are displayed at selected sites, including non-shared differences.
d <- d[d$Position %in% sites, ]
write.csv(d, file.path(output_dir, "plotted_comparisons.csv"), row.names = FALSE)
d$Position <- factor(d$Position, levels = sites)
d$Label <- factor(d$Label, levels = rev(row_order))
d$Group <- factor(d$Group, levels = c("Baleen"))

# Plot categories distinguish amino-acid substitutions, gaps and unknowns.
# Exact residues are retained in labels and the exported CSV tables.
d$Status <- ifelse(d$Base %in% c("No coverage", "Unknown AA", "Both gaps", "Match", "RG gap"),
                   d$Base, ifelse(d$Parent_AA == "-", "Parent gap",
                                  ifelse(d$Parent_AA == "*" | d$RG_AA == "*",
                                         "Stop-related difference", "Amino-acid difference")))
# Color substitutions by the amino acid present in the RG.
d$Fill <- ifelse(
  d$Status == "Amino-acid difference",
  paste0("RG: ", d$RG_AA),
  d$Status
)

aa_labels <- sort(unique(
  d$Fill[d$Status == "Amino-acid difference"]
))

aa_colors <- setNames(
  grDevices::hcl.colors(length(aa_labels), palette = "Pastel 1"),
  aa_labels
)
# strong_colors <- c(
#   "#0072B2", "#E69F00", "#006400", "#FF69B4",
#   "#6A3D9A", "#8B4513", "#00BFC4", "#B2182B",
#   "#999900", "#000080", "#A6CEE3", "#FDBF6F",
#   "#B15928", "#CAB2D6", "#FFFF99", "#666666",
#   "#00A651", "#F0027F", "#80B1D3", "#BC80BD"
# )
# 
# aa_colors <- setNames(
#   strong_colors[seq_along(aa_labels)],
#   aa_labels
# )

fill_colors <- c(
  aa_colors,
  "Match" = "#F2F2F2",
  "RG gap" = "#969696",
  "Parent gap" = "#BCBDDC",
  "Both gaps" = "#FFFFFF",
  "No coverage" = "#DDD7EE",
  "Unknown AA" = "#FFF2B2",
  "Stop-related difference" = "#AD526B"
)
p <- ggplot(d, aes(Position, Label, fill = Fill)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Change), size = 2.5) +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(
    values = fill_colors,
    labels = function(x) sub("^RG: ", "", x),
    name = "Amino acid / status"
  ) +
  labs(
    title = "Shared baleen PCNA parent-retrogene amino-acid differences",
    subtitle = paste(length(screened_positions), "alignment columns screened.",
                     length(sites), "shown, with exact difference in at least",
                     min_species, "species"),
    x = "Protein alignment column (selected positions only)", y = NULL,
    caption = paste("Text: parent AA > RG AA. '-' = gap; '.' = matching residues; '--' = both gaps.")
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom")
print(p)
ggsave(file.path(output_dir, "baleen_protein_differences_all_pastel.png"),
       p, width = max(12, length(sites) * 0.35 + 4), height = 8,
       dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(output_dir, "baleen_protein_differences_all.pdf"),
       p, width = max(10, length(sites) * 0.45 + 4), height = 5,
       bg = "white", limitsize = FALSE)
message("Results saved in: ", normalizePath(output_dir, winslash = "/"))
