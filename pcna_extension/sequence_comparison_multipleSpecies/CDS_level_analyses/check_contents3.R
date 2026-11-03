# Run once if needed:
# install.packages(c("ape", "ggplot2", "BiocManager"))
# BiocManager::install(c("ggtree", "Biostrings"))

library(ape)
library(ggtree)
library(ggplot2)
library(Biostrings)

folder <- "C:/Users/angel/Downloads/Research project summer"

tree <- read.tree(file.path(
  folder, "pcna_rooted_cow_ncbinarwhal.treefile"
))

aln <- readDNAStringSet(file.path(
  folder,
  "pcnaRG_plus_parents_ncbinarwhal_all_hippo_cow_aligned.fasta"
))

# Match FASTA names to IQ-TREE labels.
names(aln) <- sub("\\s.*$", "", names(aln))
names(aln) <- gsub(":", "_", names(aln), fixed = TRUE)

# Check that all sequences match and have equal aligned lengths.
stopifnot(
  !anyDuplicated(names(aln)),
  setequal(tree$tip.label, names(aln)),
  length(unique(width(aln))) == 1L
)

# Write a temporary plotting copy; the original FASTA is unchanged.

plot_fasta <- tempfile(fileext = ".fasta")

# Convert alignment to a matrix: sequences × alignment columns.
mat <- as.matrix(aln)

# Keep columns containing at least two different A/C/G/T bases.
# Gaps and ambiguous bases do not count as differences here.
variable <- apply(mat, 2, function(column) {
  bases <- column[column %in% c("A", "C", "G", "T")]
  length(unique(bases)) > 1L
})

if (!any(variable)) stop("No variable nucleotide sites found.")

aln <- DNAStringSet(
  setNames(
    apply(mat[, variable, drop = FALSE], 1, paste0, collapse = ""),
    rownames(mat)
  )
)
writeXStringSet(aln, plot_fasta)

p <- ggtree(tree) +
  geom_tiplab(size = 3)

tree_width <- diff(range(p$data$x, na.rm = TRUE))

fig <- msaplot(
  p,
  fasta = plot_fasta,
  offset = tree_width,
  width = 2
)

print(fig)

ggsave(
  file.path(folder, "pcna_tree_alignment.pdf"),
  plot = fig, width = 16, height = 9, units = "in"
)

unlink(plot_fasta)