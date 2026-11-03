# Install once if needed:
# install.packages(c("ape", "ggplot2", "BiocManager"))
# 
BiocManager::install("ggtree")

library(ape)
library(ggtree)
library(ggplot2)

# Replace with your file paths (use forward slashes on Windows).
tree_file <- "C:/Users/angel/Downloads/Research project summer/pcna_rooted_cow_ncbinarwhal.treefile"
alignment_file <- "C:/Users/angel/Downloads/Research project summer/pcnaRG_plus_parents_ncbinarwhal_all_hippo_cow_aligned.fasta"

tree <- read.tree(tree_file)

p <- ggtree(tree) +
  geom_tiplab(size = 3)

# Set spacing relative to the tree's plotted width.
tree_width <- diff(range(p$data$x, na.rm = TRUE))

fig <- msaplot(
  p,
  fasta = alignment_file,
  offset = 0.3 * tree_width,
  width = 2
)

print(fig)

ggsave("tree_alignment.pdf", plot = fig,
       width = 14, height = 8, units = "in")