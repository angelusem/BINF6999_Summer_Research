setwd("C:/Users/angel/epi2melabs/instances/wf-amplicon_barcoded_addnl_r2_belugandblue_docmdefaults")

# install once if needed:
# BiocManager::install("Biostrings")

library(Biostrings)

files <- list.files(
  "output",
  pattern = "medaka.consensus.fasta$",
  recursive = TRUE,
  full.names = TRUE
)

seqs <- do.call(c, lapply(files, readDNAStringSet))
names(seqs) <- sub("^output/([^/]+)/.*", "\\1", gsub("\\\\", "/", files))

seqs
width(seqs)
alphabetFrequency(seqs, baseOnly = FALSE)
# to view the sequences as a table
consensus_df <- data.frame(
  sample = names(seqs),
  length = width(seqs),
  sequence = as.character(seqs)
)

View(consensus_df)
# to create a matrix that gives the number of differing positions between each consensus sequence pair
consensus_df$N_count <- letterFrequency(seqs, "N")
consensus_df

consensus_strings <- as.character(seqs)
outer(consensus_strings, consensus_strings, Vectorize(function(a, b) {
  sum(strsplit(a, "")[[1]] != strsplit(b, "")[[1]])
}))
######
# basc qc
width(seqs)
letterFrequency(seqs, letters = c("A", "C", "G", "T", "N", "-"))
alphabetFrequency(seqs)
# pairwise comparison 
consensus_strings <- as.character(seqs)

diff_matrix <- outer(consensus_strings, consensus_strings, Vectorize(function(a, b) {
  sum(strsplit(a, "")[[1]] != strsplit(b, "")[[1]])
}))

diff_matrix
# multiple sequence alignment 
# install once if needed:

BiocManager::install("msa")

library(msa)

alignment <- msa(seqs)
print(alignment)

#
alignment_matrix <- as.matrix(alignment)

variable_sites <- apply(alignment_matrix, 2, function(x) length(unique(x)) > 1)

which(variable_sites)
# to see the actual bases at these positions
alignment_matrix[, variable_sites]
