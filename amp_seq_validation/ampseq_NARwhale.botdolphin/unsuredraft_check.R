B_denovo <- seqs[["NARW01_B_denovo"]]
pcna_like <- seqs[["REF_PCNA_like"]]

core_match <- matchPattern(
  pcna_like,
  B_denovo
)

core_match
# 116–901

B_core <- subseq(
  B_denovo,
  start = start(core_match)[1],
  end = end(core_match)[1]
)

B_core == pcna_like
# TRUE