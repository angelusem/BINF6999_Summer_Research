# Amplicon Sequencing Validation of PCNA Copy Number in Beluga and Blue Whale
The hypothesis is that beluga has 3 PCNA copies (1 parent gene \+ 2 retrogenes) and blue whale has 2 copies (1 parent gene \+ 1 retrogene). ddPCR on samples of beluga whale (B12) and blue whale confirmed the copy numbers and amplicon sequencing (ONT) was performed to characterise the individual copies at the sequence level. Analysis was run through the wf-amplicon workflow on Nextflow, through the AAC at the UoG. When performing variant calling, this first run used a consensus sequence from the provided samples as a reference for each amplicon to map against: “do any reads differ (and where) from the average of all reads” \- unsure about utility.

Objective: re-run amplicon workflow and variant calling using reference genome sequences: “how does individual B12 \+ our blue whale sample differ from reference genome”

PROGRESS LOGS

*21 May 2026:*

**Established key biological concepts:**

* **Parent gene for PCNA \=** intron-containing, on chromosome 14/15 based on related species. not mapped to chromosome in beluga assembly (ASM228892v3). Located on unplaced scaffold NW\_022098051.1  
* **Retrogene: c**opy derived from reverse transcription of processed mRNA back into the genome; has no introns (all exons fused continuously). identified by BLAST of exon sequence against genome  
* **Paralog:** gene related to another within the same organism due to duplication. so in beluga, for example, all three PCNA copies are paralogs of each other. Paralogs are distinct from alleles (same locus, two chromosomes) and orthologs (same gene across species)  
* Amplicon \= PCR product

**Genomic information available:**

**A. Beluga (Delphinapterus leucas)**

* **Assembly:** ASM228892v3, RefSeq accession GCF\_002288925.2  
* GCA\_002288925.3 was also downloaded initially but is the GenBank version with different scaffold identifiers, so not used. All work uses GCF (RefSeq) version to match NW\_ scaffold accessions returned by BLAST.  
* Genome has 5,904 unplaced scaffolds. PCNA copies not on named chromosomes

| Copy | Scaffold | Coordinates | Strand | % Identity to Parent |
| ----- | ----- | ----- | ----- | ----- |
| Parent gene | NW\_022098051.1 | 7,066,550–7,072,506 | Plus | 100% |
| Retrogene 1 | NW\_022098045.1 | 8,167,997–8,168,781 | **Minus** | 98.73% |
| Retrogene 2 | NW\_022098076.1 | 94,206,118–94,206,903 | **Minus** | 96.31% |

Retrogene 1 vs Retrogene 2: 95% identity, same length (786 bp). Both retrogenes confirmed intronless (single continuous BLAST hit each vs 6 separate exon hits for parent gene). Coordinates validated by NCBI BLAST of beluga PCNA CDS (XM\_022556493.2) against ASM228892v3.

**B. Blue Whale (Balaenoptera musculus)**

* **Assembly:** mBalMus1.pri.v3, RefSeq accession GCF\_009873245.2  
* Chromosome-level assembly. 23 chromosomes \+ 82 unplaced scaffolds \+ MT (106 sequences total) \*\* MT \= sequence corresponding to the circular mitochondrial genome. unplaced scaffold \= Contiguous sequences (contigs) that could not be mapped to the primary chromosomes due to highly repetitive regions, sequencing gaps, or structural variations  
* NCBI web BLAST returned NC\_041915028.1 for the retrogene. this accession does NOT exist in the downloaded GCF\_009873245.2 file. This is because web BLAST used a different assembly patch version. Local BLAST (makeblastdb \+ blastn) against the downloaded file confirmed the retrogene is on **NC\_045791.1 (chromosome 7\)**.

| Copy | Scaffold | Coordinates | Strand | % Identity to BW CDS |
| ----- | ----- | ----- | ----- | ----- |
| Parent gene | NC\_045799.1 (chr 15\) | 32,113,534–32,120,503 | Minus | 100% |
| Retrogene | NC\_045791.1 (chr 7\) | 111,914,528–111,915,313 | **Minus** | 97.84% |

**Blue whale has only 1 retrogene, not 2\.** Retrogene 2 either arose after divergence from blue whale, or was lost in blue whale. The BlueWhaleRetrogene2B sample (designed with beluga Retrogene 2 flanking primers) likely amplified non-specifically.its mapping rate in the pipeline will confirm this.

**Workflow:**

ddPCR established copy number (confirmed \~3 copies in beluga), which justified amplicon sequencing to characterise each copy individually: identifying paralog-specific variants (PSVs), confirming retrogene architecture (intron absence), and assessing functional integrity of each copy.

Two primer set strategies were used:

**A primer sets \-** designed within conserved exon sequence shared across all copies:

* PCNA-A Parent: spans introns and exons  
* PCNA-A Retrogenes: same primers for both retrogenes based on their 95% homology (?), sensitive to annealing temperature, which was optimized to preferentially amplify short intron-free retrogene products. This amplifies **both** retrogenes simultaneously in one pool

**B primer sets:** one primer in unique genomic flanking region (150 bp outside each retrogene insertion site) \+ one primer in exon body:

* Retrogene1B: specific to Retrogene 1 insertion site; \~1086 bp product since exon only \~786 bp  
* Retrogene2B: specific to Retrogene 2 insertion site; \~1086 bp product

The flanking regions are unique to each insertion site, making these assays copy-specific.

wf- amplicon pipeline: input raw reads → aligned to a short target reference → Medaka for polished consensus sequence → variant calling


\*\* Issue \#1 reference.fasta needed requires extracted sequences from the reference genome. Script used is saved inside folder as nextflow_ref.sh
\*\* Issue \#2 multiple amplicons per barcode are supported but only if they are different enough. It is possible that with 95% sequence homology, retrogene1 (R1) and 2 (R2) are too similar to be pooled together in ARetrogenes sample.
* if ARetrogene primer set amplifies both R1 and R2 simultaneously from the same rxn then the resulting pool of reads is a mix of 2 distinct templates that are 95% similar to each other ⇒ different positions are \~ 5% of the 786 bp so \~ 38 positions, at these positions Medaka will see roughly a 50/50 split of 2 alleles. Haploid variant caller interprets this 50/50 split as noise and wont call 2 alleles confidently versus at the 95% where R1 and R2 agree it will be confident that there is no variant  
* options:  
  * use 2 references \= \*\* did not work with wf-amplicon, does not support 2 references mapped to same sample  
  * amplicon “demultiplexing” by known paralog specific variants from correctly referenced Retrogene1A and 1B  
  * use a variant caller capable of resolving 2 distinct haplotypes from a mixed read pool

*28 May 2026:* troubleshooting pipeline, running into errors relating to sample.csv input that stop its progress

*2 June 2026:* re-ran pipeline successfully with sample.csv. barcodes were assigned uniquely to each. For ARetrogenes sample, retrogene1 exon sequence was used. need to look into downstream analysis.

**Interpretation of results:**

**Variant calling:** no results for retrogene with flanking regions and blue whale amplicons

**Next steps**

* investigating the utility of running a separate analysis for variant calling using minimap2 for alignment and then Clair3 for variant calling  
* exploring multi-mapping and high sequence identity between parent gene and retrogene- are there specific tools for my purposes  
* confirm that cutadapt recommended by prof. Taylor would not be suited for purposes of “correct” variant calling with the sample that was intended to “pull out” both retrogenes 1 and 2 in beluga whale (used the same primers for both of them)


