# AI Validation Log

This document records significant uses of generative AI during the BINF*6999 project and how the output was handled. As the researcher, I acknowledge that I remain responsible for all code, analyses, biological interpretations, and written claims.

## Principles

- AI-generated code is reviewed for comprehension, sense and utility for own purposes before execution.
- Commands are tested before use on primary data.
- File paths, reference accessions, coordinate systems, sequence orientation, and sample labels are independently and manually verified.
- Biological claims are checked against literature.
- Numerical outputs are regenerated from source data and compared with independent checks where practical.
- AI output that cannot be validated is not used as evidence for a conclusion.

| Date | AI-assisted task | Problem or risk identified | How I checked or corrected it | Outcome |
|---|---|---|---|---|
| 2026-08-13 | Drafted repository organization, a root README, and a Conda environment file | The first environment setup failed because `environment.yml` had not been placed in the repository. Commands that followed were therefore run in the base environment. | I checked the terminal output, confirmed the missing file, and changed the procedure to create and verify the file before building and activating the environment. | Corrected; environment still needs to be tested successfully. |
| Approx. 2026-06-11 | Investigated whether RetroScan could help distinguish PCNA parent genes and retrocopies | AI could not confirm the intended RetroScan tool and suggested RetroSeq, which detects non-reference transposable-element or viral insertions and did not clearly match my analysis goal. | I read the relevant paper and compared the tool's stated purpose with my need to examine highly similar PCNA copies. I did not substitute RetroSeq solely because of the name similarity. | AI suggestion not adopted. |
| Date unknown, likely June-July 2026 | Generated candidate tools for inspecting ORFs, translations, variants, and sequence alignments | A broad tool list did not establish that every program was appropriate for my data or available annotations. Consequence annotation tools such as VEP, snpEff, or `bcftools csq` require suitable genomic coordinates and annotations. | I manually shortlisted tools based on their documented purpose, familiarity, prominence, and compatibility with the available sequence and annotation data. | Used as a starting list only. |
| 2026-06-29 | Planned retrieval of skin or proliferative-tissue RNA-seq data for bottlenose dolphin | AI initially suggested querying all 575 runs, which was unnecessarily broad and ignored the SRA study hierarchy. | After further review, I used the BioProject/SRA study -> BioSample -> run hierarchy and screened studies before individual runs. | Corrected to a study-first retrieval strategy. |
| Approx. 2026-07 to 2026-08-13 | Designed PCNA sequence-identity heatmap comparisons | AI suggested using the NCBI retrogene sequence to define the comparison region. This risked a circular and biased comparison because the query and comparison region would be defined by the same reference sequence. | I independently defined the North Atlantic right whale B PCNA-homologous core by orienting and aligning the primer-trimmed de novo consensus to the parental PCNA CDS before comparing sequences. | Corrected to a less circular comparison. |
