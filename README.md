# HCL_manuscript
Scripts and data from the manuscript entitled "Multiparametric profiling of circulating extracellular vesicles supports their potential as liquid-biopsy biomarkers in classic Hairy Cell Leukemia".

This repository has two folders:
- Scripts, that contains the scripts for running all the analyses (sorted from 1 to 14);
- Counts, which has the count matrices (after QC, trimming and alignment with kallisto - as in scripts 1, 2 and 3 from Scripts folder).

"data_Filippo.RData" has gene-level counts.
The three "txi_....RData" files have transcript-level counts, and can be aggregated as follows in R:

```load("scripts/RData/txi_abundance.RData")
load("scripts/RData/txi_counts.RData")
load("scripts/RData/txi_len.RData")
txi$counts = txi_counts
txi$length = txi_len```
