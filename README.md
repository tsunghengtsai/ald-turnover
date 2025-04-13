# ald-turnover

This repository hosts a set of R scripts for reproducing the MS-based LFQ and dynamic analyses and figures presented in the manuscript: 

"Chronic alcohol consumption reprograms hepatic metabolism through organelle-specific acetylation in mice."

## Datasets

The following supplementary data files are required.

- `ALD_MS.xlsx`
- `ALD_NonMS.xlsx`

## Analyses

### MS-based differential and dynamic proteomics analyses

Run the source file `ALD25_MS.R` to reproduce the following analyses:

- Differential analysis of total protein abundance
- Differential analysis of acetylation level
- Differential analysis of native protein turnover
- Differential analysis of acetylated protein turnover
- Associations between protein abundance, acetylation level, and protein turnover
- Functional enrichment analysis

With these analyses, the following figures are produced: 2B, 3D-3F, 4D-4G, S1B, S2, S6, S7B, S9

### Others 

Run the source file `ALD25_other.R` to reproduce the other analyses and figures in the manuscript: 1B-1C, 2A, 3A, 3C, 4A-4C, 5B, 5D-5E, 6C, 7A-7F, 8, S7A, S8, S10, S11, S12
