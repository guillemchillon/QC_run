# RNA-seq Bulk Analysis using R: QC Run for STXBP1 haploinsufficient mice

## Project Overview

This project involves processing RNA-seq data on Respublica before performing statistical analysis in RStudio.

## Directory Structure

- `data/`: 
  - `counts/`: Raw count matrices for different samples.
  - `metadata/`: Experimental design files.
  - `multiqc/`: MultiQC quality control reports.
- `docs/`: Placeholder for additional documentation (e.g., notes, detailed methodology).
- `renv/`: R environment setup for reproducibility.
- `results/`: Post-trimming QC.
  - `figures/`: Plots genereated during the analysis.
  - `reports/`: Output reports.
  - `tables/`: Processed data tables.
  - `versions/`: Archives of result versions.
- `scripts/`: R scripts for analysis.
- `renv.lock`: Captures exact versions of R pacakges.
- `README.md`: Project documentation.

## Workflow Summary

1. **Data Import**: step1_importData.R. Load counts matrices and metadata.
2. **Data Wrangling**: step2_dataWrangling.R. Perform filtering and normalization.
3. **PCA Analysis**: step3_PCA.R. Conduct principal component analysis.
4. **Differential Gene Analysis**: step4_diffGenes.R. Identify differentially expressed genes using DESeq2. Generate DEG tables and visualizations (not finished).

## Usage Instructions

1. Clone repository: git clone https://github.com/guillemchillon/QC_run.git
2. Use the renv.lock file to restore the exact package versions: renv::restore()

## Contact

For questions, contact Guillem Chillon at chillonbog@chop.edu.

