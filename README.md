# MASLD Sex-Stratified Multi-Omics Analysis

Reproducible R workflow for sex-stratified multi-omics analysis of metabolic dysfunction-associated steatotic liver disease (MASLD). The workflow integrates liver transcriptomics, lipidomics, metabolomics, and serum metabolomics to evaluate molecular modules, histological associations, cross-omics relationships, and mediation patterns.

> **Repository status:** Analysis code and documentation only. Individual-level phenotype and omics data are not included because their public release is restricted.

## Overview

This repository contains the analysis workflow used to investigate molecular patterns associated with liver histology in men and women. The primary R Markdown document covers:

- phenotypic characterization;
- preprocessing and quality control;
- weighted gene co-expression network analysis using WGCNA;
- module eigengene estimation;
- module-histology association testing;
- sex-stratified correlation analyses;
- chord diagrams of significant omics-histology associations;
- cross-omics correlation analysis;
- differential abundance and enrichment analyses;
- mediation analysis;
- selected predictive modeling and visualization workflows;
- layer-specific size-matched sensitivity analysis.

The repository is intended to support methodological transparency and reproducibility without distributing restricted participant-level data.

## Main analysis file

The complete workflow is provided in:

```text
EDC_MASLD_complete_analysis.Rmd
```


## Study design represented in the code

The workflow evaluates four omics layers:

| Abbreviation | Omics layer |
|---|---|
| SMC | Serum metabolomics |
| MC | Liver metabolomics |
| LC | Liver lipidomics |
| TC | Liver transcriptomics |

Module-histology associations are evaluated for three binary histological traits:

- steatosis;
- inflammation;
- fibrosis.

Associations used in the chord-diagram analysis are based on Spearman correlations with false discovery rate adjustment performed separately within each omics layer. The discovery threshold is FDR < 0.10.

## Repository structure

A recommended local structure is shown below. Restricted data directories are listed for reproducibility but must not be committed to a public repository.

```text
masld-sex-stratified-multiomics/
├── README.md
├── EDC_MASLD_complete_analysis.Rmd
├── .gitignore
└── Figures/                               
```

The folder name `Analyisis` is retained because it matches the directory referenced by the current workflow.

## Data availability and privacy

Individual-level clinical and omics data are not included in this repository. The analysis expects locally available input files, including phenotype, metabolomics, lipidomics, serum metabolomics, and transcriptomics data.

Do not commit any of the following unless public sharing has been explicitly authorized:

- participant identifiers;
- phenotype or clinical metadata;
- raw or normalized omics matrices;
- sample-level intermediate objects;
- spreadsheets containing individual-level observations;
- rendered reports that display sample IDs or participant-level tables.

Before publishing generated Markdown, HTML, CSV, Excel, figure, or cache files, review them for participant identifiers, row names, local file paths, and individual-level values.

## Software requirements

The analysis is written in R and uses packages from CRAN and Bioconductor. Major packages include:

- `tidyverse`
- `WGCNA`
- `psych`
- `circlize`
- `limma`
- `DESeq2`
- `clusterProfiler`
- `org.Hs.eg.db`
- `mediation`
- `vegan`
- `gtsummary`
- `ggplot2`
- `pheatmap`
- `MetaboAnalystR`

The complete package list is provided in the package-loading section of the R Markdown document. Package availability can vary by R and Bioconductor version. For long-term reproducibility, record the versions used for the final analysis with `sessionInfo()` or an `renv.lock` file.

## Running the analysis

### 1. Clone the repository

```bash
git clone YOUR_REPOSITORY_URL
cd masld-sex-stratified-multiomics
```

### 2. Prepare the local data structure

Place authorized input files in the folders expected by the R Markdown workflow. These files should remain local and should be protected by `.gitignore`.

### 3. Render the GitHub document

From R:

```r
rmarkdown::render(
  input = "EDC_MASLD_complete_analysis.Rmd",
  output_format = "github_document",
  params = list(
    project_dir = ".",
    results_dir = "results"
  )
)
```


### 4. Review generated files

Before committing generated reports or outputs, inspect them for:

- sample identifiers;
- participant-level values;
- local file paths;
- restricted filenames or metadata;
- unexpectedly large files.


## Reproducibility notes

- The size-matched analysis uses the fixed random seed `20260811`.
- Female participants are sampled without replacement.
- Matching is performed separately within each omics layer.
- The resampling analysis uses 10,000 iterations.
- Spearman correlation is used for module-histology associations.
- FDR adjustment is performed separately within each omics layer.
- Significant associations are defined using FDR < 0.10.
- The empirical lower-tail probability uses a plus-one correction.

## Responsible interpretation

The analysis is exploratory and uses a discovery-oriented FDR threshold. Results should be interpreted together with the study design, data availability, multiple-testing procedure, and lack of independent replication where applicable.

Sex-stratified findings should be described as associations within the analyzed cohort. The size-matched sensitivity result does not establish a biological sex difference in network architecture.

## Citation

If you use this workflow, cite the associated manuscript when its full citation becomes available. A `CITATION.cff` file can be added after the publication details and persistent identifier have been finalized.

Suggested temporary citation:

```text
Gupta S. MASLD Sex-Stratified Multi-Omics Analysis. Analysis code repository.
```

## Author

**Shashank Gupta**  
Örebro University  
[GitHub profile](https://github.com/shashank-KU)  
[Örebro University profile](https://www.oru.se/english/employee/shashank_gupta)

## License

No open-source license has been assigned yet. Until a license is added, reuse, redistribution, and modification are not granted beyond rights provided by applicable law or the hosting platform's terms.

## Contact and issues

For questions about the analysis workflow, open an issue in this repository without including participant-level information, restricted data, credentials, or confidential study details.
