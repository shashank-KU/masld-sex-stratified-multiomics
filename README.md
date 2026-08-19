<p align="center">
  <img src="assets/banner.svg" alt="MASLD Sex-Stratified Multi-Omics Analysis" width="100%" />
</p>

<p align="center">
  <a href="https://www.r-project.org/"><img alt="R" src="https://img.shields.io/badge/R-4.4.1-276DC3?logo=r&logoColor=white" /></a>
  <a href="LICENSE"><img alt="License: MIT" src="https://img.shields.io/badge/License-MIT-yellow.svg" /></a>
  <a href="https://github.com/shashank-KU/masld-sex-stratified-multiomics"><img alt="Repository" src="https://img.shields.io/badge/GitHub-analysis%20workflow-181717?logo=github&logoColor=white" /></a>
  <a href="#-data-availability-and-privacy"><img alt="Data status" src="https://img.shields.io/badge/data-restricted-critical" /></a>
  <a href="MASLD_complete_analysis.Rmd"><img alt="Reproducibility" src="https://img.shields.io/badge/workflow-R%20Markdown-75AADB" /></a>
  <img alt="Maintenance" src="https://img.shields.io/badge/Maintained%3F-yes-brightgreen.svg" />
</p>

<p align="center">
  <img alt="Last commit" src="https://img.shields.io/github/last-commit/shashank-KU/masld-sex-stratified-multiomics" />
  <img alt="Repo size" src="https://img.shields.io/github/repo-size/shashank-KU/masld-sex-stratified-multiomics" />
  <img alt="Open issues" src="https://img.shields.io/github/issues/shashank-KU/masld-sex-stratified-multiomics" />
  <img alt="Stars" src="https://img.shields.io/github/stars/shashank-KU/masld-sex-stratified-multiomics?style=social" />
</p>

<p align="center">
A reproducible R workflow for <b>sex-stratified multi-omics analysis</b> of metabolic dysfunction-associated steatotic liver disease (MASLD). The workflow integrates liver transcriptomics, liver lipidomics, liver metabolomics, and serum metabolomics to evaluate molecular modules, histological associations, cross-omics relationships, mediation patterns, and sample-size sensitivity.
</p>

> **Repository scope:** This public repository contains analysis code, documentation, and approved figures only. Individual-level phenotype and omics data are **not** included because public release is restricted.

---

## Contents

- [🧬 Scientific scope](#-scientific-scope)
- [🔬 Analysis overview](#-analysis-overview)
- [📁 Repository structure](#-repository-structure)
- [🔧 Software environment](#-software-environment)
- [🚀 Running the workflow](#-running-the-workflow)
- [📊 Expected outputs](#-expected-outputs)
- [🔒 Data availability and privacy](#-data-availability-and-privacy)
- [✅ Reproducibility and quality control](#-reproducibility-and-quality-control)
- [⚠ Responsible interpretation](#-responsible-interpretation)
- [📄 Citation](#-citation)
- [📜 License](#-license)
- [👤 Author](#-author)
- [📧 Contact](#-contact)
- [🙏 Acknowledgment](#-acknowledgment)

---

## 🧬 Scientific scope

The workflow evaluates molecular patterns associated with three liver histological features: **steatosis**, **inflammation**, and **fibrosis**.

<p align="center">
  <img src="assets/omics-layers.svg" alt="Four integrated omics layers: TC, LC, MC, SMC" width="100%" />
</p>

| Abbreviation | Omics layer | Biological source |
|:---:|---|---|
| **SMC** | Serum metabolomics | Serum |
| **MC** | Metabolomics | Liver biopsy |
| **LC** | Lipidomics | Liver biopsy |
| **TC** | Transcriptomics | Liver biopsy |

The repository supports the following analytical objectives:

1. characterize clinical and histological features;
2. preprocess and quality-control multi-omics datasets;
3. identify omics modules using weighted correlation network analysis;
4. calculate module eigengenes and hub features;
5. evaluate module-histology associations;
6. perform sex-stratified analyses;
7. visualize significant associations using heatmaps and chord diagrams;
8. evaluate cross-omics correlations;
9. conduct differential and pathway-enrichment analyses;
10. assess mediation patterns across molecular layers;
11. evaluate sensitivity to unequal sample size and layer-specific data availability.

---

## 🔬 Analysis overview

```mermaid
flowchart LR
    A[Restricted local input data] --> B[Preprocessing and quality control]
    B --> C1[Liver transcriptomics]
    B --> C2[Liver lipidomics]
    B --> C3[Liver metabolomics]
    B --> C4[Serum metabolomics]
    C1 --> D[WGCNA modules and eigengenes]
    C2 --> D
    C3 --> D
    C4 --> D
    D --> E[Module-histology associations]
    D --> F[Cross-omics integration]
    D --> G[Mediation analysis]
    E --> H[Sex-stratified interpretation]
    E --> I[Layer-specific size-matched sensitivity analysis]
    F --> J[Figures and summary outputs]
    G --> J
    H --> J
    I --> J
```

### Core analytical methods

| Category | Method / tool | Description |
|---|---|---|
| Network analysis | **WGCNA** (biweight midcorrelation, signed) | Weighted co-expression network construction |
| Module associations | Spearman correlation | Module eigengenes vs. histological traits |
| Multiple testing | FDR adjustment | Threshold defined in the analysis script and manuscript |
| Transcriptomics | STAR → RSEM → DESeq2 | Alignment, quantification, normalization, differential analysis |
| Pathway analysis | clusterProfiler / MetaboAnalyst | Transcriptomic enrichment and metabolic pathway analysis |
| Mediation | `mediation` R package | Regression-based mediation analysis |
| Visualization | ggplot2, circlize, pheatmap, ggalluvial, cowplot, patchwork | Statistical and integrative visualization |

### Advanced governance-aware workflow

The full workflow separates restricted participant-level inputs from reproducible analytical steps and approved public outputs. Dashed arrows below indicate validation or sensitivity checks rather than the main analytical path.

<details>
<summary><b>📈 Show the full workflow diagram</b> — governance → preprocessing → networks → sex-stratified associations → integration → sensitivity → release</summary>

```mermaid
flowchart TB

    %% =========================
    %% INPUT AND GOVERNANCE LAYER
    %% =========================
    subgraph G0["Data governance and restricted inputs"]
        direction LR
        P1["Clinical and histological metadata<br/>age, sex, BMI, steatosis,<br/>inflammation, fibrosis"]
        P2["Liver transcriptomics<br/>raw RNA-seq reads"]
        P3["Liver lipidomics<br/>UHPLC-QTOFMS data"]
        P4["Liver metabolomics<br/>polar and semipolar features"]
        P5["Serum metabolomics<br/>systemic metabolic features"]
        PRIV["Restricted local environment<br/>participant-level data are not public"]

        P1 --> PRIV
        P2 --> PRIV
        P3 --> PRIV
        P4 --> PRIV
        P5 --> PRIV
    end

    %% =========================
    %% PREPROCESSING LAYER
    %% =========================
    subgraph G1["Preprocessing and quality control"]
        direction LR

        TPRE["Transcriptomics preprocessing<br/>STAR alignment<br/>RSEM quantification<br/>low-expression filtering<br/>DESeq2 transformation"]

        LPRE["Lipidomics preprocessing<br/>MSConvert to mzML<br/>MZmine feature processing<br/>internal-standard normalization<br/>log2 transformation and scaling"]

        MPRE["Liver metabolomics preprocessing<br/>MSConvert to mzML<br/>MZmine feature processing<br/>half-minimum imputation<br/>normalization and scaling"]

        SPRE["Serum metabolomics preprocessing<br/>feature filtering<br/>half-minimum imputation<br/>normalization and scaling"]

        QC["Quality-control gate<br/>sample identifiers aligned<br/>duplicate identifiers excluded<br/>missingness and feature variation checked"]
    end

    PRIV --> TPRE
    PRIV --> LPRE
    PRIV --> MPRE
    PRIV --> SPRE

    TPRE --> QC
    LPRE --> QC
    MPRE --> QC
    SPRE --> QC

    %% =========================
    %% OMICS NETWORK LAYER
    %% =========================
    subgraph G2["Layer-specific network analysis"]
        direction LR

        TSFT["TC soft-threshold selection<br/>signed network and bicor"]
        LSFT["LC soft-threshold selection<br/>signed network and bicor"]
        MSFT["MC soft-threshold selection<br/>signed network and bicor"]
        SSFT["SMC soft-threshold selection<br/>signed network and bicor"]

        TW["TC WGCNA<br/>modules, eigengenes,<br/>hub genes"]
        LW["LC WGCNA<br/>modules, eigengenes,<br/>hub lipids"]
        MW["MC WGCNA<br/>modules, eigengenes,<br/>hub metabolites"]
        SW["SMC WGCNA<br/>modules, eigengenes,<br/>hub metabolites"]

        TSFT --> TW
        LSFT --> LW
        MSFT --> MW
        SSFT --> SW
    end

    QC --> TSFT
    QC --> LSFT
    QC --> MSFT
    QC --> SSFT

    %% =========================
    %% ASSOCIATION LAYER
    %% =========================
    subgraph G3["Sex-stratified module-histology analysis"]
        direction TB

        SPLIT["Prespecified sex stratification<br/>women and men analyzed separately"]

        TRAITS["Histological traits<br/>steatosis, inflammation, fibrosis"]

        COR["Module-trait association testing<br/>Spearman correlation"]

        FDR["Within-layer multiple-testing control<br/>FDR threshold defined in the analysis"]

        SIG["Significant module-histology associations"]

        SPLIT --> COR
        TRAITS --> COR
        COR --> FDR
        FDR --> SIG
    end

    TW --> SPLIT
    LW --> SPLIT
    MW --> SPLIT
    SW --> SPLIT

    P1 --> TRAITS

    %% =========================
    %% INTEGRATION LAYER
    %% =========================
    subgraph G4["Multi-omics integration and interpretation"]
        direction LR

        CROSS["Cross-omics correlations<br/>TC-LC, TC-MC, LC-MC,<br/>and serum-linked comparisons"]

        DIFF["Differential analysis<br/>histology-associated molecular features"]

        TENRICH["Transcriptomic enrichment<br/>KEGG analysis with clusterProfiler"]

        MENRICH["Metabolic pathway analysis<br/>MetaboAnalyst"]

        MED["Mediation analysis<br/>transcriptomic module to mediator module<br/>to binary histological outcome"]

        VIS["Integrative visualization<br/>heatmaps, chord diagrams,<br/>alluvial plots, multipanel figures"]
    end

    TW --> CROSS
    LW --> CROSS
    MW --> CROSS
    SW --> CROSS

    SIG --> CROSS
    SIG --> DIFF
    TW --> TENRICH
    LW --> MENRICH
    MW --> MENRICH
    SW --> MENRICH
    CROSS --> MED

    CROSS --> VIS
    DIFF --> VIS
    TENRICH --> VIS
    MENRICH --> VIS
    MED --> VIS

    %% =========================
    %% SENSITIVITY AND VALIDATION
    %% =========================
    subgraph G5["Sensitivity analysis and internal validation"]
        direction TB

        AVAIL["Layer-specific sample availability<br/>SMC: 57 men and 152 women<br/>MC: 52 men and 137 women<br/>LC: 57 men and 154 women<br/>TC: 56 men and 154 women"]

        MVALID["Male-count validation<br/>original and recalculated counts agree<br/>SMC 0, MC 0, LC 16, TC 7<br/>total 23"]

        RESAMPLE["10,000 female resampling iterations<br/>sampling without replacement<br/>matched to male availability within each layer"]

        RECALC["Recalculate Spearman associations<br/>and within-layer FDR in every iteration"]

        DIST["Female matched distribution<br/>median 28, IQR 22 to 35<br/>2.5th to 97.5th percentile: 12 to 49"]

        EMP["Empirical comparison<br/>69.12% of iterations exceeded 23<br/>one-sided lower-tail probability 0.3089"]

        CAUTION["Interpretation gate<br/>apparent sex difference remains unresolved<br/>after matching sample size and availability"]

        AVAIL --> RESAMPLE
        MVALID --> RESAMPLE
        RESAMPLE --> RECALC
        RECALC --> DIST
        DIST --> EMP
        EMP --> CAUTION
    end

    SIG -. "audit observed counts" .-> MVALID
    QC -. "verify layer availability" .-> AVAIL
    FDR -. "repeat identical test family" .-> RECALC
    CAUTION -. "qualifies interpretation" .-> VIS

    %% =========================
    %% OUTPUT AND RELEASE LAYER
    %% =========================
    subgraph G6["Outputs, disclosure review, and release"]
        direction LR

        OUT1["Analysis outputs<br/>summary tables and figures"]
        OUT2["Reproducibility records<br/>R Markdown, session information,<br/>random seeds, software versions"]
        REVIEW["Disclosure review<br/>remove identifiers, local paths,<br/>restricted tables, caches, and logs"]
        PUBLIC["Public repository<br/>code, documentation,<br/>approved figures only"]
        ARCHIVE["Versioned archival release<br/>GitHub release and Zenodo DOI"]

        OUT1 --> REVIEW
        OUT2 --> REVIEW
        REVIEW --> PUBLIC
        PUBLIC --> ARCHIVE
    end

    VIS --> OUT1
    CAUTION --> OUT1
    QC --> OUT2
    G2 --> OUT2
    G5 --> OUT2

    %% =========================
    %% STYLING
    %% =========================
    classDef restricted fill:#FDECEC,stroke:#B42318,color:#7A271A,stroke-width:1.5px;
    classDef preprocessing fill:#FFF4E5,stroke:#B54708,color:#7A2E0E,stroke-width:1.2px;
    classDef network fill:#E8F3FF,stroke:#175CD3,color:#1849A9,stroke-width:1.2px;
    classDef analysis fill:#EEF4FF,stroke:#3538CD,color:#2D31A6,stroke-width:1.2px;
    classDef integration fill:#F4EBFF,stroke:#7F56D9,color:#53389E,stroke-width:1.2px;
    classDef validation fill:#ECFDF3,stroke:#039855,color:#027A48,stroke-width:1.2px;
    classDef caution fill:#FFF1F3,stroke:#E31B54,color:#C01048,stroke-width:1.5px;
    classDef output fill:#F2F4F7,stroke:#475467,color:#344054,stroke-width:1.2px;

    class P1,P2,P3,P4,P5,PRIV restricted;
    class TPRE,LPRE,MPRE,SPRE,QC preprocessing;
    class TSFT,LSFT,MSFT,SSFT,TW,LW,MW,SW network;
    class SPLIT,TRAITS,COR,FDR,SIG analysis;
    class CROSS,DIFF,TENRICH,MENRICH,MED,VIS integration;
    class AVAIL,MVALID,RESAMPLE,RECALC,DIST,EMP validation;
    class CAUTION caution;
    class OUT1,OUT2,REVIEW,PUBLIC,ARCHIVE output;
```

</details>

**Workflow interpretation**

1. **Restricted participant-level data remain local.** The public repository begins at the reproducible-code and approved-output layer, not at the raw-data layer.
2. **Each omics layer is preprocessed and networked separately.** Module eigengenes provide a common representation for downstream integration.
3. **Sex-stratified associations are evaluated consistently.** Steatosis, inflammation, and fibrosis are tested against module eigengenes using Spearman correlation with within-layer FDR control.
4. **Cross-omics and mediation analyses connect molecular layers.** Enrichment and visualization provide biological interpretation of module-level findings.
5. **The size-matched analysis is a sensitivity analysis of detected association counts.** It does not rebuild the WGCNA networks within every resample and therefore does not directly test intrinsic network-density differences.
6. **A disclosure-review gate separates internal analysis from public release.** Only code, documentation, and approved non-sensitive outputs are published.

### Main analysis file

The complete workflow is provided in [`MASLD_complete_analysis.Rmd`](MASLD_complete_analysis.Rmd) — the authoritative analysis source in this repository.

---

## 📁 Repository structure

```text
masld-sex-stratified-multiomics/
├── README.md
├── LICENSE
├── MASLD_complete_analysis.Rmd
├── Figures/                  
└── assets/                   # README graphics (banner, diagrams)
    ├── banner.svg
    └── omics-layers.svg
```

Restricted data directories are intentionally excluded from the public repository. The workflow expects authorized input files to be available locally in the directory structure referenced by the R Markdown document.

---

## 🔧 Software environment

The analysis was conducted in **R 4.4.1**.

<p align="left">
  <img alt="WGCNA" src="https://img.shields.io/badge/WGCNA-1.73-1849A9" />
  <img alt="DESeq2" src="https://img.shields.io/badge/DESeq2-1.44.0-B54708" />
  <img alt="clusterProfiler" src="https://img.shields.io/badge/clusterProfiler-4.12.6-175CD3" />
  <img alt="mediation" src="https://img.shields.io/badge/mediation-4.5.0-7F56D9" />
  <img alt="circlize" src="https://img.shields.io/badge/circlize-0.4.16-039855" />
  <img alt="ggplot2" src="https://img.shields.io/badge/ggplot2-3.5.1-3538CD" />
  <img alt="MetaboAnalystR" src="https://img.shields.io/badge/MetaboAnalystR-4.0.0-C11574" />
  <img alt="POMA" src="https://img.shields.io/badge/POMA-1.14.0-475467" />
</p>

<details>
<summary><b>Show full package table</b></summary>

<br/>

| Package or software | Version used | Primary purpose |
|---|---:|---|
| WGCNA | 1.73 | Network construction and module eigengenes |
| psych | 2.4.6.26 | Correlation analysis |
| limma | 3.60.4 | Linear-model analysis |
| DESeq2 | 1.44.0 | RNA-seq normalization and differential analysis |
| clusterProfiler | 4.12.6 | Functional enrichment |
| mediation | 4.5.0 | Mediation analysis |
| circlize | 0.4.16 | Chord diagrams |
| ggplot2 | 3.5.1 | Statistical visualization |
| POMA | 1.14.0 | Metabolomics data handling and preprocessing |
| ggalluvial | 0.12.5 | Alluvial visualizations |
| patchwork | 1.2.0 | Multipanel figure assembly |
| cowplot | 1.1.3 | Plot alignment and figure assembly |
| pheatmap | 1.0.12 | Heatmaps |
| MetaboAnalystR | 4.0.0 | Metabolomics analysis support |

Additional software used upstream includes ProteoWizard/MSConvert, MZmine, STAR, RSEM, and MetaboAnalyst. Exact versions and citations are reported in the manuscript Key Resources Table where available.

</details>

For exact dependency and platform information, add a `sessionInfo.txt` file generated from the final analysis environment:

```r
writeLines(
  capture.output(sessionInfo()),
  "sessionInfo.txt"
)
```

For stronger environment reproducibility, create an `renv.lock` file from the final working R environment:

```r
install.packages("renv")
renv::init()
renv::snapshot()
```

---

## 🚀 Running the workflow

### 1. Clone the repository

```bash
git clone https://github.com/shashank-KU/masld-sex-stratified-multiomics.git
cd masld-sex-stratified-multiomics
```

### 2. Prepare authorized local inputs

Place authorized data files in the local folders expected by `MASLD_complete_analysis.Rmd`. Do not place restricted participant-level files in the Git working tree unless the files are explicitly excluded and repository access is appropriately controlled.

### 3. Review paths and parameters

Open `MASLD_complete_analysis.Rmd` and verify:

- input directory paths;
- output directory paths;
- software-specific file locations;
- analysis parameters;
- availability of all required R and Bioconductor packages.

### 4. Render the analysis

From R:

```r
rmarkdown::render(
  input = "MASLD_complete_analysis.Rmd"
)
```

Or from a shell with R available:

```bash
Rscript -e 'rmarkdown::render("MASLD_complete_analysis.Rmd")'
```

### 5. Review generated content before committing

Inspect all rendered reports, tables, figures, caches, and logs for:

- participant or sample identifiers;
- individual-level clinical or omics values;
- restricted filenames;
- local user paths;
- temporary files;
- unexpectedly large outputs.

---

## 📊 Expected outputs

Depending on the available inputs and enabled analysis chunks, the workflow can generate:

- cohort summary tables;
- soft-threshold diagnostic plots;
- WGCNA module and eigengene summaries;
- module-histology correlation heatmaps;
- sex-stratified chord diagrams;
- cross-omics correlation outputs;
- differential-analysis tables;
- pathway-enrichment plots;
- mediation summaries and alluvial plots;
- layer-specific size-matched sensitivity distributions;
- manuscript and supplementary figures.

Approved figures can be stored in the `Figures/` directory. Generated participant-level or intermediate outputs must remain outside the public repository.

---

## 🔒 Data availability and privacy

Individual-level clinical and omics datasets are not distributed through this repository. Publicly available code does not confer access to restricted study data.

<details>
<summary><b>Show what must never be committed, and the recommended <code>.gitignore</code></b></summary>

<br/>

Do not commit:

- participant identifiers or sample keys;
- clinical metadata;
- raw or normalized omics matrices;
- sample-level intermediate objects;
- spreadsheets containing individual-level observations;
- unreviewed rendered reports;
- R caches or serialized workspaces containing restricted data;
- access credentials, tokens, or local configuration secrets.

A suitable `.gitignore` should exclude restricted data and derived sample-level files. At minimum, consider:

```gitignore
# Restricted data
Analyisis/
Analysis_Polars_Semipolars/
Analysis_Serum_Polars_Semipolars/
Analysis_PFAS/
Transcriptomics/
data/
raw_data/
private_data/

# R state and caches
.Rhistory
.RData
.Ruserdata
.Rproj.user/
*_cache/
*_files/

# Potentially sensitive outputs
results/
*.rds
*.RDS
*.xlsx
*.csv

# Operating-system and temporary files
.DS_Store
Thumbs.db
~$*
```

Review every exception before using `git add -f`.

</details>

---

## ✅ Reproducibility and quality control

The workflow includes reproducibility controls designed for the reported analyses:

- explicit package and software versions;
- fixed random seeds for stochastic analyses;
- validation of sample overlap and uniqueness;
- layer-specific accounting for omics availability;
- recalculation checks against observed association counts;
- raw resampling-distribution export where enabled;
- multiple-testing correction within defined analysis families;
- separation of public code from restricted data.

**Recommended additions before creating an archival release:**

- `sessionInfo.txt`
- `renv.lock`
- `CITATION.cff`
- a tagged GitHub release
- a permanent Zenodo archive and DOI
- a concise changelog documenting the manuscript-associated release

---

## ⚠ Responsible interpretation

> This workflow supports exploratory and discovery-oriented multi-omics analyses. Results should be interpreted in the context of the study design, unequal group sizes, layer-specific missingness, multiple testing, and the availability of independent validation.

Sex-stratified findings describe associations in the analyzed cohort. The size-matched sensitivity analysis evaluates the number of detected module-histology associations after accounting for sample size and layer-specific availability. **The analysis does not, by itself, establish an intrinsic biological difference in network architecture.**

### Code availability statement

The complete analysis workflow is publicly available in this repository:
[MASLD sex-stratified multi-omics analysis repository](https://github.com/shashank-KU/masld-sex-stratified-multiomics)

For publication, cite a versioned archival DOI once the manuscript-associated release has been deposited in Zenodo. The GitHub repository can remain the actively maintained version.

---

## 📄 Citation

Until the associated manuscript and archival DOI are available, use:

```text
Gupta S. MASLD Sex-Stratified Multi-Omics Analysis. GitHub repository.
https://github.com/shashank-KU/masld-sex-stratified-multiomics
```

After publication, replace the temporary citation with the full manuscript citation and archived software DOI. A `CITATION.cff` file should be added at that stage so GitHub can display the preferred citation.

---

## 📜 License

This repository is distributed under the [MIT License](LICENSE). The license applies to repository code and documentation only. It does not grant access to, or permission to redistribute, restricted participant-level data.

---

## 👤 Author

**Shashank Gupta** · Örebro University

<p align="left">
  <a href="https://github.com/shashank-KU"><img alt="GitHub" src="https://img.shields.io/badge/GitHub-shashank--KU-181717?logo=github&logoColor=white" /></a>
  <a href="https://www.oru.se/english/employee/shashank_gupta"><img alt="Örebro University" src="https://img.shields.io/badge/%C3%96rebro%20University-profile-003D7C" /></a>
</p>

---

## 📧 Contact

For questions about the analysis workflow, open a [GitHub issue](https://github.com/shashank-KU/masld-sex-stratified-multiomics/issues). Do not include participant-level information, restricted data, credentials, or confidential study details in public issues.

---

## 🙏 Acknowledgment

When reusing or adapting this workflow, cite the associated manuscript and the archived software release when available. Please also cite the original software publications listed in the manuscript Key Resources Table.
