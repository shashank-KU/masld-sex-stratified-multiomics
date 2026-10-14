# data/

Input files are not included in this repository. The participant-level data cannot be shared because the consent obtained does not cover public release.

To run the code, place your own copies of the files below in this folder (or set `MASLD_DATA_DIR` to the folder that contains them). Folder names are written exactly as the code reads them.

| Path (relative to `data/`) | Content |
|---|---|
| `Analyisis/Metadata.xlsx` | clinical and histology data |
| `Analyisis/exposures.xlsx` | hepatic lipidomics |
| `Analyisis/description.xlsx` | lipid annotation |
| `Analysis_Polars_Semipolars/exposure_polars_semipolars.xlsx` | hepatic polar metabolomics |
| `Analysis_Serum_Polars_Semipolars/` | serum polar metabolomics |
| `Transcriptomics/MAFALDA_STAR/genes.expected_count.txt` | liver RNA-seq counts |
| `PRS_analysis/pprs.xlsx` | partitioned polygenic risk scores |
| `Analyisis/GSE135251_raw_counts_GRCh38.p13_NCBI.tsv.gz` | GEO GSE135251 counts (public) |
| `Analyisis/ST001711_AN002786.txt` | EPoS serum metabolomics (public, Metabolomics Workbench) |

Public datasets are downloaded by the code or can be obtained from GEO (GSE162694, GSE135251) and Metabolomics Workbench (ST000915, ST001711).

The Supplementary Tables workbook (`gupta_MASLD_SupplementaryTables*.xlsx`) is optional: if absent, the code rebuilds it from the `Table_S*.csv` files in each figure's `Source_data/` folder.
