# Tables for Fig. 1 (and Fig. S1): hepatic transcriptome

First row of each CSV = column headers (caption removed; captions below).

**Table S1** (`tables/Table_S1.csv`)

Table S1. Differentially expressed genes associated with steatosis, lobular inflammation and fibrosis (grade/stage 0 vs ≥1) in 208 participants (153 women, 55 men). DESeq2 Wald test with age and BMI (standardized) and sex as covariates; log2FoldChange is apeglm-shrunk; padj, Benjamini–Hochberg FDR (≤ 0.05 shown). Steatosis 1,532 (987 up, 545 down); inflammation 1,852 (1,094 up, 758 down); fibrosis 93 (65 up, 28 down).

**Table S2a** (`tables/Table_S2a.csv`)

Table S2a. Overlap of differentially expressed genes (FDR ≤ 0.05; Table S1) across steatosis, lobular inflammation and fibrosis (1 = significant for that feature). Upregulated: 257 steatosis only, 364 inflammation only, 12 fibrosis only, 680 steatosis + inflammation, 3 steatosis + fibrosis, 3 inflammation + fibrosis, 47 all three. Downregulated: 87, 299, 10, 442, 1, 2 and 15, respectively (Fig. 1b).

**Table S2b** (`tables/Table_S2b.csv`)

Table S2b. KEGG pathway enrichment of up- and downregulated DEGs (Table S1) for each histological feature, tested with clusterProfiler against all tested genes as background (Benjamini–Hochberg FDR ≤ 0.10). Fig. 1c shows the top 10 pathways per feature and direction with ≥ 3 genes.
