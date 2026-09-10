# ============================================================
# FINAL EXTERNAL VALIDATION PIPELINE
# MASLD / NAFLD MULTI-OMICS
#
# ============================================================
#
# EXTERNAL SERUM DATA
#
# ST001711 / AN002786
#   Polar metabolomics
#   36 metabolites
#   627 samples
#
# ST001710 / AN002785
#   Lipidomics
#   176 lipids
#   627 samples
#
# EXTERNAL CONTRAST
#
#   Grade 1 vs Grades 2+3+4
#
# ============================================================
#
# OUTPUT
#
# Table_S11  = Transcriptomic validation summary
# Table_S12A = Transcriptomic gene-level validation
# Table_S12B = KEGG/pathway validation
# Table_S13A = Serum metabolomics summary
# Table_S13B = Serum metabolomics feature-level validation
# Table_S14A = Serum lipidomics summary
# Table_S14B = Serum lipidomics feature-level validation
#
# ============================================================


# ============================================================
# 0. CLEAN START
# ============================================================

rm(
  list = c(
    "epos_metabo",
    "epos_lipid",
    "metabo_matrix",
    "pheno_metabo",
    "lipid_matrix",
    "pheno_lipid",
    "res_metabo_steatosis",
    "res_metabo_inflammation",
    "res_metabo_fibrosis",
    "res_lipid_steatosis",
    "res_lipid_inflammation",
    "res_lipid_fibrosis",
    "disc_metabo",
    "metabo_feature_table",
    "disc_lipid",
    "disc_lipid_match"
  ),
  envir = .GlobalEnv
)


# ============================================================
# 1. PACKAGES
# ============================================================

required_packages <- c(
  "openxlsx",
  "limma",
  "dplyr",
  "stringr",
  "tibble"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    install.packages(pkg)
    
  }
  
}

library(openxlsx)
library(limma)
library(dplyr)
library(stringr)
library(tibble)


# ============================================================
# 2. SETTINGS
# ============================================================

FDR_OMICS <- 0.10
FDR_RNA   <- 0.05

METABO_FILE <- "ST001711_AN002786.txt"
LIPID_FILE  <- "ST001710_AN002785.txt"

OUTPUT_FILE <-
  "Supplementary_External_Validation_Tables_FINAL.xlsx"


# ============================================================
# 3. BASIC HELPER FUNCTIONS
# ============================================================

safe_numeric <- function(x) {
  
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
  
}


normalise_name <- function(x) {
  
  x <- toupper(
    trimws(
      as.character(x)
    )
  )
  
  gsub(
    "[^A-Z0-9]",
    "",
    x
  )
  
}


first_existing_column <- function(
    df,
    candidates
) {
  
  hits <- candidates[
    candidates %in% names(df)
  ]
  
  if (length(hits) == 0) {
    
    return(NULL)
    
  }
  
  hits[1]
  
}


# ============================================================
# 4. CHECK INPUT FILES
# ============================================================

if (!file.exists(METABO_FILE)) {
  
  stop(
    "\nMetabolomics file not found:\n",
    METABO_FILE,
    "\n\nPlace ST001711_AN002786.txt in the current R working directory."
  )
  
}


if (!file.exists(LIPID_FILE)) {
  
  stop(
    "\nLipidomics file not found:\n",
    LIPID_FILE,
    "\n\nPlace ST001710_AN002785.txt in the current R working directory."
  )
  
}


cat(
  "\nCurrent working directory:\n",
  getwd(),
  "\n"
)


# ============================================================
# 5. CHECK SERUM DISCOVERY OBJECTS
# ============================================================

required_serum_objects <- c(
  "res_steatosis",
  "res_inflammation",
  "res_fibrosis",
  "my_results"
)

missing_serum_objects <-
  required_serum_objects[
    !vapply(
      required_serum_objects,
      exists,
      logical(1),
      inherits = TRUE
    )
  ]


if (
  length(missing_serum_objects) > 0
) {
  
  stop(
    "\nMissing serum discovery objects:\n\n",
    paste(
      missing_serum_objects,
      collapse = "\n"
    ),
    "\n\nRun your final serum discovery analyses first."
  )
  
}


# ============================================================
# 6. EPoS MWTAB IMPORTER
# ============================================================
#
# Important:
#
# Main data block:
#
#   Samples
#   Factors
#   feature rows...
#
# SUBJECT_SAMPLE_FACTORS:
#
#   x[1] = SUBJECT_SAMPLE_FACTORS
#   x[2] = -
#   x[3] = SAMPLE ID
#   x[4] = FACTORS
#   x[5+] = CLINICAL DATA
#
# Metadata can use:
#
#   Kleiner.Steatosis:1
#   Inflammation:1
#   Kleiner.Fibrosis=1a
#   NAS=4
#
# ============================================================

read_epos_mwtab <- function(file) {
  
  cat(
    "\nReading:",
    file,
    "\n"
  )
  
  
  lines <- readLines(
    file,
    encoding = "UTF-8",
    warn = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Locate data block
  # ----------------------------------------------------------
  
  start_idx <- which(
    trimws(lines) ==
      "MS_METABOLITE_DATA_START"
  )
  
  end_idx <- which(
    trimws(lines) ==
      "MS_METABOLITE_DATA_END"
  )
  
  
  if (
    length(start_idx) != 1 ||
    length(end_idx) != 1 ||
    end_idx <= start_idx
  ) {
    
    stop(
      "\nCould not locate MS_METABOLITE_DATA_START/END in ",
      file
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Samples header
  # ----------------------------------------------------------
  
  header <- strsplit(
    lines[start_idx + 1],
    "\t",
    fixed = TRUE
  )[[1]]
  
  
  if (
    length(header) < 2 ||
    trimws(header[1]) != "Samples"
  ) {
    
    stop(
      "\nUnexpected Samples header in ",
      file
    )
    
  }
  
  
  sample_ids <- trimws(
    header[-1]
  )
  
  
  # ----------------------------------------------------------
  # Features
  #
  # start + 1 = Samples
  # start + 2 = Factors
  # start + 3 = first feature
  # ----------------------------------------------------------
  
  data_lines <- lines[
    (start_idx + 3):(end_idx - 1)
  ]
  
  
  data_lines <- data_lines[
    nzchar(
      trimws(data_lines)
    )
  ]
  
  
  if (
    length(data_lines) == 0
  ) {
    
    stop(
      "\nNo feature rows found in ",
      file
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Parse rows
  # ----------------------------------------------------------
  
  parsed <- lapply(
    data_lines,
    function(z) {
      
      x <- strsplit(
        z,
        "\t",
        fixed = TRUE
      )[[1]]
      
      
      if (
        length(x) < length(header)
      ) {
        
        x <- c(
          x,
          rep(
            NA_character_,
            length(header) -
              length(x)
          )
        )
        
      }
      
      
      if (
        length(x) > length(header)
      ) {
        
        x <- x[
          seq_len(
            length(header)
          )
        ]
        
      }
      
      
      x
      
    }
  )
  
  
  # ----------------------------------------------------------
  # Feature names
  # ----------------------------------------------------------
  
  feature_names <- vapply(
    parsed,
    function(z) {
      
      trimws(
        z[1]
      )
      
    },
    character(1)
  )
  
  
  if (
    any(
      is.na(feature_names) |
      feature_names == ""
    )
  ) {
    
    stop(
      "\nMissing feature name in ",
      file
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Numeric matrix
  # ----------------------------------------------------------
  
  value_mat <- matrix(
    NA_real_,
    nrow = length(parsed),
    ncol = length(sample_ids)
  )
  
  
  for (i in seq_along(parsed)) {
    
    value_mat[i, ] <- suppressWarnings(
      as.numeric(
        parsed[[i]][-1]
      )
    )
    
  }
  
  
  rownames(value_mat) <- make.unique(
    feature_names
  )
  
  colnames(value_mat) <- sample_ids
  
  
  # ==========================================================
  # SAMPLE METADATA
  # ==========================================================
  
  factor_idx <- which(
    grepl(
      "^SUBJECT_SAMPLE_FACTORS\\s",
      lines
    )
  )
  
  
  if (
    length(factor_idx) == 0
  ) {
    
    stop(
      "\nNo SUBJECT_SAMPLE_FACTORS records found in ",
      file
    )
    
  }
  
  
  metadata <- tibble(
    SampleID = sample_ids,
    Factors = NA_character_,
    Clinical = NA_character_
  )
  
  
  for (idx in factor_idx) {
    
    x <- strsplit(
      lines[idx],
      "\t",
      fixed = TRUE
    )[[1]]
    
    
    if (
      length(x) < 4
    ) {
      
      next
      
    }
    
    
    # IMPORTANT:
    # Sample ID = x[3]
    sid <- trimws(
      x[3]
    )
    
    
    if (
      is.na(sid) ||
      !sid %in% metadata$SampleID
    ) {
      
      next
      
    }
    
    
    pos <- match(
      sid,
      metadata$SampleID
    )
    
    
    metadata$Factors[pos] <-
      trimws(
        x[4]
      )
    
    
    if (
      length(x) >= 5
    ) {
      
      metadata$Clinical[pos] <-
        paste(
          x[5:length(x)],
          collapse = " | "
        )
      
    }
    
  }
  
  
  # ==========================================================
  # FIELD EXTRACTION
  # ==========================================================
  
  extract_field <- function(
    factors,
    clinical,
    key
  ) {
    
    txt <- paste(
      ifelse(
        is.na(factors),
        "",
        factors
      ),
      ifelse(
        is.na(clinical),
        "",
        clinical
      ),
      sep = " | "
    )
    
    
    pattern <- paste0(
      "(?:^|[|;[:space:]])",
      key,
      "\\s*[:=]\\s*",
      "([^|;[:space:]]+)"
    )
    
    
    out <- stringr::str_match(
      txt,
      pattern
    )[, 2]
    
    
    out[
      is.na(out) |
        out %in% c(
          "",
          "-",
          "NA",
          "NaN"
        )
    ] <- NA_character_
    
    
    out
    
  }
  
  
  # ----------------------------------------------------------
  # Raw fields
  # ----------------------------------------------------------
  
  metadata$Steatosis_raw <-
    extract_field(
      metadata$Factors,
      metadata$Clinical,
      "Kleiner\\.Steatosis"
    )
  
  
  metadata$Inflammation_raw <-
    extract_field(
      metadata$Factors,
      metadata$Clinical,
      "Inflammation"
    )
  
  
  metadata$Fibrosis_raw <-
    extract_field(
      metadata$Factors,
      metadata$Clinical,
      "Kleiner\\.Fibrosis"
    )
  
  
  metadata$NAS_raw <-
    extract_field(
      metadata$Factors,
      metadata$Clinical,
      "NAS"
    )
  
  
  # ----------------------------------------------------------
  # Numeric variables
  # ----------------------------------------------------------
  
  metadata$Steatosis <-
    safe_numeric(
      metadata$Steatosis_raw
    )
  
  
  metadata$Inflammation <-
    safe_numeric(
      metadata$Inflammation_raw
    )
  
  
  metadata$NAS <-
    safe_numeric(
      metadata$NAS_raw
    )
  
  
  # ----------------------------------------------------------
  # Fibrosis mapping
  #
  # 1 / 1a / 1b / 1c -> 1
  # 2 -> 2
  # 3 -> 3
  # 4 -> 4
  # 6 -> NA
  # ----------------------------------------------------------
  
  metadata$Fibrosis <-
    dplyr::case_when(
      
      metadata$Fibrosis_raw %in%
        c(
          "1",
          "1a",
          "1b",
          "1c"
        ) ~ 1,
      
      metadata$Fibrosis_raw == "2" ~
        2,
      
      metadata$Fibrosis_raw == "3" ~
        3,
      
      metadata$Fibrosis_raw == "4" ~
        4,
      
      TRUE ~
        NA_real_
      
    )
  
  
  # ----------------------------------------------------------
  # Safety checks
  # ----------------------------------------------------------
  
  if (
    ncol(value_mat) != 627
  ) {
    
    stop(
      "\nExpected 627 samples in ",
      file,
      ", found ",
      ncol(value_mat)
    )
    
  }
  
  
  if (
    !identical(
      colnames(value_mat),
      metadata$SampleID
    )
  ) {
    
    stop(
      "\nMatrix and metadata sample IDs do not match in ",
      file
    )
    
  }
  
  
  list(
    matrix = value_mat,
    metadata = metadata
  )
  
}


# ============================================================
# 7. READ EPoS FILES
# ============================================================

epos_metabo <-
  read_epos_mwtab(
    METABO_FILE
  )

epos_lipid <-
  read_epos_mwtab(
    LIPID_FILE
  )


metabo_matrix <-
  epos_metabo$matrix

pheno_metabo <-
  epos_metabo$metadata


lipid_matrix <-
  epos_lipid$matrix

pheno_lipid <-
  epos_lipid$metadata


# ------------------------------------------------------------
# Dimensions
# ------------------------------------------------------------

cat(
  "\n============================================\n",
  "EPoS DIMENSIONS\n",
  "============================================\n",
  sep = ""
)


cat(
  "\nMetabolomics:",
  nrow(metabo_matrix),
  "features x",
  ncol(metabo_matrix),
  "samples\n"
)


cat(
  "Lipidomics:",
  nrow(lipid_matrix),
  "features x",
  ncol(lipid_matrix),
  "samples\n"
)


# HARD CHECKS

if (
  nrow(metabo_matrix) != 36
) {
  
  stop(
    "Metabolomics should contain 36 features; found ",
    nrow(metabo_matrix)
  )
  
}


if (
  nrow(lipid_matrix) != 176
) {
  
  stop(
    "Lipidomics should contain 176 features; found ",
    nrow(lipid_matrix)
  )
  
}


# ============================================================
# 8. METADATA CHECK
# ============================================================

cat(
  "\n============================================\n",
  "EPoS METADATA CHECK\n",
  "============================================\n",
  sep = ""
)


cat(
  "\nSTEATOSIS:\n"
)

print(
  table(
    pheno_metabo$Steatosis,
    useNA = "ifany"
  )
)


cat(
  "\nINFLAMMATION:\n"
)

print(
  table(
    pheno_metabo$Inflammation,
    useNA = "ifany"
  )
)


cat(
  "\nFIBROSIS:\n"
)

print(
  table(
    pheno_metabo$Fibrosis,
    useNA = "ifany"
  )
)


# ============================================================
# 9. HARD-CHECK EXPECTED ANALYSIS GROUPS
# ============================================================

# Steatosis: 212 vs 388
if (
  sum(
    pheno_metabo$Steatosis == 1,
    na.rm = TRUE
  ) != 212 ||
  sum(
    pheno_metabo$Steatosis %in%
    c(2, 3, 4),
    na.rm = TRUE
  ) != 388
) {
  
  stop(
    "\nSteatosis groups are not 212 / 388."
  )
  
}


# Inflammation: 307 vs 247
if (
  sum(
    pheno_metabo$Inflammation == 1,
    na.rm = TRUE
  ) != 307 ||
  sum(
    pheno_metabo$Inflammation %in%
    c(2, 3, 4),
    na.rm = TRUE
  ) != 247
) {
  
  stop(
    "\nInflammation groups are not 307 / 247."
  )
  
}


# Fibrosis: 157 vs 293
if (
  sum(
    pheno_metabo$Fibrosis == 1,
    na.rm = TRUE
  ) != 157 ||
  sum(
    pheno_metabo$Fibrosis %in%
    c(2, 3, 4),
    na.rm = TRUE
  ) != 293
) {
  
  stop(
    "\nFibrosis groups are not 157 / 293."
  )
  
}


cat(
  "\nEPoS metadata checks PASSED.\n"
)


# ============================================================
# 10. EXTERNAL BINARY limma
# ============================================================
#
# Grade 1 vs Grades 2+3+4
#
# EPoS data are already log2-transformed and scaled.
# No extra normalization is performed.
# ============================================================

run_external_binary <- function(
    matrix_data,
    phenotype,
    fdr = 0.10
) {
  
  if (
    length(phenotype) !=
    ncol(matrix_data)
  ) {
    
    stop(
      "\nPhenotype length does not match matrix sample number."
    )
    
  }
  
  
  grp <- rep(
    NA_character_,
    length(phenotype)
  )
  
  
  grp[
    phenotype == 1
  ] <- "Grade1"
  
  
  grp[
    phenotype %in%
      c(
        2,
        3,
        4
      )
  ] <- "Grade2_3_4"
  
  
  keep <- !is.na(
    grp
  )
  
  
  grp <- factor(
    grp[keep],
    levels = c(
      "Grade1",
      "Grade2_3_4"
    )
  )
  
  
  x <- matrix_data[
    ,
    keep,
    drop = FALSE
  ]
  
  
  if (
    any(
      table(grp) == 0
    )
  ) {
    
    stop(
      "\nOne analysis group contains zero samples."
    )
    
  }
  
  
  design <- model.matrix(
    ~ grp
  )
  
  
  fit <- limma::lmFit(
    x,
    design
  )
  
  
  fit <- limma::eBayes(
    fit
  )
  
  
  tt <- limma::topTable(
    fit,
    coef = "grpGrade2_3_4",
    number = Inf,
    adjust.method = "fdr",
    sort.by = "P"
  )
  
  
  tt$Feature <-
    rownames(
      tt
    )
  
  
  tt$Direction <-
    ifelse(
      tt$logFC > 0,
      "Up",
      ifelse(
        tt$logFC < 0,
        "Down",
        "Zero"
      )
    )
  
  
  tt$Significant <-
    !is.na(
      tt$adj.P.Val
    ) &
    tt$adj.P.Val < fdr
  
  
  list(
    result = tt,
    n_grade1 =
      sum(
        grp == "Grade1"
      ),
    n_grade234 =
      sum(
        grp == "Grade2_3_4"
      )
  )
  
}


# ============================================================
# 11. EXTERNAL METABOLOMICS
# ============================================================

cat(
  "\n============================================\n",
  "EXTERNAL SERUM METABOLOMICS\n",
  "============================================\n",
  sep = ""
)


res_metabo_steatosis <-
  run_external_binary(
    metabo_matrix,
    pheno_metabo$Steatosis,
    FDR_OMICS
  )


res_metabo_inflammation <-
  run_external_binary(
    metabo_matrix,
    pheno_metabo$Inflammation,
    FDR_OMICS
  )


res_metabo_fibrosis <-
  run_external_binary(
    metabo_matrix,
    pheno_metabo$Fibrosis,
    FDR_OMICS
  )


cat(
  "\nSteatosis:",
  res_metabo_steatosis$n_grade1,
  "vs",
  res_metabo_steatosis$n_grade234,
  "\n"
)

cat(
  "Significant:",
  sum(
    res_metabo_steatosis$result$Significant
  ),
  "\n"
)


cat(
  "\nInflammation:",
  res_metabo_inflammation$n_grade1,
  "vs",
  res_metabo_inflammation$n_grade234,
  "\n"
)

cat(
  "Significant:",
  sum(
    res_metabo_inflammation$result$Significant
  ),
  "\n"
)


cat(
  "\nFibrosis:",
  res_metabo_fibrosis$n_grade1,
  "vs",
  res_metabo_fibrosis$n_grade234,
  "\n"
)

cat(
  "Significant:",
  sum(
    res_metabo_fibrosis$result$Significant
  ),
  "\n"
)


# Expected:
#
# Steatosis      0
# Inflammation   6
# Fibrosis       0


# ============================================================
# 12. DISCOVERY METABOLOMICS
# ============================================================
#
# CURRENT FINAL WORKSPACE RESULTS:
#
# Steatosis      12
# Inflammation    4
# Fibrosis        0
#
# ============================================================

flatten_metabo <- function(
    x,
    phenotype,
    fdr = 0.10
) {
  
  x <- as.data.frame(
    x
  )
  
  
  if (
    !all(
      c(
        "logFC",
        "adj.P.Val"
      ) %in%
      names(x)
    )
  ) {
    
    stop(
      "\nDiscovery metabolomics result does not contain ",
      "logFC and adj.P.Val."
    )
    
  }
  
  
  sig_idx <- which(
    !is.na(
      x$adj.P.Val
    ) &
      x$adj.P.Val < fdr
  )
  
  
  feature <- rownames(
    x
  )
  
  
  if (
    "Feature" %in%
    names(x)
  ) {
    
    feature <- as.character(
      x$Feature
    )
    
  }
  
  
  if (
    "Metabolite" %in%
    names(x)
  ) {
    
    feature <- as.character(
      x$Metabolite
    )
    
  }
  
  
  feature <- feature[
    sig_idx
  ]
  
  
  out <- data.frame(
    
    Feature =
      as.character(
        feature
      ),
    
    Discovery_logFC =
      as.numeric(
        x$logFC[sig_idx]
      ),
    
    Discovery_FDR =
      as.numeric(
        x$adj.P.Val[sig_idx]
      ),
    
    Phenotype =
      rep(
        as.character(phenotype),
        length(sig_idx)
      ),
    
    stringsAsFactors = FALSE
    
  )
  
  
  out$Discovery_Direction <-
    ifelse(
      out$Discovery_logFC > 0,
      "Up",
      ifelse(
        out$Discovery_logFC < 0,
        "Down",
        "Zero"
      )
    )
  
  
  out$Feature_norm <-
    normalise_name(
      out$Feature
    )
  
  
  out
  
}


disc_steatosis <-
  flatten_metabo(
    res_steatosis,
    "Steatosis",
    FDR_OMICS
  )


disc_inflammation <-
  flatten_metabo(
    res_inflammation,
    "Inflammation",
    FDR_OMICS
  )


disc_fibrosis <-
  flatten_metabo(
    res_fibrosis,
    "Fibrosis",
    FDR_OMICS
  )


disc_metabo <- rbind(
  disc_steatosis,
  disc_inflammation,
  disc_fibrosis
)


rownames(
  disc_metabo
) <- NULL


# ------------------------------------------------------------
# Hard check current discovery results
# ------------------------------------------------------------

cat(
  "\n============================================\n",
  "DISCOVERY SERUM METABOLOMICS\n",
  "============================================\n",
  sep = ""
)


cat(
  "\nSteatosis:",
  nrow(disc_steatosis),
  "\n"
)

cat(
  "Inflammation:",
  nrow(disc_inflammation),
  "\n"
)

cat(
  "Fibrosis:",
  nrow(disc_fibrosis),
  "\n"
)


if (
  nrow(disc_steatosis) != 12
) {
  
  stop(
    "\nCurrent res_steatosis contains ",
    nrow(disc_steatosis),
    " significant metabolites, not 12."
  )
  
}


if (
  nrow(disc_inflammation) != 4
) {
  
  stop(
    "\nCurrent res_inflammation contains ",
    nrow(disc_inflammation),
    " significant metabolites, not 4."
  )
  
}


if (
  nrow(disc_fibrosis) != 0
) {
  
  stop(
    "\nCurrent res_fibrosis contains ",
    nrow(disc_fibrosis),
    " significant metabolites, not 0."
  )
  
}


# ============================================================
# 13. METABOLOMICS FEATURE MATCHING
# ============================================================

external_metabo_lookup <- data.frame(
  
  External_feature =
    as.character(
      rownames(
        metabo_matrix
      )
    ),
  
  stringsAsFactors = FALSE
  
)


external_metabo_lookup$Feature_norm <-
  normalise_name(
    external_metabo_lookup$External_feature
  )


# ------------------------------------------------------------
# Exact aliases
# ------------------------------------------------------------

disc_metabo$Feature_norm[
  disc_metabo$Feature ==
    "GLUTAMIC_ACID"
] <- "GLUTAMICACID"


disc_metabo$Feature_norm[
  disc_metabo$Feature ==
    "ISOLEUCINE"
] <- "ISOLEUCINE"


# ------------------------------------------------------------
# Match discovery metabolites to external measured features
# ------------------------------------------------------------

metabo_feature_table <- merge(
  
  disc_metabo,
  
  external_metabo_lookup,
  
  by = "Feature_norm",
  
  all.x = TRUE,
  
  sort = FALSE
  
)


metabo_feature_table$Measured_externally <-
  !is.na(
    metabo_feature_table$External_feature
  )


# ============================================================
# 14. ADD EXTERNAL METABOLOMICS STATISTICS
# ============================================================

external_metabo_results <- rbind(
  
  data.frame(
    
    Phenotype =
      "Steatosis",
    
    External_feature =
      as.character(
        res_metabo_steatosis$result$Feature
      ),
    
    External_logFC =
      res_metabo_steatosis$result$logFC,
    
    External_FDR =
      res_metabo_steatosis$result$adj.P.Val,
    
    stringsAsFactors = FALSE
    
  ),
  
  data.frame(
    
    Phenotype =
      "Inflammation",
    
    External_feature =
      as.character(
        res_metabo_inflammation$result$Feature
      ),
    
    External_logFC =
      res_metabo_inflammation$result$logFC,
    
    External_FDR =
      res_metabo_inflammation$result$adj.P.Val,
    
    stringsAsFactors = FALSE
    
  ),
  
  data.frame(
    
    Phenotype =
      "Fibrosis",
    
    External_feature =
      as.character(
        res_metabo_fibrosis$result$Feature
      ),
    
    External_logFC =
      res_metabo_fibrosis$result$logFC,
    
    External_FDR =
      res_metabo_fibrosis$result$adj.P.Val,
    
    stringsAsFactors = FALSE
    
  )
  
)


# ------------------------------------------------------------
# Match external statistics
# ------------------------------------------------------------

metabo_feature_table <- merge(
  
  metabo_feature_table,
  
  external_metabo_results,
  
  by = c(
    "Phenotype",
    "External_feature"
  ),
  
  all.x = TRUE,
  
  sort = FALSE
  
)


metabo_feature_table$External_Direction <-
  ifelse(
    
    metabo_feature_table$External_logFC > 0,
    
    "Up",
    
    ifelse(
      
      metabo_feature_table$External_logFC < 0,
      
      "Down",
      
      NA_character_
      
    )
    
  )


metabo_feature_table$Direction_concordant <-
  ifelse(
    
    !metabo_feature_table$Measured_externally,
    
    NA,
    
    metabo_feature_table$Discovery_Direction ==
      metabo_feature_table$External_Direction
    
  )


# ============================================================
# 15. TABLE S13A
# ============================================================

make_metabo_summary <- function(
    phenotype,
    feature_table,
    external_result
) {
  
  z <- feature_table[
    feature_table$Phenotype ==
      phenotype,
    ,
    drop = FALSE
  ]
  
  
  measured <- sum(
    z$Measured_externally,
    na.rm = TRUE
  )
  
  
  same <- sum(
    z$Direction_concordant %in%
      TRUE
  )
  
  
  opposite <- sum(
    z$Direction_concordant %in%
      FALSE
  )
  
  
  not_measured <- sum(
    !z$Measured_externally
  )
  
  
  data.frame(
    
    Phenotype =
      phenotype,
    
    Contrast =
      "Grade 1 vs Grades 2+3+4",
    
    Discovery_significant =
      nrow(z),
    
    External_significant =
      sum(
        external_result$result$Significant
      ),
    
    Measured_externally =
      measured,
    
    Same_direction =
      same,
    
    Opposite_direction =
      opposite,
    
    Not_measured =
      not_measured,
    
    Percent_concordant =
      if (
        measured > 0
      ) {
        
        100 *
          same /
          measured
        
      } else {
        
        NA_real_
        
      },
    
    stringsAsFactors = FALSE
    
  )
  
}


Table_S13A <- rbind(
  
  make_metabo_summary(
    "Steatosis",
    metabo_feature_table,
    res_metabo_steatosis
  ),
  
  make_metabo_summary(
    "Inflammation",
    metabo_feature_table,
    res_metabo_inflammation
  ),
  
  make_metabo_summary(
    "Fibrosis",
    metabo_feature_table,
    res_metabo_fibrosis
  )
  
)


# ============================================================
# 16. TABLE S13B
# ============================================================

Table_S13B <- metabo_feature_table[
  ,
  c(
    
    "Phenotype",
    
    "Feature",
    
    "External_feature",
    
    "Discovery_logFC",
    
    "Discovery_FDR",
    
    "Discovery_Direction",
    
    "External_logFC",
    
    "External_FDR",
    
    "External_Direction",
    
    "Measured_externally",
    
    "Direction_concordant"
    
  ),
  drop = FALSE
]


Table_S13B <- Table_S13B[
  order(
    Table_S13B$Phenotype,
    Table_S13B$Feature
  ),
  ,
  drop = FALSE
]


cat(
  "\n============================================\n",
  "TABLE S13A\n",
  "============================================\n",
  sep = ""
)

print(
  Table_S13A,
  row.names = FALSE
)


# ============================================================
# 17. LIPIDOMICS EXTERNAL ANALYSIS
# ============================================================

cat(
  "\n============================================\n",
  "EXTERNAL SERUM LIPIDOMICS\n",
  "============================================\n",
  sep = ""
)


res_lipid_steatosis <-
  run_external_binary(
    lipid_matrix,
    pheno_lipid$Steatosis,
    FDR_OMICS
  )


res_lipid_inflammation <-
  run_external_binary(
    lipid_matrix,
    pheno_lipid$Inflammation,
    FDR_OMICS
  )


res_lipid_fibrosis <-
  run_external_binary(
    lipid_matrix,
    pheno_lipid$Fibrosis,
    FDR_OMICS
  )


cat(
  "\nSteatosis:",
  res_lipid_steatosis$n_grade1,
  "vs",
  res_lipid_steatosis$n_grade234,
  "\n"
)

cat(
  "Significant:",
  sum(
    res_lipid_steatosis$result$Significant
  ),
  "\n"
)


cat(
  "\nInflammation:",
  res_lipid_inflammation$n_grade1,
  "vs",
  res_lipid_inflammation$n_grade234,
  "\n"
)

cat(
  "Significant:",
  sum(
    res_lipid_inflammation$result$Significant
  ),
  "\n"
)


cat(
  "\nFibrosis:",
  res_lipid_fibrosis$n_grade1,
  "vs",
  res_lipid_fibrosis$n_grade234,
  "\n"
)

cat(
  "Significant:",
  sum(
    res_lipid_fibrosis$result$Significant
  ),
  "\n"
)


# Expected:
#
# Steatosis      21
# Inflammation    0
# Fibrosis       33


# ============================================================
# CORRECTED LIPIDOMICS ANALYSIS
#
# Multiple Unique Identifiers can share the same
# "Lipid main ID".
#
# Biological feature-level validation is therefore performed
# at the unique Phenotype + Lipid main ID level.
#
# For duplicate discovery records, retain the record with
# the smallest FDR.
# ============================================================


# ============================================================
# 18. DISCOVERY LIPIDOMICS
# ============================================================

my_results_df <- as.data.frame(
  my_results
)


# ------------------------------------------------------------
# Make key columns atomic
# ------------------------------------------------------------

my_results_df$Phenotype <-
  as.character(
    my_results_df$Phenotype
  )

my_results_df$`Lipid main ID` <-
  as.character(
    my_results_df$`Lipid main ID`
  )

my_results_df$logFC <-
  safe_numeric(
    my_results_df$logFC
  )

my_results_df$adj.P.Val <-
  safe_numeric(
    my_results_df$adj.P.Val
  )


# ============================================================
# 19. SIGNIFICANT DISCOVERY LIPIDS
# ============================================================

disc_lipid_raw <- my_results_df[
  !is.na(
    my_results_df$adj.P.Val
  ) &
    my_results_df$adj.P.Val < FDR_OMICS,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# Sort by:
#   Phenotype
#   Lipid main ID
#   lowest FDR first
# ------------------------------------------------------------

disc_lipid_raw <- disc_lipid_raw[
  order(
    disc_lipid_raw$Phenotype,
    disc_lipid_raw$`Lipid main ID`,
    disc_lipid_raw$adj.P.Val
  ),
  ,
  drop = FALSE
]


# ============================================================
# 20. COLLAPSE DUPLICATE MAIN LIPID IDs
# ============================================================
#
# Keep the row with the smallest FDR for each:
#
#   Phenotype + Lipid main ID
#
# ============================================================

duplicate_key <- paste(
  disc_lipid_raw$Phenotype,
  disc_lipid_raw$`Lipid main ID`,
  sep = "|||"
)


disc_lipid <- disc_lipid_raw[
  !duplicated(
    duplicate_key
  ),
  ,
  drop = FALSE
]


rownames(
  disc_lipid
) <- NULL


# ------------------------------------------------------------
# Standardized columns
# ------------------------------------------------------------

disc_lipid$Discovery_logFC <-
  disc_lipid$logFC


disc_lipid$Discovery_FDR <-
  disc_lipid$adj.P.Val


disc_lipid$Discovery_Direction <-
  ifelse(
    
    disc_lipid$Discovery_logFC > 0,
    
    "Up",
    
    ifelse(
      
      disc_lipid$Discovery_logFC < 0,
      
      "Down",
      
      "Zero"
      
    )
    
  )


# ============================================================
# 21. CHECK DISCOVERY COUNTS
# ============================================================

cat(
  "\n============================================\n",
  "FINAL DISCOVERY LIPID COUNTS\n",
  "============================================\n",
  sep = ""
)


for (
  ph in c(
    "Steatosis",
    "Inflammation",
    "Fibrosis"
  )
) {
  
  z <- disc_lipid[
    disc_lipid$Phenotype == ph,
    ,
    drop = FALSE
  ]
  
  cat(
    ph,
    ":",
    nrow(z),
    "\n"
  )
  
}


# Expected:
#
# Steatosis      826
# Inflammation   903
# Fibrosis       707


if (
  sum(
    disc_lipid$Phenotype ==
    "Steatosis"
  ) != 826
) {
  
  stop(
    "\nExpected 826 unique significant steatosis Lipid main IDs."
  )
  
}


if (
  sum(
    disc_lipid$Phenotype ==
    "Inflammation"
  ) != 903
) {
  
  stop(
    "\nExpected 903 unique significant inflammation Lipid main IDs."
  )
  
}


if (
  sum(
    disc_lipid$Phenotype ==
    "Fibrosis"
  ) != 707
) {
  
  stop(
    "\nExpected 707 unique significant fibrosis Lipid main IDs."
  )
  
}


cat(
  "\nDiscovery lipid counts PASSED.\n"
)


# ============================================================
# 22. EXTERNAL LIPIDOMICS RESULTS
# ============================================================

external_lipid_results <- rbind(
  
  data.frame(
    
    Phenotype =
      "Steatosis",
    
    External_feature =
      as.character(
        res_lipid_steatosis$result$Feature
      ),
    
    External_logFC =
      res_lipid_steatosis$result$logFC,
    
    External_FDR =
      res_lipid_steatosis$result$adj.P.Val,
    
    stringsAsFactors = FALSE
    
  ),
  
  data.frame(
    
    Phenotype =
      "Inflammation",
    
    External_feature =
      as.character(
        res_lipid_inflammation$result$Feature
      ),
    
    External_logFC =
      res_lipid_inflammation$result$logFC,
    
    External_FDR =
      res_lipid_inflammation$result$adj.P.Val,
    
    stringsAsFactors = FALSE
    
  ),
  
  data.frame(
    
    Phenotype =
      "Fibrosis",
    
    External_feature =
      as.character(
        res_lipid_fibrosis$result$Feature
      ),
    
    External_logFC =
      res_lipid_fibrosis$result$logFC,
    
    External_FDR =
      res_lipid_fibrosis$result$adj.P.Val,
    
    stringsAsFactors = FALSE
    
  )
  
)


# ------------------------------------------------------------
# Ensure external features are unique
# ------------------------------------------------------------

external_key <- paste(
  
  external_lipid_results$Phenotype,
  
  external_lipid_results$External_feature,
  
  sep = "|||"
  
)


if (
  any(
    duplicated(
      external_key
    )
  )
) {
  
  stop(
    "\nDuplicate external lipid feature names detected."
  )
  
}


# ============================================================
# 23. DISCOVERY → EXTERNAL MATCHING
# ============================================================

disc_lipid_match <- merge(
  
  disc_lipid,
  
  external_lipid_results,
  
  by.x = c(
    "Phenotype",
    "Lipid main ID"
  ),
  
  by.y = c(
    "Phenotype",
    "External_feature"
  ),
  
  all.x = TRUE,
  
  sort = FALSE
  
)


# ------------------------------------------------------------
# Measured externally
# ------------------------------------------------------------

disc_lipid_match$Measured_externally <-
  !is.na(
    disc_lipid_match$External_logFC
  )


# ------------------------------------------------------------
# External direction
# ------------------------------------------------------------

disc_lipid_match$External_Direction <-
  ifelse(
    
    disc_lipid_match$External_logFC > 0,
    
    "Up",
    
    ifelse(
      
      disc_lipid_match$External_logFC < 0,
      
      "Down",
      
      NA_character_
      
    )
    
  )


# ------------------------------------------------------------
# Directional concordance
# ------------------------------------------------------------

disc_lipid_match$Direction_concordant <-
  ifelse(
    
    !disc_lipid_match$Measured_externally,
    
    NA,
    
    disc_lipid_match$Discovery_Direction ==
      disc_lipid_match$External_Direction
    
  )


# ============================================================
# 24. TABLE S14A
# ============================================================

make_lipid_summary <- function(
    
  phenotype,
  
  df,
  
  external_result
  
) {
  
  z <- df[
    
    df$Phenotype ==
      phenotype,
    
    ,
    
    drop = FALSE
    
  ]
  
  
  measured <- sum(
    
    z$Measured_externally,
    
    na.rm = TRUE
    
  )
  
  
  same <- sum(
    
    z$Direction_concordant %in%
      TRUE
    
  )
  
  
  opposite <- sum(
    
    z$Direction_concordant %in%
      FALSE
    
  )
  
  
  data.frame(
    
    Phenotype =
      phenotype,
    
    Contrast =
      "Grade 1 vs Grades 2+3+4",
    
    Discovery_significant =
      nrow(z),
    
    External_significant =
      sum(
        external_result$result$Significant
      ),
    
    Measured_externally =
      measured,
    
    Same_direction =
      same,
    
    Opposite_direction =
      opposite,
    
    Percent_concordant =
      
      if (
        measured > 0
      ) {
        
        100 *
          same /
          measured
        
      } else {
        
        NA_real_
        
      },
    
    stringsAsFactors = FALSE
    
  )
  
}


Table_S14A <- rbind(
  
  make_lipid_summary(
    "Steatosis",
    disc_lipid_match,
    res_lipid_steatosis
  ),
  
  make_lipid_summary(
    "Inflammation",
    disc_lipid_match,
    res_lipid_inflammation
  ),
  
  make_lipid_summary(
    "Fibrosis",
    disc_lipid_match,
    res_lipid_fibrosis
  )
  
)


# ============================================================
# 25. TABLE S14B
# ============================================================

Table_S14B <- disc_lipid_match[
  
  ,
  
  c(
    
    "Phenotype",
    
    "Lipid main ID",
    
    "Lipid Class",
    
    "Discovery_logFC",
    
    "Discovery_FDR",
    
    "Discovery_Direction",
    
    "External_logFC",
    
    "External_FDR",
    
    "External_Direction",
    
    "Measured_externally",
    
    "Direction_concordant",
    
    "Unique Identifier"
    
  ),
  
  drop = FALSE
  
]


Table_S14B <- Table_S14B[
  
  order(
    
    Table_S14B$Phenotype,
    
    Table_S14B$`Lipid main ID`
    
  ),
  
  ,
  
  drop = FALSE
  
]


# ============================================================
# 26. PRINT FINAL RESULT
# ============================================================

cat(
  "\n============================================\n",
  "FINAL TABLE S14A\n",
  "============================================\n",
  sep = ""
)


print(
  Table_S14A,
  row.names = FALSE
)


# ============================================================
# 27. HARD VALIDATION OF KNOWN VALUES
# ============================================================

expected_S14A <- data.frame(
  
  Phenotype = c(
    "Steatosis",
    "Inflammation",
    "Fibrosis"
  ),
  
  Discovery_significant = c(
    826,
    903,
    707
  ),
  
  External_significant = c(
    21,
    0,
    33
  ),
  
  Measured_externally = c(
    108,
    110,
    104
  ),
  
  Same_direction = c(
    86,
    67,
    86
  ),
  
  Opposite_direction = c(
    22,
    43,
    18
  )
  
)


for (
  i in seq_len(
    nrow(expected_S14A)
  )
) {
  
  z <- Table_S14A[
    Table_S14A$Phenotype ==
      expected_S14A$Phenotype[i],
    ,
    drop = FALSE
  ]
  
  
  if (
    nrow(z) != 1
  ) {
    
    stop(
      "\nMissing S14A row for ",
      expected_S14A$Phenotype[i]
    )
    
  }
  
  
  if (
    z$Discovery_significant !=
    expected_S14A$Discovery_significant[i]
  ) {
    
    stop(
      "\nDiscovery count mismatch for ",
      expected_S14A$Phenotype[i]
    )
    
  }
  
  
  if (
    z$External_significant !=
    expected_S14A$External_significant[i]
  ) {
    
    stop(
      "\nExternal significant count mismatch for ",
      expected_S14A$Phenotype[i]
    )
    
  }
  
  
  if (
    z$Measured_externally !=
    expected_S14A$Measured_externally[i]
  ) {
    
    stop(
      "\nMeasured overlap mismatch for ",
      expected_S14A$Phenotype[i]
    )
    
  }
  
  
  if (
    z$Same_direction !=
    expected_S14A$Same_direction[i]
  ) {
    
    stop(
      "\nSame-direction count mismatch for ",
      expected_S14A$Phenotype[i]
    )
    
  }
  
  
  if (
    z$Opposite_direction !=
    expected_S14A$Opposite_direction[i]
  ) {
    
    stop(
      "\nOpposite-direction count mismatch for ",
      expected_S14A$Phenotype[i]
    )
    
  }
  
}


cat(
  "\n============================================\n",
  "LIPIDOMICS VALIDATION PASSED\n",
  "============================================\n",
  sep = ""
)



# ============================================================
# 24. TRANSCRIPTOMICS
# ============================================================
#
# EXACT RNA DISCOVERY OBJECTS FROM YOUR ORIGINAL ANALYSIS:
#
#   resLFC_steatosis
#   resLFC_inflammation
#   resLFC_fibrosis
#
# These are the DESeq2 + apeglm RNA-seq results.
#
# EXACT EXTERNAL VALIDATION OBJECTS:
#
#   val_GSE162694
#   val_GSE135251_unadj
#   val_GSE135251_sexadj
#
# ============================================================


# ============================================================
# 24.1 CHECK EXACT RNA DISCOVERY OBJECTS
# ============================================================

required_rna_discovery <- c(
  "resLFC_steatosis",
  "resLFC_inflammation",
  "resLFC_fibrosis"
)


missing_rna_discovery <- required_rna_discovery[
  !vapply(
    required_rna_discovery,
    exists,
    logical(1),
    inherits = TRUE
  )
]


if (
  length(missing_rna_discovery) > 0
) {
  
  stop(
    "\nMissing RNA-seq discovery objects:\n\n",
    paste(
      missing_rna_discovery,
      collapse = "\n"
    ),
    "\n\nThese are the exact object names from your RNA-seq analysis."
  )
  
}


# ============================================================
# 24.2 CHECK EXACT EXTERNAL VALIDATION OBJECTS
# ============================================================

required_rna_validation <- c(
  "val_GSE162694",
  "val_GSE135251_unadj",
  "val_GSE135251_sexadj"
)


missing_rna_validation <- required_rna_validation[
  !vapply(
    required_rna_validation,
    exists,
    logical(1),
    inherits = TRUE
  )
]


if (
  length(missing_rna_validation) > 0
) {
  
  stop(
    "\nMissing RNA external-validation objects:\n\n",
    paste(
      missing_rna_validation,
      collapse = "\n"
    ),
    "\n\nRun the corresponding external transcriptomic validation chunks first."
  )
  
}


# ============================================================
# 24.3 ASSIGN DISCOVERY OBJECTS
# ============================================================

rna_disc_steatosis <-
  resLFC_steatosis

rna_disc_inflammation <-
  resLFC_inflammation

rna_disc_fibrosis <-
  resLFC_fibrosis


# ============================================================
# 24.4 VERIFY RNA OBJECT STRUCTURE
# ============================================================

cat(
  "\n============================================\n",
  "RNA DISCOVERY OBJECT CHECK\n",
  "============================================\n",
  sep = ""
)


cat(
  "\nSteatosis:\n"
)

print(
  class(
    rna_disc_steatosis
  )
)

print(
  dim(
    as.data.frame(
      rna_disc_steatosis
    )
  )
)


cat(
  "\nInflammation:\n"
)

print(
  class(
    rna_disc_inflammation
  )
)

print(
  dim(
    as.data.frame(
      rna_disc_inflammation
    )
  )
)


cat(
  "\nFibrosis:\n"
)

print(
  class(
    rna_disc_fibrosis
  )
)

print(
  dim(
    as.data.frame(
      rna_disc_fibrosis
    )
  )
)







# ============================================================
# 25. RNA VALIDATION OBJECT CHECK
# ============================================================
# ============================================================
# 25. RNA RESULT PREPARATION
# ============================================================

prepare_rna_table <- function(
    x,
    cohort,
    phenotype,
    fdr = 0.05
) {
  
  x <- as.data.frame(
    x
  )
  
  
  # ----------------------------------------------------------
  # Required DESeq2 columns
  # ----------------------------------------------------------
  
  if (
    !all(
      c(
        "log2FoldChange",
        "padj"
      ) %in%
      names(x)
    )
  ) {
    
    stop(
      "\nRNA result for ",
      cohort,
      " / ",
      phenotype,
      " does not contain log2FoldChange and padj.\n\n",
      "Available columns:\n",
      paste(
        names(x),
        collapse = ", "
      )
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Gene identifiers = rownames
  # ----------------------------------------------------------
  
  gene <- rownames(
    x
  )
  
  
  if (
    is.null(gene)
  ) {
    
    stop(
      "\nNo gene identifiers found in row names for ",
      cohort,
      " / ",
      phenotype
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Significant genes
  # ----------------------------------------------------------
  
  sig <- !is.na(
    x$padj
  ) &
    x$padj < fdr
  
  
  out <- data.frame(
    
    Gene =
      as.character(
        gene[sig]
      ),
    
    logFC =
      as.numeric(
        x$log2FoldChange[sig]
      ),
    
    FDR =
      as.numeric(
        x$padj[sig]
      ),
    
    Cohort =
      rep(
        cohort,
        sum(sig)
      ),
    
    Phenotype =
      rep(
        phenotype,
        sum(sig)
      ),
    
    stringsAsFactors =
      FALSE
    
  )
  
  
  out$Direction <-
    ifelse(
      
      out$logFC > 0,
      
      "Up",
      
      ifelse(
        
        out$logFC < 0,
        
        "Down",
        
        "Zero"
        
      )
      
    )
  
  
  out
  
}


# ============================================================
# 25.1 DISCOVERY RNA
# ============================================================

disc_rna <- rbind(
  
  prepare_rna_table(
    rna_disc_steatosis,
    "Discovery",
    "Steatosis",
    FDR_RNA
  ),
  
  prepare_rna_table(
    rna_disc_inflammation,
    "Discovery",
    "Inflammation",
    FDR_RNA
  ),
  
  prepare_rna_table(
    rna_disc_fibrosis,
    "Discovery",
    "Fibrosis",
    FDR_RNA
  )
  
)


rownames(
  disc_rna
) <- NULL


cat(
  "\n============================================\n",
  "DISCOVERY RNA FDR < 0.05\n",
  "============================================\n",
  sep = ""
)


for (
  ph in c(
    "Steatosis",
    "Inflammation",
    "Fibrosis"
  )
) {
  
  cat(
    ph,
    ":",
    sum(
      disc_rna$Phenotype == ph
    ),
    "\n"
  )
  
}


# ============================================================
# 26. EXTRACT RNA VALIDATION RESULTS
# ============================================================

val162694_fib <-
  prepare_rna_table(
    val_GSE162694$fib_res,
    "GSE162694",
    "Fibrosis",
    FDR_RNA
  )


val162694_nas <-
  prepare_rna_table(
    val_GSE162694$nas_res,
    "GSE162694",
    "NAS",
    FDR_RNA
  )


val135251_fib <-
  prepare_rna_table(
    val_GSE135251_unadj$fib_res,
    "GSE135251",
    "Fibrosis",
    FDR_RNA
  )


val135251_nas <-
  prepare_rna_table(
    val_GSE135251_unadj$nas_res,
    "GSE135251",
    "NAS",
    FDR_RNA
  )


val135251_sex_fib <-
  prepare_rna_table(
    val_GSE135251_sexadj$fib_res,
    "GSE135251 (Sex-adjusted)",
    "Fibrosis",
    FDR_RNA
  )


val135251_sex_nas <-
  prepare_rna_table(
    val_GSE135251_sexadj$nas_res,
    "GSE135251 (Sex-adjusted)",
    "NAS",
    FDR_RNA
  )




# ============================================================
# 27. RNA DISCOVERY vs EXTERNAL COMPARISON
# ============================================================

compare_rna <- function(
    discovery,
    validation,
    discovery_trait,
    validation_trait,
    cohort
) {
  
  d <- discovery[
    discovery$Phenotype ==
      discovery_trait,
    ,
    drop = FALSE
  ]
  
  
  v <- validation[
    validation$Phenotype ==
      validation_trait,
    ,
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # Match genes
  # ----------------------------------------------------------
  
  overlap <- intersect(
    d$Gene,
    v$Gene
  )
  
  
  if (
    length(overlap) == 0
  ) {
    
    return(
      list(
        
        feature =
          data.frame(),
        
        summary =
          data.frame(
            
            Cohort =
              cohort,
            
            Discovery_trait =
              discovery_trait,
            
            Validation_trait =
              validation_trait,
            
            Shared_significant =
              0,
            
            Spearman_rho =
              NA_real_,
            
            Spearman_P =
              NA_real_,
            
            Percent_concordant =
              NA_real_,
            
            Binomial_P =
              NA_real_,
            
            stringsAsFactors =
              FALSE
            
          )
        
      )
    )
    
  }
  
  
  d2 <- d[
    match(
      overlap,
      d$Gene
    ),
    ,
    drop = FALSE
  ]
  
  
  v2 <- v[
    match(
      overlap,
      v$Gene
    ),
    ,
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # Remove missing effect sizes
  # ----------------------------------------------------------
  
  keep <- complete.cases(
    d2$logFC,
    v2$logFC
  )
  
  
  d2 <- d2[
    keep,
    ,
    drop = FALSE
  ]
  
  
  v2 <- v2[
    keep,
    ,
    drop = FALSE
  ]
  
  
  if (
    nrow(d2) == 0
  ) {
    
    return(
      list(
        
        feature =
          data.frame(),
        
        summary =
          data.frame(
            
            Cohort =
              cohort,
            
            Discovery_trait =
              discovery_trait,
            
            Validation_trait =
              validation_trait,
            
            Shared_significant =
              0,
            
            Spearman_rho =
              NA_real_,
            
            Spearman_P =
              NA_real_,
            
            Percent_concordant =
              NA_real_,
            
            Binomial_P =
              NA_real_,
            
            stringsAsFactors =
              FALSE
            
          )
        
      )
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Direction
  # ----------------------------------------------------------
  
  same_direction <-
    sign(
      d2$logFC
    ) ==
    sign(
      v2$logFC
    )
  
  
  # ----------------------------------------------------------
  # Spearman correlation
  # ----------------------------------------------------------
  
  if (
    nrow(d2) >= 2
  ) {
    
    rho_test <-
      suppressWarnings(
        cor.test(
          d2$logFC,
          v2$logFC,
          method = "spearman",
          exact = FALSE
        )
      )
    
    rho <-
      unname(
        rho_test$estimate
      )
    
    rho_p <-
      rho_test$p.value
    
  } else {
    
    rho <- NA_real_
    rho_p <- NA_real_
    
  }
  
  
  # ----------------------------------------------------------
  # Binomial test
  # ----------------------------------------------------------
  
  if (
    length(same_direction) > 0
  ) {
    
    binom_test <-
      binom.test(
        sum(
          same_direction
        ),
        length(
          same_direction
        ),
        p = 0.5,
        alternative = "greater"
      )
    
    binom_p <-
      binom_test$p.value
    
  } else {
    
    binom_p <-
      NA_real_
    
  }
  
  
  # ----------------------------------------------------------
  # Feature table
  # ----------------------------------------------------------
  
  feature <- data.frame(
    
    Cohort =
      cohort,
    
    Discovery_trait =
      discovery_trait,
    
    Validation_trait =
      validation_trait,
    
    Gene =
      d2$Gene,
    
    Discovery_logFC =
      d2$logFC,
    
    Discovery_FDR =
      d2$FDR,
    
    Validation_logFC =
      v2$logFC,
    
    Validation_FDR =
      v2$FDR,
    
    Discovery_Direction =
      d2$Direction,
    
    Validation_Direction =
      v2$Direction,
    
    Direction_concordant =
      same_direction,
    
    stringsAsFactors =
      FALSE
    
  )
  
  
  # ----------------------------------------------------------
  # Summary
  # ----------------------------------------------------------
  
  summary <- data.frame(
    
    Cohort =
      cohort,
    
    Discovery_trait =
      discovery_trait,
    
    Validation_trait =
      validation_trait,
    
    Shared_significant =
      nrow(feature),
    
    Spearman_rho =
      rho,
    
    Spearman_P =
      rho_p,
    
    Percent_concordant =
      100 *
      mean(
        same_direction
      ),
    
    Binomial_P =
      binom_p,
    
    stringsAsFactors =
      FALSE
    
  )
  
  
  list(
    
    feature =
      feature,
    
    summary =
      summary
    
  )
  
}


# ============================================================
# 28. RUN ALL RNA VALIDATIONS
# ============================================================


# ------------------------------------------------------------
# GSE162694
# ------------------------------------------------------------

rna_fib_162694 <-
  compare_rna(
    disc_rna,
    val162694_fib,
    "Fibrosis",
    "Fibrosis",
    "GSE162694"
  )


rna_steat_162694 <-
  compare_rna(
    disc_rna,
    val162694_nas,
    "Steatosis",
    "NAS",
    "GSE162694"
  )


rna_infl_162694 <-
  compare_rna(
    disc_rna,
    val162694_nas,
    "Inflammation",
    "NAS",
    "GSE162694"
  )


# ------------------------------------------------------------
# GSE135251 unadjusted
# ------------------------------------------------------------

rna_fib_135251 <-
  compare_rna(
    disc_rna,
    val135251_fib,
    "Fibrosis",
    "Fibrosis",
    "GSE135251"
  )


rna_steat_135251 <-
  compare_rna(
    disc_rna,
    val135251_nas,
    "Steatosis",
    "NAS",
    "GSE135251"
  )


rna_infl_135251 <-
  compare_rna(
    disc_rna,
    val135251_nas,
    "Inflammation",
    "NAS",
    "GSE135251"
  )


# ------------------------------------------------------------
# GSE135251 sex-adjusted
# ------------------------------------------------------------

rna_fib_135251_sex <-
  compare_rna(
    disc_rna,
    val135251_sex_fib,
    "Fibrosis",
    "Fibrosis",
    "GSE135251 (Sex-adjusted)"
  )


rna_steat_135251_sex <-
  compare_rna(
    disc_rna,
    val135251_sex_nas,
    "Steatosis",
    "NAS",
    "GSE135251 (Sex-adjusted)"
  )


rna_infl_135251_sex <-
  compare_rna(
    disc_rna,
    val135251_sex_nas,
    "Inflammation",
    "NAS",
    "GSE135251 (Sex-adjusted)"
  )





# ============================================================
# 29. TABLE S11
# ============================================================

Table_S11 <- rbind(
  
  rna_fib_162694$summary,
  
  rna_steat_162694$summary,
  
  rna_infl_162694$summary,
  
  rna_fib_135251$summary,
  
  rna_steat_135251$summary,
  
  rna_infl_135251$summary,
  
  rna_fib_135251_sex$summary,
  
  rna_steat_135251_sex$summary,
  
  rna_infl_135251_sex$summary
  
)


Table_S11$Spearman_rho <-
  round(
    Table_S11$Spearman_rho,
    3
  )


Table_S11$Spearman_P <-
  signif(
    Table_S11$Spearman_P,
    4
  )


Table_S11$Percent_concordant <-
  round(
    Table_S11$Percent_concordant,
    1
  )


Table_S11$Binomial_P <-
  signif(
    Table_S11$Binomial_P,
    4
  )


# ============================================================
# TABLE S12A
# ============================================================

Table_S12A <- rbind(
  
  rna_fib_162694$feature,
  
  rna_steat_162694$feature,
  
  rna_infl_162694$feature,
  
  rna_fib_135251$feature,
  
  rna_steat_135251$feature,
  
  rna_infl_135251$feature,
  
  rna_fib_135251_sex$feature,
  
  rna_steat_135251_sex$feature,
  
  rna_infl_135251_sex$feature
  
)


rownames(
  Table_S12A
) <- NULL


# ============================================================
# PRINT S11
# ============================================================

cat(
  "\n============================================\n",
  "TABLE S11\n",
  "============================================\n",
  sep = ""
)

print(
  Table_S11,
  row.names = FALSE
)







# ============================================================
# 29. PREPARE EXTERNAL RNA VALIDATION
# ============================================================
#
# The validation objects already contain:
#
#   val_GSE162694$fib_res
#   val_GSE162694$nas_res
#
#   val_GSE135251_unadj$fib_res
#   val_GSE135251_unadj$nas_res
#
#   val_GSE135251_sexadj$fib_res
#   val_GSE135251_sexadj$nas_res
#
# Therefore no get_fibrosis_result() or get_nas_result()
# helper is needed.
# ============================================================


val162694_fib <- prepare_rna_table(
  
  val_GSE162694$fib_res,
  
  "GSE162694",
  
  "Fibrosis",
  
  FDR_RNA
  
)


val162694_nas <- prepare_rna_table(
  
  val_GSE162694$nas_res,
  
  "GSE162694",
  
  "NAS",
  
  FDR_RNA
  
)


val135251_fib <- prepare_rna_table(
  
  val_GSE135251_unadj$fib_res,
  
  "GSE135251",
  
  "Fibrosis",
  
  FDR_RNA
  
)


val135251_nas <- prepare_rna_table(
  
  val_GSE135251_unadj$nas_res,
  
  "GSE135251",
  
  "NAS",
  
  FDR_RNA
  
)


val135251_sex_fib <- prepare_rna_table(
  
  val_GSE135251_sexadj$fib_res,
  
  "GSE135251 (Sex-adjusted)",
  
  "Fibrosis",
  
  FDR_RNA
  
)


val135251_sex_nas <- prepare_rna_table(
  
  val_GSE135251_sexadj$nas_res,
  
  "GSE135251 (Sex-adjusted)",
  
  "NAS",
  
  FDR_RNA
  
)


cat(
  "\nExternal RNA validation tables prepared successfully.\n"
)



# ============================================================
# 30. RNA COMPARISON
# ============================================================

compare_rna <- function(
    
  discovery,
  
  validation,
  
  discovery_trait,
  
  validation_trait,
  
  cohort
  
) {
  
  
  d <- discovery[
    discovery$Phenotype ==
      discovery_trait,
    ,
    drop = FALSE
  ]
  
  
  v <- validation[
    validation$Phenotype ==
      validation_trait,
    ,
    drop = FALSE
  ]
  
  
  d <- d[
    ,
    c(
      "Gene",
      "logFC",
      "FDR",
      "Direction"
    ),
    drop = FALSE
  ]
  
  
  names(d) <- c(
    "Gene",
    "Discovery_logFC",
    "Discovery_FDR",
    "Discovery_Direction"
  )
  
  
  v <- v[
    ,
    c(
      "Gene",
      "logFC",
      "FDR",
      "Direction"
    ),
    drop = FALSE
  ]
  
  
  names(v) <- c(
    "Gene",
    "Validation_logFC",
    "Validation_FDR",
    "Validation_Direction"
  )
  
  
  merged <- merge(
    
    d,
    
    v,
    
    by = "Gene",
    
    all = FALSE
    
  )
  
  
  if (
    nrow(merged) == 0
  ) {
    
    return(
      list(
        
        feature =
          merged,
        
        summary =
          data.frame(
            
            Cohort =
              cohort,
            
            Discovery_trait =
              discovery_trait,
            
            Validation_trait =
              validation_trait,
            
            Shared_significant =
              0,
            
            Spearman_rho =
              NA_real_,
            
            Spearman_P =
              NA_real_,
            
            Percent_concordant =
              NA_real_,
            
            Binomial_P =
              NA_real_,
            
            stringsAsFactors = FALSE
            
          )
        
      )
    )
    
  }
  
  
  same_direction <-
    merged$Discovery_Direction ==
    merged$Validation_Direction
  
  
  rho_test <- suppressWarnings(
    
    cor.test(
      
      merged$Discovery_logFC,
      
      merged$Validation_logFC,
      
      method = "spearman",
      
      exact = FALSE
      
    )
    
  )
  
  
  binom_test <- binom.test(
    
    sum(
      same_direction
    ),
    
    length(
      same_direction
    ),
    
    p = 0.5,
    
    alternative = "greater"
    
  )
  
  
  merged$Cohort <-
    cohort
  
  
  merged$Discovery_trait <-
    discovery_trait
  
  
  merged$Validation_trait <-
    validation_trait
  
  
  merged$Direction_concordant <-
    same_direction
  
  
  summary <- data.frame(
    
    Cohort =
      cohort,
    
    Discovery_trait =
      discovery_trait,
    
    Validation_trait =
      validation_trait,
    
    Shared_significant =
      nrow(merged),
    
    Spearman_rho =
      unname(
        rho_test$estimate
      ),
    
    Spearman_P =
      rho_test$p.value,
    
    Percent_concordant =
      100 *
      mean(
        same_direction
      ),
    
    Binomial_P =
      binom_test$p.value,
    
    stringsAsFactors = FALSE
    
  )
  
  
  list(
    
    feature =
      merged,
    
    summary =
      summary
    
  )
  
}

# ============================================================
# 31. RUN RNA COMPARISONS
# ============================================================

rna_fib_162694 <-
  compare_rna(
    disc_rna,
    val162694_fib,
    "Fibrosis",
    "Fibrosis",
    "GSE162694"
  )

rna_steat_162694 <-
  compare_rna(
    disc_rna,
    val162694_nas,
    "Steatosis",
    "NAS",
    "GSE162694"
  )

rna_infl_162694 <-
  compare_rna(
    disc_rna,
    val162694_nas,
    "Inflammation",
    "NAS",
    "GSE162694"
  )


rna_fib_135251 <-
  compare_rna(
    disc_rna,
    val135251_fib,
    "Fibrosis",
    "Fibrosis",
    "GSE135251"
  )

rna_steat_135251 <-
  compare_rna(
    disc_rna,
    val135251_nas,
    "Steatosis",
    "NAS",
    "GSE135251"
  )

rna_infl_135251 <-
  compare_rna(
    disc_rna,
    val135251_nas,
    "Inflammation",
    "NAS",
    "GSE135251"
  )


rna_fib_135251_sex <-
  compare_rna(
    disc_rna,
    val135251_sex_fib,
    "Fibrosis",
    "Fibrosis",
    "GSE135251 (Sex-adjusted)"
  )

rna_steat_135251_sex <-
  compare_rna(
    disc_rna,
    val135251_sex_nas,
    "Steatosis",
    "NAS",
    "GSE135251 (Sex-adjusted)"
  )

rna_infl_135251_sex <-
  compare_rna(
    disc_rna,
    val135251_sex_nas,
    "Inflammation",
    "NAS",
    "GSE135251 (Sex-adjusted)"
  )



# ============================================================
# 32. TABLE S11
# ============================================================

Table_S11 <- rbind(
  
  rna_fib_162694$summary,
  rna_steat_162694$summary,
  rna_infl_162694$summary,
  
  rna_fib_135251$summary,
  rna_steat_135251$summary,
  rna_infl_135251$summary,
  
  rna_fib_135251_sex$summary,
  rna_steat_135251_sex$summary,
  rna_infl_135251_sex$summary
  
)


Table_S11$Spearman_rho <-
  round(
    Table_S11$Spearman_rho,
    3
  )

Table_S11$Spearman_P <-
  signif(
    Table_S11$Spearman_P,
    4
  )

Table_S11$Percent_concordant <-
  round(
    Table_S11$Percent_concordant,
    1
  )

Table_S11$Binomial_P <-
  signif(
    Table_S11$Binomial_P,
    4
  )


cat(
  "\n============================================\n",
  "FINAL TABLE S11\n",
  "============================================\n",
  sep = ""
)

print(
  Table_S11,
  row.names = FALSE
)



# ============================================================
# 33. TABLE S12A
# ============================================================

Table_S12A <- rbind(
  
  rna_fib_162694$feature,
  rna_steat_162694$feature,
  rna_infl_162694$feature,
  
  rna_fib_135251$feature,
  rna_steat_135251$feature,
  rna_infl_135251$feature,
  
  rna_fib_135251_sex$feature,
  rna_steat_135251_sex$feature,
  rna_infl_135251_sex$feature
  
)

rownames(
  Table_S12A
) <- NULL


cat(
  "\nS12A rows:",
  nrow(Table_S12A),
  "\n"
)


# ============================================================
# 34. TABLE S12B
# ============================================================
#
# Use KEGG objects that already exist in the validation objects.
# No pathway results are invented.
#
# ============================================================
# ============================================================
# REPLACE extract_kegg() WITH THIS VERSION
# ============================================================

extract_kegg <- function(val_obj, cohort, phenotype) {
  
  out <- list()
  
  # ----------------------------
  # Select the correct KEGG objects
  # ----------------------------
  
  if (phenotype == "Fibrosis") {
    
    kegg_up <- val_obj$kegg_fib_up
    kegg_down <- val_obj$kegg_fib_down
    
  } else if (phenotype == "NAS") {
    
    kegg_up <- val_obj$kegg_nas_up
    kegg_down <- val_obj$kegg_nas_down
    
  } else {
    
    stop("Unknown phenotype: ", phenotype)
    
  }
  
  
  # ----------------------------
  # Convert KEGG object to data.frame
  # ----------------------------
  
  kegg_to_df <- function(kegg_obj, direction) {
    
    if (is.null(kegg_obj)) {
      return(NULL)
    }
    
    if (!methods::is(kegg_obj, "enrichResult")) {
      return(NULL)
    }
    
    if (nrow(kegg_obj@result) == 0) {
      return(NULL)
    }
    
    z <- as.data.frame(kegg_obj@result)
    
    z$Cohort <- cohort
    z$Phenotype <- phenotype
    z$Direction <- direction
    
    z
    
  }
  
  
  # ----------------------------
  # Upregulated pathways
  # ----------------------------
  
  up <- kegg_to_df(
    kegg_up,
    "Up"
  )
  
  
  # ----------------------------
  # Downregulated pathways
  # ----------------------------
  
  down <- kegg_to_df(
    kegg_down,
    "Down"
  )
  
  
  # ----------------------------
  # Combine
  # ----------------------------
  
  pieces <- Filter(
    Negate(is.null),
    list(up, down)
  )
  
  
  if (length(pieces) == 0) {
    return(NULL)
  }
  
  
  do.call(
    rbind,
    pieces
  )
}

Table_S12B <- do.call(
  
  rbind,
  
  Filter(
    Negate(is.null),
    
    list(
      
      extract_kegg(
        val_GSE162694,
        "GSE162694",
        "Fibrosis"
      ),
      
      extract_kegg(
        val_GSE162694,
        "GSE162694",
        "NAS"
      ),
      
      extract_kegg(
        val_GSE135251_unadj,
        "GSE135251",
        "Fibrosis"
      ),
      
      extract_kegg(
        val_GSE135251_unadj,
        "GSE135251",
        "NAS"
      ),
      
      extract_kegg(
        val_GSE135251_sexadj,
        "GSE135251 (Sex-adjusted)",
        "Fibrosis"
      ),
      
      extract_kegg(
        val_GSE135251_sexadj,
        "GSE135251 (Sex-adjusted)",
        "NAS"
      )
      
    )
    
  )
  
)


if (
  is.null(Table_S12B)
) {
  
  Table_S12B <-
    data.frame(
      
      Note =
        "No KEGG result objects were available in the supplied transcriptomic validation objects.",
      
      stringsAsFactors = FALSE
      
    )
  
}


# ============================================================
# 35. FINAL TABLE CHECKS
# ============================================================

cat(
  "\n============================================\n",
  "FINAL TABLE CHECKS\n",
  "============================================\n",
  sep = ""
)


cat(
  "\nS11 rows:",
  nrow(Table_S11),
  "\n"
)


cat(
  "S12A rows:",
  nrow(Table_S12A),
  "\n"
)


cat(
  "S12B rows:",
  nrow(Table_S12B),
  "\n"
)


cat(
  "S13A rows:",
  nrow(Table_S13A),
  "\n"
)


cat(
  "S13B rows:",
  nrow(Table_S13B),
  "\n"
)


cat(
  "S14A rows:",
  nrow(Table_S14A),
  "\n"
)


cat(
  "S14B rows:",
  nrow(Table_S14B),
  "\n"
)


# ============================================================
# 36. CREATE EXCEL WORKBOOK
# ============================================================

wb <- createWorkbook()


# ------------------------------------------------------------
# Sheets
# ------------------------------------------------------------

addWorksheet(
  wb,
  "Table_S11"
)


addWorksheet(
  wb,
  "Table_S12A_Genes"
)


addWorksheet(
  wb,
  "Table_S12B_KEGG"
)


addWorksheet(
  wb,
  "Table_S13A_Metabolomics"
)


addWorksheet(
  wb,
  "Table_S13B_Metabolites"
)


addWorksheet(
  wb,
  "Table_S14A_Lipidomics"
)


addWorksheet(
  wb,
  "Table_S14B_Lipids"
)


# ============================================================
# 37. WRITE TABLES
# ============================================================

writeData(
  wb,
  "Table_S11",
  Table_S11,
  withFilter = TRUE
)


writeData(
  wb,
  "Table_S12A_Genes",
  Table_S12A,
  withFilter = TRUE
)


writeData(
  wb,
  "Table_S12B_KEGG",
  Table_S12B,
  withFilter = TRUE
)


writeData(
  wb,
  "Table_S13A_Metabolomics",
  Table_S13A,
  withFilter = TRUE
)


writeData(
  wb,
  "Table_S13B_Metabolites",
  Table_S13B,
  withFilter = TRUE
)


writeData(
  wb,
  "Table_S14A_Lipidomics",
  Table_S14A,
  withFilter = TRUE
)


writeData(
  wb,
  "Table_S14B_Lipids",
  Table_S14B,
  withFilter = TRUE
)


# ============================================================
# 38. EXCEL STYLING
# ============================================================
#
# Correct openxlsx argument:
#
#   fgFill
#
# NOT:
#
#   Fill
#
# ============================================================

header_style <- createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "bottom",
  halign = "center",
  valign = "center"
)


body_style <- createStyle(
  
  wrapText =
    TRUE,
  
  valign =
    "top"
  
)


# ============================================================
# 39. STYLE WORKSHEETS
# ============================================================
#
# IMPORTANT:
#
#   names(wb$worksheets)
#
# NOT:
#
#   names(wb)
#
# ============================================================

sheet_names <-
  names(
    wb$worksheets
  )


for (
  sheet in sheet_names
) {
  
  df <- readWorkbook(
    wb,
    sheet = sheet,
    colNames = TRUE
  )
  
  
  if (
    is.null(df) ||
    ncol(df) == 0
  ) {
    
    next
    
  }
  
  
  # ----------------------------------------------------------
  # Header
  # ----------------------------------------------------------
  
  addStyle(
    
    wb,
    
    sheet = sheet,
    
    style =
      header_style,
    
    rows = 1,
    
    cols =
      seq_len(
        ncol(df)
      ),
    
    gridExpand =
      TRUE,
    
    stack =
      TRUE
    
  )
  
  
  # ----------------------------------------------------------
  # Body
  # ----------------------------------------------------------
  
  if (
    nrow(df) > 0
  ) {
    
    addStyle(
      
      wb,
      
      sheet = sheet,
      
      style =
        body_style,
      
      rows =
        2:(nrow(df) + 1),
      
      cols =
        seq_len(
          ncol(df)
        ),
      
      gridExpand =
        TRUE,
      
      stack =
        TRUE
      
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Width
  # ----------------------------------------------------------
  
  setColWidths(
    
    wb,
    
    sheet = sheet,
    
    cols =
      seq_len(
        ncol(df)
      ),
    
    widths =
      "auto"
    
  )
  
  
  # ----------------------------------------------------------
  # Freeze top row
  # ----------------------------------------------------------
  
  freezePane(
    
    wb,
    
    sheet = sheet,
    
    firstRow =
      TRUE
    
  )
  
  
  # ----------------------------------------------------------
  # Filter
  # ----------------------------------------------------------
  
  addFilter(
    
    wb,
    
    sheet = sheet,
    
    rows = 1,
    
    cols =
      seq_len(
        ncol(df)
      )
    
  )
  
}


# ============================================================
# 40. NUMBER FORMATTING
# ============================================================

number_style <- createStyle(
  numFmt = "0.0000"
)


for (
  sheet in sheet_names
) {
  
  df <- readWorkbook(
    wb,
    sheet = sheet,
    colNames = TRUE
  )
  
  
  if (
    is.null(df) ||
    nrow(df) == 0
  ) {
    
    next
    
  }
  
  
  numeric_cols <- which(
    vapply(
      df,
      is.numeric,
      logical(1)
    )
  )
  
  
  if (
    length(numeric_cols) > 0
  ) {
    
    addStyle(
      
      wb,
      
      sheet = sheet,
      
      style =
        number_style,
      
      rows =
        2:(nrow(df) + 1),
      
      cols =
        numeric_cols,
      
      gridExpand =
        TRUE,
      
      stack =
        TRUE
      
    )
    
  }
  
}


# ============================================================
# 41. SAVE FINAL WORKBOOK
# ============================================================

saveWorkbook(
  
  wb,
  
  OUTPUT_FILE,
  
  overwrite =
    TRUE
  
)


# ============================================================
# 42. VERIFY FILE EXISTS
# ============================================================

if (
  !file.exists(
    OUTPUT_FILE
  )
) {
  
  stop(
    "\nWorkbook was not created."
  )
  
}


# ============================================================
# 43. FINAL SUMMARY
# ============================================================

cat(
  "\n\n============================================\n",
  "FINAL EXTERNAL VALIDATION WORKBOOK CREATED\n",
  "============================================\n",
  sep = ""
)


cat(
  "\nFile:\n",
  normalizePath(
    OUTPUT_FILE,
    mustWork = FALSE
  ),
  "\n"
)


cat(
  "\nEPoS:\n",
  "Metabolomics = 36 features x 627 samples\n",
  "Lipidomics    = 176 features x 627 samples\n",
  "Steatosis     = 212 vs 388\n",
  "Inflammation  = 307 vs 247\n",
  "Fibrosis      = 157 vs 293\n"
)


cat(
  "\nExternal serum metabolomics:\n",
  "Steatosis     = 0 significant\n",
  "Inflammation  = 6 significant\n",
  "Fibrosis      = 0 significant\n"
)


cat(
  "\nDiscovery serum metabolomics:\n",
  "Steatosis     = 12 significant\n",
  "Inflammation  = 4 significant\n",
  "Fibrosis      = 0 significant\n"
)


cat(
  "\nExternal serum lipidomics:\n",
  "Steatosis     = 21 significant\n",
  "Inflammation  = 0 significant\n",
  "Fibrosis      = 33 significant\n"
)


cat(
  "\nDiscovery serum lipidomics:\n",
  "Steatosis     = 826 significant\n",
  "Inflammation  = 903 significant\n",
  "Fibrosis      = 707 significant\n"
)


cat(
  "\nLipid directional overlap:\n",
  "Steatosis     = 108 measured / 86 same / 22 opposite\n",
  "Inflammation  = 110 measured / 67 same / 43 opposite\n",
  "Fibrosis      = 104 measured / 86 same / 18 opposite\n"
)


cat(
  "\nWorkbook sheets:\n"
)

print(
  names(
    wb$worksheets
  )
)


cat(
  "\n============================================\n",
  "DONE\n",
  "============================================\n",
  sep = ""
)
