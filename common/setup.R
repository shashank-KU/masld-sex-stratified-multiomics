# common/setup.R -- shared setup for every analysis script in this repository.
# Sourced by the first chunk of each <Folder>/<Folder>.Rmd. Defines paths, thresholds and helpers.

# ---- Repository root (folder that contains common/) ----------------------
if (!exists("REPO_ROOT")) {
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "common", "setup.R")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, "common", "setup.R"))) stop("Run from inside the repository (common/setup.R not found).")
  REPO_ROOT <- d
}

# ---- Input data (not in the repository, see data/README.md) -------------
# Order: environment variable MASLD_DATA_DIR, then <repo>/data. The folder must contain Analyisis/ etc.
.cands <- c(Sys.getenv("MASLD_DATA_DIR", unset = NA), file.path(REPO_ROOT, "data"))
.cands <- .cands[!is.na(.cands) & nzchar(.cands)]
.hit <- .cands[dir.exists(.cands)]
if (!length(.hit)) stop("Input folder not found. Put the data in <repo>/data or set\n  Sys.setenv(MASLD_DATA_DIR = \"/path/to/folder/with/Analyisis\")")
project_dir <- normalizePath(.hit[1])

# ---- Results (git-ignored): tables, figures, cached R workspaces --------
.res <- Sys.getenv("MASLD_RESULTS_DIR", unset = file.path(REPO_ROOT, "results"))
dir.create(file.path(.res, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(.res, "cache"),  showWarnings = FALSE, recursive = TRUE)
results_dir <- normalizePath(.res)
OUT_DIR <- results_dir
knitr::opts_knit$set(root.dir = results_dir)   # relative output paths land in results/
setwd(results_dir)
message("Repository: ", REPO_ROOT, "\nInput folder: ", project_dir, "\nResults folder: ", results_dir)

# ---- Supplementary Tables (published) ---------------------------------------
# Some steps read the published tables (Table S3 curated lipid names, Table S11d-e AUCs, Table S12).
# Order: MASLD_SUPP_XLSX, any gupta_MASLD_SupplementaryTables*.xlsx in data/ or the repo root,
# otherwise a workbook is rebuilt from the Table_S*.csv files in the <Figure>/Source_data/ folders.
build_supp_from_csv <- function(out) {
  f <- list.files(REPO_ROOT, pattern = "^Table[_ ]S[0-9]+[ab]?\\.csv$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  f <- f[!grepl("/results/", f)]
  f <- f[!duplicated(basename(f))]
  if (!length(f)) return(NULL)
  wb <- openxlsx::createWorkbook()
  for (p in f) {
    sh <- sub("^Table[_ ]S", "Table S", sub("\\.csv$", "", basename(p), ignore.case = TRUE), ignore.case = TRUE)
    x <- utils::read.csv(p, header = FALSE, colClasses = "character", check.names = FALSE, na.strings = character(0), encoding = "UTF-8")
    openxlsx::addWorksheet(wb, sh)
    openxlsx::writeData(wb, sh, sh, startRow = 1)            # row 1 = caption row of the published workbook
    for (j in seq_along(x)) {
      v <- x[[j]]; nv <- suppressWarnings(as.numeric(v)); isn <- !is.na(nv) & nzchar(v)
      if (all(isn | !nzchar(v))) { v2 <- nv; v2[!nzchar(v)] <- NA; openxlsx::writeData(wb, sh, v2, startCol = j, startRow = 2, colNames = FALSE) }
      else {
        v2 <- v; v2[!nzchar(v2)] <- NA
        openxlsx::writeData(wb, sh, v2, startCol = j, startRow = 2, colNames = FALSE)
        if (any(isn) && sum(isn) < 2000) for (i in which(isn)) openxlsx::writeData(wb, sh, nv[i], startCol = j, startRow = 1 + i, colNames = FALSE)
      }
    }
  }
  openxlsx::saveWorkbook(wb, out, overwrite = TRUE)
  out
}
SUPP_XLSX <- Sys.getenv("MASLD_SUPP_XLSX", unset = "")
if (!nzchar(SUPP_XLSX) || !file.exists(SUPP_XLSX)) {
  .f <- list.files(c(project_dir, REPO_ROOT, file.path(REPO_ROOT, "data")), pattern = "^gupta_MASLD_SupplementaryTables.*\\.xlsx$", full.names = TRUE)
  SUPP_XLSX <- if (length(.f)) .f[1] else {
    .o <- file.path(results_dir, "cache", "Supplementary_Tables_from_csv.xlsx")
    if (file.exists(.o) || !is.null(tryCatch(build_supp_from_csv(.o), error = function(e) NULL))) .o else file.path(project_dir, "gupta_MASLD_SupplementaryTables.xlsx")
  }
}
message("Supplementary Tables: ", SUPP_XLSX)

# ---- Significance thresholds (Benjamini-Hochberg FDR) -------------------
FDR_GENE  <- 0.05   # DESeq2 genes (Table S1, Figs 1, 3, 6) and external replication (Tables S5-S7)
FDR_OMICS <- 0.10   # limma lipids and metabolites (Tables S3, S4, S10a, S10b, S11, S12, S19, S20)
FDR_KEGG  <- 0.10   # KEGG enrichment, at least 3 genes (Table S2b, Table S7)

# ---- Tables written by the code are collected in manuscript order ----------
SUPP <- list()
add_table <- function(name, x) {
  SUPP[[name]] <<- x
  if (is.data.frame(x)) {
    write.csv(x, file.path(results_dir, "tables", paste0(gsub(" ", "_", name), ".csv")), row.names = FALSE)
  } else for (lab in names(x)) if (is.data.frame(x[[lab]]))
    write.csv(x[[lab]], file.path(results_dir, "tables", paste0(gsub(" ", "_", name), "_", gsub("[^A-Za-z0-9]+", "_", lab), ".csv")), row.names = FALSE)
  invisible(x)
}

# ---- Source Data (Nature style: one workbook per figure, one sheet per panel) ---
# Written to <Figure folder>/Source_data/ ; SD_DIR is set by each figure script.
sd_panel <- function(name, desc, data) list(name = name, desc = desc, data = as.data.frame(data))
write_source_data <- function(file, panels) {
  dir.create(SD_DIR, showWarnings = FALSE, recursive = TRUE)
  wb <- openxlsx::createWorkbook()
  for (p in panels) {
    sh <- substr(p$name, 1, 31)
    openxlsx::addWorksheet(wb, sh)
    openxlsx::writeData(wb, sh, p$desc, startRow = 1)
    openxlsx::writeData(wb, sh, p$data, startRow = 3)
  }
  openxlsx::saveWorkbook(wb, file.path(SD_DIR, file), overwrite = TRUE)
  message("Wrote source data: ", file.path(SD_DIR, file))
}

# dotplot of an enrichment result; returns NULL (instead of an error) when nothing is enriched
safe_dotplot <- function(x, ...) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) return(invisible(NULL))
  enrichplot::dotplot(x, ...)
}

# ---- Running the figure scripts in dependency order -----------------------
# Analyses build on each other (e.g. Fig. 3 reuses the DESeq2 results of Fig. 1). Each script ends with
# save_state(<folder>); the next script calls load_prereq(<previous folder>), which loads that workspace,
# or runs the previous script first if it has not been run yet.
save_state <- function(tag) {
  save(list = ls(globalenv(), all.names = TRUE), file = file.path(results_dir, "cache", paste0(tag, ".RData")), envir = globalenv())
  message("Saved workspace: ", tag)
}
load_prereq <- function(tag) {
  f <- file.path(results_dir, "cache", paste0(tag, ".RData"))
  if (!file.exists(f)) {
    rmd <- file.path(REPO_ROOT, tag, paste0(tag, ".Rmd"))
    message("Running prerequisite: ", tag)
    r <- tempfile(fileext = ".R"); knitr::purl(rmd, output = r, quiet = TRUE)
    owd <- setwd(dirname(rmd)); on.exit(setwd(owd), add = TRUE)
    source(r, local = FALSE)
  } else {
    load(f, envir = globalenv()); message("Loaded workspace: ", tag)
  }
  invisible(TRUE)
}

# ---- Packages -----------------------------------------------------------------
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(gtsummary)
  library(ggplot2)
  library(tidyverse)
  library(RColorBrewer)
  library(POMA)
  library(ggtext)
  library(magrittr)
  library(tibble)
  library(cowplot)
  library(patchwork)
  library(ggpubr)
  library(ggstatsplot)
  library(mediation)
  library(limma)
  library(ggcorrplot)
  library(stats)
  library(corrplot)
  library(GO.db)
  library(GenomicFeatures)
  library(AnnotationDbi)
  library(DESeq2)
  library(dbplyr)
  library(dplyr)
  library(readr)
  library(gridExtra)
  library(grid)
  library(stringr)
  theme_set(theme_bw())
})
theme_set(theme_bw())
