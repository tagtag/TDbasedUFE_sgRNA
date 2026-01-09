#!/usr/bin/env Rscript

# jacks_auc.R
# - Read precomputed JACKS gene effect tables from JACKS.zip
# - Compute ROC-AUC per cell line (standard benchmark style)
# - Summarize AUC distribution per dataset
#
# Usage:
#   Rscript jacks_auc.R JACKS.zip essential.txt nonessential.txt outdir
#
# Input assumptions:
#   - ZIP contains files like: JACKS/avana_gene_JACKS_results.txt
#   - Each *_gene_JACKS_results.txt is TSV:
#       column 1: Gene (gene symbol)
#       remaining columns: cell lines
#       values: gene effect (more negative = more essential)
#
# Notes for reproducibility:
#   - No randomness here, but we record sessionInfo().
#   - We force a consistent score direction: score = -effect (higher = more essential).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript jacks_auc.R JACKS.zip essential.txt nonessential.txt [outdir]\n")
  quit(status = 1)
}

zip_path <- args[1]
ess_path <- args[2]
non_path <- args[3]
outdir  <- ifelse(length(args) >= 4, args[4], "jacks_auc_out")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Packages ----
need_pkgs <- c("data.table", "pROC")
for (p in need_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(sprintf("Package '%s' is not installed. Install it via install.packages('%s').", p, p))
  }
}
library(data.table)

# ---- Read gene lists ----
read_gene_list <- function(path) {
  x <- fread(path, header = FALSE, sep = "\n", data.table = FALSE)[[1]]
  x <- trimws(as.character(x))
  x <- x[nzchar(x)]
  unique(x)
}

ess_genes <- read_gene_list(ess_path)
non_genes <- read_gene_list(non_path)

# sanity check
overlap <- intersect(ess_genes, non_genes)
if (length(overlap) > 0) {
  warning(sprintf("Essential/Nonessential lists overlap (%d genes). They will be removed from BOTH sets.", length(overlap)))
  ess_genes <- setdiff(ess_genes, overlap)
  non_genes <- setdiff(non_genes, overlap)
}

# ---- Unzip to temp dir ----
unz_dir <- file.path(outdir, "unzipped_jacks")
dir.create(unz_dir, showWarnings = FALSE, recursive = TRUE)
unzip(zip_path, exdir = unz_dir)

# ---- Find gene effect result files (exclude pval + rand by default) ----
all_files <- list.files(unz_dir, recursive = TRUE, full.names = TRUE)

is_gene_effect <- grepl("_gene_JACKS_results\\.txt$", basename(all_files), ignore.case = TRUE)
not_pval       <- !grepl("pval", basename(all_files), ignore.case = TRUE)
not_rand       <- !grepl("rand", basename(all_files), ignore.case = TRUE)

gene_files <- all_files[is_gene_effect & not_pval & not_rand]

if (length(gene_files) == 0) {
  stop("No *_gene_JACKS_results.txt files found (excluding pval/rand). Check ZIP contents.")
}

# ---- AUC computation helper ----
compute_auc_per_cellline <- function(dt, ess, non) {
  # dt: data.table with columns: Gene + cell lines
  stopifnot("Gene" %in% names(dt))
  genes <- dt[["Gene"]]

  keep <- genes %in% c(ess, non)
  dt2 <- dt[keep]

  if (nrow(dt2) == 0) {
    return(data.table(cell_line = character(), auc = numeric(), n_used = integer()))
  }

  y <- ifelse(dt2[["Gene"]] %in% ess, 1, 0)

  # cell lines = all but Gene
  cl_cols <- setdiff(names(dt2), "Gene")

  res <- lapply(cl_cols, function(cl) {
    effect <- dt2[[cl]]
    score  <- -effect  # higher = more essential (consistent direction)

    ok <- !is.na(score) & !is.na(y)
    yy <- y[ok]
    ss <- score[ok]

    # need both classes present
    if (length(unique(yy)) < 2) {
      return(list(cell_line = cl, auc = NA_real_, n_used = length(yy)))
    }

    # pROC: force direction so that controls(0) < cases(1) in score
    roc_obj <- pROC::roc(response = yy, predictor = ss, direction = "<", quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))

    list(cell_line = cl, auc = auc_val, n_used = length(yy))
  })

  rbindlist(res)
}

# ---- Main loop ----
all_auc <- list()
all_sum <- list()

for (f in gene_files) {
  dataset <- sub("_gene_JACKS_results\\.txt$", "", basename(f), ignore.case = TRUE)
  message(sprintf("[INFO] Reading %s", f))

  dt <- fread(f, sep = "\t", header = TRUE, data.table = TRUE)
  if (!("Gene" %in% names(dt))) {
    stop(sprintf("File %s does not have a 'Gene' column.", f))
  }

  # Compute AUC per cell line
  auc_dt <- compute_auc_per_cellline(dt, ess_genes, non_genes)
  auc_dt[, dataset := dataset]
  all_auc[[dataset]] <- auc_dt

  # Summary
  sum_dt <- auc_dt[!is.na(auc), .(
    n_cell_lines = .N,
    mean_auc     = mean(auc),
    median_auc   = median(auc),
    sd_auc       = sd(auc),
    min_auc      = min(auc),
    max_auc      = max(auc),
    mean_n_used  = mean(n_used)
  )]
  if (nrow(sum_dt) == 0) {
    sum_dt <- data.table(
      n_cell_lines = 0, mean_auc = NA_real_, median_auc = NA_real_, sd_auc = NA_real_,
      min_auc = NA_real_, max_auc = NA_real_, mean_n_used = NA_real_
    )
  }
  sum_dt[, dataset := dataset]
  all_sum[[dataset]] <- sum_dt
}

auc_out <- rbindlist(all_auc, use.names = TRUE, fill = TRUE)
sum_out <- rbindlist(all_sum, use.names = TRUE, fill = TRUE)

# ---- Write outputs ----
fwrite(auc_out, file.path(outdir, "JACKS_AUC_per_cellline.csv"))
fwrite(sum_out, file.path(outdir, "JACKS_AUC_summary_by_dataset.csv"))

# record session info for reproducibility
sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo())
sink()

message(sprintf("[DONE] Wrote:\n  - %s\n  - %s\n  - %s",
                file.path(outdir, "JACKS_AUC_per_cellline.csv"),
                file.path(outdir, "JACKS_AUC_summary_by_dataset.csv"),
                file.path(outdir, "sessionInfo.txt")))
