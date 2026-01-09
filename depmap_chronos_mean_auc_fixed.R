#!/usr/bin/env Rscript

# depmap_chronos_mean_auc_fixed.R
# DepMap (Chronos) mean->AUC:
#  1) (optional) restrict to ModelIDs listed in model_ids.txt
#  2) for each gene: mean gene effect across models
#  3) score is FIXED: score = mean_effect  (based on your observed medians)
#  4) ROC-AUC vs essential/nonessential lists

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("Usage: Rscript depmap_chronos_mean_auc_fixed.R CRISPRGeneEffect.csv essential.txt nonessential.txt outdir [model_ids.txt]\n")
  quit(status = 1)
}

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
})

depmap_csv <- args[1]
ess_path   <- args[2]
non_path   <- args[3]
outdir     <- args[4]
models_path <- ifelse(length(args) >= 5, args[5], NA_character_)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- helpers ----------
read_gene_list <- function(path) {
  x <- fread(path, header = FALSE, sep = "\n", data.table = FALSE, fill = TRUE)[[1]]
  x <- as.character(x)
  x <- trimws(x)
  x <- x[nzchar(x)]
  x <- sub("^\ufeff", "", x)              # BOM
  x <- x[!grepl("^#", x)]                 # comments
  x <- gsub('^"|"$', "", x)
  x <- gsub("^'|'$", "", x)
  x <- sub("[\t, ].*$", "", x)            # first token only
  x <- sub("\\s*\\(.*\\)$", "", x)        # drop " (xxxx)"
  x <- toupper(x)
  x <- x[!is.na(x) & x != "" & x != "NA"]
  unique(x)
}

col_to_symbol <- function(colname) {
  # "TP53 (7157)" -> "TP53"
  toupper(sub(" \\(.*\\)$", "", colname))
}

read_model_ids <- function(path) {
  if (is.na(path)) return(NULL)
  if (!file.exists(path)) stop(sprintf("model_ids.txt not found: %s", path))
  if (file.info(path)$size == 0) stop("model_ids.txt is empty (size 0). Create it first, or omit this argument.")
  v <- fread(path, header = FALSE, data.table = FALSE)[[1]]
  v <- trimws(as.character(v))
  v <- v[nzchar(v)]
  if (length(v) == 0) stop("model_ids.txt contained no valid IDs.")
  unique(v)
}

# ---------- read gene sets ----------
ess <- read_gene_list(ess_path)
non <- read_gene_list(non_path)

overlap <- intersect(ess, non)
if (length(overlap) > 0) {
  warning(sprintf("Essential/Nonessential lists overlap (%d genes). Removing from BOTH.", length(overlap)))
  ess <- setdiff(ess, overlap)
  non <- setdiff(non, overlap)
}
cat(sprintf("[INFO] essential=%d nonessential=%d\n", length(ess), length(non)))
wanted_genes <- sort(unique(c(ess, non)))

# ---------- optional model restriction ----------
models <- read_model_ids(models_path)
if (!is.null(models)) cat(sprintf("[INFO] restricting to %d ModelIDs\n", length(models)))

# ---------- read header and decide ID column + gene columns ----------
hdr <- fread(depmap_csv, nrows = 0)
cols <- names(hdr)

# include V1 for the common "blank header first column" case
id_candidates <- c("ModelID","DepMap_ID","Model_ID","depmap_id","model_id","V1")
id_col <- intersect(id_candidates, cols)
if (length(id_col) == 0) {
  id_col <- cols[1]
  warning(sprintf("Could not find ModelID/DepMap_ID/V1 column. Using first column: %s", id_col))
} else {
  id_col <- id_col[1]
}
cat(sprintf("[INFO] using ID column: %s\n", id_col))

gene_symbols <- col_to_symbol(cols)
gene_cols <- cols[gene_symbols %in% wanted_genes]
if (length(gene_cols) == 0) stop("No gene columns matched your essential/nonessential lists. Check symbols/file.")
cat(sprintf("[INFO] matched %d gene columns (from %d wanted genes)\n", length(gene_cols), length(wanted_genes)))

select_cols <- unique(c(id_col, gene_cols))

# ---------- load only required columns ----------
dt <- fread(depmap_csv, select = select_cols, showProgress = TRUE)

# Restrict models if provided
if (!is.null(models)) {
  dt[, (id_col) := as.character(get(id_col))]
  dt <- dt[get(id_col) %in% models]
}
n_models_used <- nrow(dt)
cat(sprintf("[INFO] models used: %d\n", n_models_used))
if (n_models_used == 0) stop("After restricting models, no rows remain. Check ID type and model_ids.txt.")

# Convert gene columns to numeric
dt[, (gene_cols) := lapply(.SD, as.numeric), .SDcols = gene_cols]

# ---------- compute mean effect per gene-column across models ----------
means_list <- dt[, lapply(.SD, mean, na.rm = TRUE), .SDcols = gene_cols]
means_vec <- as.numeric(means_list[1])
names(means_vec) <- names(means_list)

# Aggregate if multiple columns map to same symbol
sym <- col_to_symbol(names(means_vec))
gene_dt <- data.table(
  gene = sym,
  mean_effect = means_vec
)[, .(mean_effect = mean(mean_effect, na.rm = TRUE)), by = gene]

# Add labels
gene_dt[, label := fifelse(gene %in% ess, 1L,
                          fifelse(gene %in% non, 0L, NA_integer_))]
gene_dt <- gene_dt[!is.na(label) & is.finite(mean_effect)]

n_ess_used <- sum(gene_dt$label == 1L)
n_non_used <- sum(gene_dt$label == 0L)
cat(sprintf("[INFO] genes used: %d (essential=%d, nonessential=%d)\n",
            nrow(gene_dt), n_ess_used, n_non_used))
if (length(unique(gene_dt$label)) < 2) stop("Need both essential and nonessential genes after filtering.")

# ---------- FIXED score rule (as per your medians) ----------
# essentials have larger mean_effect than nonessentials in your file,
# so higher mean_effect should indicate "more essential".
gene_dt[, score := mean_effect]
score_rule <- "score = mean_effect (fixed)"

med_ess <- median(gene_dt[label == 1L, mean_effect], na.rm = TRUE)
med_non <- median(gene_dt[label == 0L, mean_effect], na.rm = TRUE)
cat(sprintf("[INFO] median(mean_effect): essential=%.6f nonessential=%.6f\n", med_ess, med_non))
cat(sprintf("[INFO] using: %s\n", score_rule))

# ---------- ROC-AUC ----------
# direction="<" means controls(0) < cases(1) in score (i.e., higher score => essential)
roc_obj <- pROC::roc(response = gene_dt$label, predictor = gene_dt$score,
                     levels = c(0, 1), direction = "<", quiet = TRUE)
auc_val <- as.numeric(pROC::auc(roc_obj))

cat(sprintf("[DONE] AUC = %.6f\n", auc_val))

# ---------- outputs ----------
fwrite(gene_dt[, .(gene, mean_effect, score, label)],
       file.path(outdir, "DepMap_Chronos_mean_to_AUC_gene_means.csv"))

summary_dt <- data.table(
  auc = auc_val,
  score_rule = score_rule,
  median_mean_effect_essential = med_ess,
  median_mean_effect_nonessential = med_non,
  n_models_used = n_models_used,
  n_gene_columns_matched = length(gene_cols),
  n_genes_used_after_merge = nrow(gene_dt),
  n_essential_used = n_ess_used,
  n_nonessential_used = n_non_used,
  id_column = id_col,
  requested_model_ids = ifelse(is.null(models), NA_integer_, length(models))
)
fwrite(summary_dt, file.path(outdir, "DepMap_Chronos_mean_to_AUC_summary.csv"))

sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo())
sink()
