# 01_olink_ttest.R — Pairwise differential expression via OlinkAnalyze::olink_ttest.

suppressPackageStartupMessages({
  library(OlinkAnalyze); library(dplyr); library(tibble)
})

# =============================================================================
# CONFIG — edit these
# =============================================================================
DATA_PATH <- "../data/your_npx_data.csv"   # long-format NPX with SampleID, Treatment, Assay, NPX
OUT_DIR   <- "."
GROUP_COL <- "Treatment"                   # name of the column that holds treatment/cohort labels

COMPARISONS <- list(
  GROUPA_vs_GROUPB = c("GroupA", "GroupB")
  # add more: e.g. GROUPC_vs_GROUPD = c("GroupC", "GroupD")
)

SUBSET_COMPARISONS <- list(
  # example:
  # SUBSET_vs_REST = list(parent = "GroupA", subset_label = "special",
  #                      subset_ids = c("SAMPLE_1","SAMPLE_2"))
)

# =============================================================================
# Pipeline
# =============================================================================
df <- read.csv(DATA_PATH, stringsAsFactors = FALSE)
df[[GROUP_COL]] <- trimws(gsub("\r", "", df[[GROUP_COL]]))   # robust to mixed line endings
df <- as_tibble(df)

run_pair <- function(d, g_pos, g_neg) {
  d2 <- d %>% filter(.data[[GROUP_COL]] %in% c(g_pos, g_neg))
  res <- olink_ttest(d2, variable = GROUP_COL)
  res$Threshold <- ifelse(res$Adjusted_pval < 0.05, "Significant", "Non-significant")
  as.data.frame(res)
}

run_subset_vs_rest <- function(d, parent, subset_label, subset_ids) {
  sub <- d %>% filter(.data[[GROUP_COL]] == parent)
  sub[[GROUP_COL]] <- ifelse(sub$SampleID %in% subset_ids, subset_label, "rest")
  res <- olink_ttest(sub, variable = GROUP_COL)
  res$Threshold <- ifelse(res$Adjusted_pval < 0.05, "Significant", "Non-significant")
  as.data.frame(res)
}

write_ttest <- function(tt_df, out_path) {
  colnames(tt_df) <- make.names(colnames(tt_df))
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  write.csv(tt_df, out_path, row.names = FALSE)
  cat(sprintf("Wrote %s  (n_sig=%d / %d)\n",
              out_path, sum(tt_df$Adjusted_pval < 0.05), nrow(tt_df)))
}

for (cmp_name in names(COMPARISONS)) {
  pair <- COMPARISONS[[cmp_name]]
  write_ttest(run_pair(df, pair[1], pair[2]),
              file.path(OUT_DIR, cmp_name, paste0(cmp_name, "_ttest_results.csv")))
}

for (cmp_name in names(SUBSET_COMPARISONS)) {
  cfg <- SUBSET_COMPARISONS[[cmp_name]]
  write_ttest(run_subset_vs_rest(df, cfg$parent, cfg$subset_label, cfg$subset_ids),
              file.path(OUT_DIR, cmp_name, paste0(cmp_name, "_ttest_results.csv")))
}
