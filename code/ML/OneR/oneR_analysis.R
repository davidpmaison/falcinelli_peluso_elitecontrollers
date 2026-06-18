# oneR_analysis.R — OneR single-feature classifier

suppressPackageStartupMessages({
  library(OneR); library(dplyr); library(tidyr); library(pROC)
})

# =============================================================================
# CONFIG — edit these
# =============================================================================
NPX_LONG_PATH <- "../../data/your_npx_long.csv"   # SampleID, Treatment, Assay, NPX
GROUP_COL     <- "Treatment"

COMPARISONS <- list(
  list(
    pos = "GroupA", neg = "GroupB",
    filter_fn = NULL, relabel_fn = NULL,
    out_csv = "posGroupA_negGroupB_machine_learning_model_performance.csv"
  )
)

K_FOLDS  <- 10
N_REPS   <- 10
N_BINS   <- 5
SEED     <- 42

# =============================================================================
# Pipeline
# =============================================================================
df <- read.csv(NPX_LONG_PATH, stringsAsFactors = FALSE)
df[[GROUP_COL]] <- trimws(gsub("\r", "", df[[GROUP_COL]]))

oneR_with_cv <- function(wide_df, pos_label, neg_label) {
  wide_df$group <- factor(wide_df$group, levels = c(pos_label, neg_label))
  od <- as.data.frame(wide_df[, -1])
  od <- od[, c(setdiff(colnames(od), "group"), "group")]
  bd <- bin(od, nbins = N_BINS, method = "content")

  m_full <- OneR(bd, verbose = FALSE); best_feat <- m_full$feature
  pred_full <- as.character(predict(m_full, bd))
  resub_acc <- mean(pred_full == as.character(bd$group)) * 100

  n <- nrow(bd); cv_acc <- c(); winners <- c()
  oof <- character(n); oof[] <- NA_character_
  for (rep_i in seq_len(N_REPS)) {
    set.seed(SEED + rep_i)
    p_i <- which(bd$group == pos_label); n_i <- which(bd$group == neg_label)
    folds <- integer(n)
    folds[p_i] <- sample(rep(seq_len(K_FOLDS), length.out = length(p_i)))
    folds[n_i] <- sample(rep(seq_len(K_FOLDS), length.out = length(n_i)))
    for (k in seq_len(K_FOLDS)) {
      tr <- bd[folds != k, ]; te <- bd[folds == k, ]
      m <- suppressWarnings(OneR(tr, verbose = FALSE)); winners <- c(winners, m$feature)
      p <- as.character(suppressWarnings(predict(m, te)))
      cv_acc <- c(cv_acc, mean(p == as.character(te$group), na.rm = TRUE))
      if (rep_i == N_REPS) oof[which(folds == k)] <- p
    }
  }

  truth <- as.character(od$group); keep <- !is.na(oof)
  tab <- table(truth = truth[keep], pred = oof[keep])
  for (lv in c(pos_label, neg_label)) {
    if (!lv %in% rownames(tab)) tab <- rbind(tab, setNames(rep(0, ncol(tab)), colnames(tab)))
    if (!lv %in% colnames(tab)) tab <- cbind(tab, setNames(rep(0, nrow(tab)), rownames(tab)))
  }
  tab <- tab[c(pos_label, neg_label), c(pos_label, neg_label)]
  tot <- sum(tab); obs <- sum(diag(tab))/tot
  ex <- sum(rowSums(tab) * colSums(tab))/tot^2
  kappa_val <- (obs - ex)/(1 - ex)
  TPR <- tab[pos_label, pos_label] / sum(tab[pos_label, ])
  FPR <- tab[neg_label, pos_label] / sum(tab[neg_label, ])
  prec <- if (sum(tab[, pos_label]) > 0) tab[pos_label, pos_label] / sum(tab[, pos_label]) else 0
  f1 <- if ((prec + TPR) > 0) 2 * prec * TPR / (prec + TPR) else 0
  auc_val <- as.numeric(pROC::auc(pROC::roc(od$group == pos_label, od[[best_feat]], quiet = TRUE)))

  list(
    summary = data.frame(
      Algorithm = "OneR",
      Accuracy = round(mean(cv_acc) * 100, 2),
      Kappa = round(kappa_val, 4),
      TP_Average = round(TPR, 3),
      FP_Average = round(FPR, 2),
      ROC_Area = round(auc_val, 3),
      F_Measure_Average = round(f1, 3)
    ),
    best_feature = best_feat,
    resub_acc = resub_acc,
    winning_features = sort(table(winners), decreasing = TRUE)
  )
}

for (cmp in COMPARISONS) {
  wide <- df %>%
    select(SampleID, !!sym(GROUP_COL), Assay, NPX) %>%
    pivot_wider(names_from = Assay, values_from = NPX, values_fn = mean)
  if (!is.null(cmp$filter_fn)) wide <- wide %>% filter(SampleID %in% cmp$filter_fn(wide))
  if (!is.null(cmp$relabel_fn)) {
    wide$group <- vapply(wide$SampleID, cmp$relabel_fn, character(1))
  } else {
    wide$group <- wide[[GROUP_COL]]
  }
  wide <- wide %>% filter(group %in% c(cmp$pos, cmp$neg))
  wide <- as.data.frame(wide[, c("SampleID","group",
                                  setdiff(colnames(wide), c("SampleID","group", GROUP_COL)))])
  res <- oneR_with_cv(wide, cmp$pos, cmp$neg)
  cat("\n", cmp$pos, "vs", cmp$neg, "— best feature:", res$best_feature, "\n")
  print(res$summary)
  write.csv(res$summary, cmp$out_csv, row.names = FALSE)
  cat("Wrote", cmp$out_csv, "\n")
}
