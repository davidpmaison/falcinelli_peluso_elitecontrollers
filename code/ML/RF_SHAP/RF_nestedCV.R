## RF_nestedCV.R -- Nested-CV Random Forest + SHAP for HIV elite controllers
## ----------------------------------------------------------------------------
## setwd("/path/to/repo") 
## ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ranger); library(pROC); library(parallel)
  library(fastshap)
})

.this_file <- local({
  a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(tryCatch(normalizePath(sub("^--file=", "", m[1])), error = function(e) ""))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    return(tryCatch(normalizePath(rstudioapi::getSourceEditorContext()$path), error = function(e) ""))
  ""
})
if (nzchar(.this_file)) {
  setwd(dirname(.this_file))
  for (up in c(".", "..", "../..", "../../..", "../../../..")) {
    if (dir.exists(file.path(up, "raw_data"))) { setwd(up); break }
  }
}

## ====================================================
DATA  <- "raw_data"
OUT   <- "."
RF_DIR    <- file.path(OUT, "results", "ML", "RF_SHAP")
SWEEP_DIR <- file.path(RF_DIR, "sweep_predictions")
MANIFEST  <- file.path(OUT, "results", "RUN_COMPLETE")

SEED      <- 42
OUTER_K   <- 10; OUTER_R <- 10; INNER_K <- 10
SENS_K    <- 5            
NTREE     <- 1000
N_BOOT    <- 1000
TOP_NS    <- c(10, 15, 20)
SHAP_NSIM <- 50
N_FOLDS   <- OUTER_K * OUTER_R   # 100 outer folds total

STRICT_SHAP <- TRUE       

USE_PARALLEL <- TRUE
NCORES       <- if (USE_PARALLEL) max(1, parallel::detectCores() - 2) else 1L
PAR_BACKEND  <- if (!USE_PARALLEL || NCORES <= 1) "serial" else
                if (.Platform$OS.type == "windows") "psock" else "fork"

min_node_grid    <- c(1, 3, 5, 10)
max_depth_grid   <- c(0, 5, 10)
sample_frac_grid <- c(1.0, 0.8, 0.6)
mtry_grid        <- function(p) unique(pmax(1, c(round(sqrt(p)), round(p / 3), round(p / 2), p)))

EEC_ALL <- c("EC02","EC08","EC12","EC19","EC41","EC46","EC55","EC57","EC58")


mergedfile <- file.path(DATA, "MergedData.csv")
demofile   <- file.path(DATA, "Olink_demographics_20230526.csv")
if (!file.exists(mergedfile) || !file.exists(demofile))
  stop("Required data file(s) not found in '", DATA, "':\n  ",
       mergedfile, "\n  ", demofile,
       "\nRun this script from the repo root (the folder that contains raw_data/).",
       call. = FALSE)

unlink(RF_DIR, recursive = TRUE)
unlink(file.path(OUT, "results", "run_config.csv"))
unlink(MANIFEST)
dir.create(SWEEP_DIR, showWarnings = FALSE, recursive = TRUE)


average_precision <- function(truth, score) {
  if (sum(truth) == 0) return(NA_real_)
  o <- order(score, decreasing = TRUE); t <- truth[o]
  tp <- cumsum(t); fp <- cumsum(1 - t)
  sum(diff(c(0, tp / sum(t))) * (tp / (tp + fp)))
}

ap_ci <- function(truth, score, n_boot, seed) {
  set.seed(seed); n <- length(truth)
  aps <- vapply(seq_len(n_boot), function(i) {
    idx <- sample.int(n, replace = TRUE)
    if (length(unique(truth[idx])) < 2) return(NA_real_)
    average_precision(truth[idx], score[idx])
  }, numeric(1))
  as.numeric(stats::quantile(aps, c(0.025, 0.975), na.rm = TRUE))
}

rank_features_ttest <- function(wide_train, protein_cols) {
  y <- wide_train$group
  pvals <- vapply(protein_cols, function(p) {
    x <- wide_train[[p]]
    a <- x[y == "pos"]; b <- x[y == "neg"]
    a <- a[!is.na(a)]; b <- b[!is.na(b)]
    if (length(a) < 2 || length(b) < 2) return(1)
    suppressWarnings(tryCatch(t.test(a, b)$p.value, error = function(e) 1))
  }, numeric(1))
  padj <- p.adjust(pvals, method = "BH")
  protein_cols[order(padj, pvals)]
}

impute_train_median <- function(train, test, cols) {
  meds <- sapply(cols, function(c) median(train[[c]], na.rm = TRUE))
  for (c in cols) {
    train[[c]][is.na(train[[c]])] <- meds[[c]]
    test[[c]][is.na(test[[c]])]   <- meds[[c]]
  }
  list(train = train, test = test)
}

eval_combo <- function(combo, tr_x, tr_y, inner_folds) {
  ufolds <- unique(inner_folds); auc_inner <- numeric(length(ufolds))
  for (ii in seq_along(ufolds)) {
    ifold <- ufolds[ii]
    itr <- which(inner_folds != ifold); ite <- which(inner_folds == ifold)
    fit <- ranger(x = tr_x[itr, , drop = FALSE], y = tr_y[itr], num.trees = NTREE,
                  mtry = combo$mtry, min.node.size = combo$min.node.size,
                  max.depth = combo$max.depth, sample.fraction = combo$sample.fraction,
                  probability = TRUE, classification = TRUE, num.threads = 1, seed = combo$seed)
    p_pos <- predict(fit, data = tr_x[ite, , drop = FALSE], num.threads = 1)$predictions[, "pos"]
    truth <- tr_y[ite] == "pos"
    auc_inner[ii] <- if (length(unique(truth)) < 2) NA_real_ else
      as.numeric(suppressMessages(pROC::auc(pROC::roc(truth, p_pos, quiet = TRUE, direction = "<"))))
  }
  mean(auc_inner, na.rm = TRUE)
}

par_eval <- function(combos, tr_x, tr_y, inner_folds) {
  if (PAR_BACKEND == "fork") {
    unlist(parallel::mclapply(seq_len(nrow(combos)), function(i)
      eval_combo(as.list(combos[i, ]), tr_x, tr_y, inner_folds), mc.cores = NCORES))
  } else if (PAR_BACKEND == "psock") {
    unlist(parallel::parLapply(.CL, seq_len(nrow(combos)),
      function(i, combos, tr_x, tr_y, inner_folds)
        eval_combo(as.list(combos[i, ]), tr_x, tr_y, inner_folds),
      combos = combos, tr_x = tr_x, tr_y = tr_y, inner_folds = inner_folds))
  } else {
    vapply(seq_len(nrow(combos)), function(i)
      eval_combo(as.list(combos[i, ]), tr_x, tr_y, inner_folds), numeric(1))
  }
}

shap_pred_fun <- function(object, newdata)
  predict(object, data = newdata, num.threads = 1)$predictions[, "pos"]

run_full_leak_free <- function(wide, top_n, label, k_outer = OUTER_K, k_inner = INNER_K) {
  protein_cols <- setdiff(colnames(wide),
                          c("SampleID","Treatment","group","Intact_10E6_FLIP","CD4abs"))
  n <- nrow(wide); truth_vec <- ifelse(wide$group == "pos", 1, 0)
  pred_sum <- numeric(n); pred_cnt <- integer(n)
  best_params_log <- list(); feat_picked <- list()
  t_start <- Sys.time()
  for (rep_i in seq_len(OUTER_R)) {
    set.seed(SEED + rep_i)
    pii <- which(wide$group == "pos"); nii <- which(wide$group == "neg")
    f <- integer(n)
    f[pii] <- sample(rep(seq_len(k_outer), length.out = length(pii)))
    f[nii] <- sample(rep(seq_len(k_outer), length.out = length(nii)))
    for (k in seq_len(k_outer)) {
      tr_idx <- which(f != k); te_idx <- which(f == k)
      tr <- wide[tr_idx, ]; te <- wide[te_idx, ]
      ranks <- rank_features_ttest(tr, protein_cols)
      top_feats <- ranks[seq_len(min(top_n, length(ranks)))]
      feat_picked[[length(feat_picked) + 1]] <- top_feats
      use_cols <- c(top_feats, "Intact_10E6_FLIP", "CD4abs")
      imp <- impute_train_median(tr[, use_cols, drop = FALSE], te[, use_cols, drop = FALSE], use_cols)
      tr_x <- imp$train; te_x <- imp$test
      tr_y <- factor(tr$group, levels = c("neg","pos"))
      set.seed(SEED + rep_i * 100 + k)
      inner_folds <- integer(nrow(tr_x))
      ip <- which(tr_y == "pos"); in_ <- which(tr_y == "neg")
      inner_folds[ip]  <- sample(rep(seq_len(k_inner), length.out = length(ip)))
      inner_folds[in_] <- sample(rep(seq_len(k_inner), length.out = length(in_)))
      p <- ncol(tr_x)
      combos <- expand.grid(mtry = mtry_grid(p), min.node.size = min_node_grid,
                            max.depth = max_depth_grid, sample.fraction = sample_frac_grid,
                            KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      combos <- unique(combos); combos$seed <- SEED + rep_i * 1000 + k * 10
      aucs <- par_eval(combos, tr_x, tr_y, inner_folds)
      best <- as.list(combos[which.max(aucs), ])
      best$mean_inner_auc <- max(aucs, na.rm = TRUE); best$rep <- rep_i; best$fold <- k
      best_params_log[[length(best_params_log) + 1]] <- best
      fit <- ranger(x = tr_x, y = tr_y, num.trees = NTREE,
                    mtry = best$mtry, min.node.size = best$min.node.size,
                    max.depth = best$max.depth, sample.fraction = best$sample.fraction,
                    probability = TRUE, classification = TRUE, num.threads = NCORES, seed = best$seed)
      p_pos <- predict(fit, data = te_x, num.threads = NCORES)$predictions[, "pos"]
      pred_sum[te_idx] <- pred_sum[te_idx] + p_pos
      pred_cnt[te_idx] <- pred_cnt[te_idx] + 1L
    }
    cat(sprintf("  rep %d/%d done (cumulative %.1f min)\n", rep_i, OUTER_R,
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  }
  pred_mean <- pred_sum / pred_cnt
  roc_obj <- pROC::roc(truth_vec, pred_mean, quiet = TRUE, direction = "<")
  auc_val <- as.numeric(pROC::auc(roc_obj))
  set.seed(SEED)                                   # reproducible bootstrap CI
  ci   <- as.numeric(pROC::ci.auc(roc_obj, conf.level = 0.95,
                                  method = "bootstrap", boot.n = N_BOOT))
  ap   <- average_precision(truth_vec, pred_mean)
  ap_c <- ap_ci(truth_vec, pred_mean, N_BOOT, SEED)
  list(label = label, n_pos = sum(wide$group == "pos"), n_neg = sum(wide$group == "neg"),
       top_n = top_n, prevalence = mean(truth_vec),
       auc = auc_val, auc_ci_lower = ci[1], auc_ci_upper = ci[3],
       aupr = ap,    aupr_ci_lower = ap_c[1], aupr_ci_upper = ap_c[2],
       pred_mean = pred_mean, truth = truth_vec, sample_ids = wide$SampleID,
       best_params_log = best_params_log,
       feat_freq = sort(table(unlist(feat_picked)), decreasing = TRUE))
}

pool_shap <- function(wide, top_n, best_params_log, k_outer = OUTER_K) {
  bp <- do.call(rbind, lapply(best_params_log, as.data.frame))
  protein_cols <- setdiff(colnames(wide),
                          c("SampleID","Treatment","group","Intact_10E6_FLIP","CD4abs"))
  n <- nrow(wide); long <- list(); n_fail <- 0L; fail_ids <- character()
  for (rep_i in seq_len(OUTER_R)) {
    set.seed(SEED + rep_i)                          # same outer folds as run_full_leak_free
    pii <- which(wide$group == "pos"); nii <- which(wide$group == "neg")
    f <- integer(n)
    f[pii] <- sample(rep(seq_len(k_outer), length.out = length(pii)))
    f[nii] <- sample(rep(seq_len(k_outer), length.out = length(nii)))
    for (k in seq_len(k_outer)) {
      tr_idx <- which(f != k); te_idx <- which(f == k)
      tr <- wide[tr_idx, ]; te <- wide[te_idx, ]
      top_feats <- rank_features_ttest(tr, protein_cols)[seq_len(min(top_n, length(protein_cols)))]
      use_cols <- c(top_feats, "Intact_10E6_FLIP", "CD4abs")
      imp <- impute_train_median(tr[, use_cols, drop = FALSE], te[, use_cols, drop = FALSE], use_cols)
      tr_x <- imp$train; te_x <- imp$test
      tr_y <- factor(tr$group, levels = c("neg","pos"))
      hp <- bp[bp$rep == rep_i & bp$fold == k, ][1, ]
      fit <- ranger(x = tr_x, y = tr_y, num.trees = NTREE,
                    mtry = hp$mtry, min.node.size = hp$min.node.size,
                    max.depth = hp$max.depth, sample.fraction = hp$sample.fraction,
                    probability = TRUE, classification = TRUE, num.threads = 1,
                    seed = SEED + rep_i * 1000 + k * 10)
      set.seed(SEED + rep_i * 1000 + k * 10)        # reproducible Monte-Carlo draws
      sv <- tryCatch(fastshap::explain(fit, X = tr_x, newdata = te_x,
                                       pred_wrapper = shap_pred_fun, nsim = SHAP_NSIM),
                     error = function(e) { message("SHAP fail: ", e$message); NULL })
      if (is.null(sv)) {
        n_fail <- n_fail + 1L
        fail_ids <- c(fail_ids, sprintf("rep%d.k%d", rep_i, k))
        next
      }
      long[[length(long) + 1]] <- do.call(rbind, lapply(use_cols, function(ff)
        data.frame(SampleID = te$SampleID, feature = ff, shap = as.numeric(sv[, ff]),
                   raw = te_x[[ff]], rep = rep_i, fold = k, stringsAsFactors = FALSE)))
    }
    cat(sprintf("  SHAP rep %d/%d\n", rep_i, OUTER_R))
  }
  n_expected <- k_outer * OUTER_R
  if (n_fail > 0) {
    msg <- sprintf("SHAP failed for %d/%d folds (%s).",
                   n_fail, n_expected, paste(fail_ids, collapse = ", "))
    if (STRICT_SHAP)
      stop(msg, " STRICT_SHAP=TRUE: aborting before the completion manifest so this run is not marked complete.",
           call. = FALSE)
    warning(msg, " rankings use incomplete coverage (STRICT_SHAP=FALSE).")
  }
  do.call(rbind, long)
}

shap_importance <- function(shap_long, pos, feat_stability) {
  tot     <- length(unique(shap_long$SampleID)) * length(unique(shap_long$rep))
  ## actual fold count covered by this shap_long (100 for primary 10-fold;
  ## SENS_K * OUTER_R for the sensitivity rerun) -- used to label forced
  ## clinical covariates correctly in the SHAP_importance table
  n_folds <- length(unique(paste(shap_long$rep, shap_long$fold)))
  agg <- do.call(rbind, lapply(split(shap_long, shap_long$feature), function(g) {
    r <- suppressWarnings(cor(g$raw, g$shap, use = "complete.obs"))
    data.frame(feature = g$feature[1],
               mean_abs_SHAP = round(sum(abs(g$shap), na.rm = TRUE) / tot, 6),
               direction = if (is.na(r)) NA
                           else if (r < 0) paste0("lower in ", pos)
                           else paste0("higher in ", pos),
               stringsAsFactors = FALSE)
  }))
  agg <- agg[order(-agg$mean_abs_SHAP), ]; agg$rank <- seq_len(nrow(agg))
  agg <- merge(agg, feat_stability[, c("feature","n_folds_selected","pct_of_folds")],
               by = "feature", all.x = TRUE)
  clin <- agg$feature %in% c("Intact_10E6_FLIP","CD4abs")     # covariates always included (forced features)
  agg$n_folds_selected[clin] <- n_folds; agg$pct_of_folds[clin] <- 100
  agg[order(agg$rank), c("rank","feature","mean_abs_SHAP","direction",
                         "n_folds_selected","pct_of_folds")]
}

## ====================================================
.CL <- NULL
if (PAR_BACKEND == "psock") {
  .CL <- parallel::makeCluster(NCORES)
  parallel::clusterEvalQ(.CL, suppressPackageStartupMessages({
    library(ranger); library(pROC) }))
  parallel::clusterExport(.CL, c("eval_combo", "NTREE"))
  reg.finalizer(globalenv(), function(e) {
    if (!is.null(.CL)) try(parallel::stopCluster(.CL), silent = TRUE)
  }, onexit = TRUE)
}
cat(sprintf("Parallel backend: %s (%d core%s)\n",
            PAR_BACKEND, NCORES, if (NCORES == 1) "" else "s"))

df_long <- read.csv(mergedfile, stringsAsFactors = FALSE)
df_long$Treatment <- trimws(gsub("\r", "", df_long$Treatment))
demo <- read.csv(demofile, stringsAsFactors = FALSE)
for (v in c("Intact_10E6_FLIP","CD4abs"))
  demo[[v]] <- suppressWarnings(as.numeric(demo[[v]]))

wide_ea <- df_long %>% filter(Treatment %in% c("EC","ART")) %>%
  select(SampleID, Treatment, Assay, NPX) %>%
  pivot_wider(names_from = Assay, values_from = NPX, values_fn = mean) %>% as.data.frame()
wide_ea$group <- ifelse(wide_ea$Treatment == "EC", "pos", "neg")
wide_ea <- merge(wide_ea, demo %>% select(SampleID, Intact_10E6_FLIP, CD4abs),
                 by = "SampleID", all.x = TRUE)

wide_ec <- df_long %>% filter(Treatment == "EC") %>%
  select(SampleID, Treatment, Assay, NPX) %>%
  pivot_wider(names_from = Assay, values_from = NPX, values_fn = mean) %>% as.data.frame()
wide_ec$group <- ifelse(wide_ec$SampleID %in% EEC_ALL, "pos", "neg")
wide_ec <- merge(wide_ec, demo %>% select(SampleID, Intact_10E6_FLIP, CD4abs),
                 by = "SampleID", all.x = TRUE)

cells <- list(
  list(comp = "EEC_vs_TEC", folder = "EEC_TEC", wide = wide_ec, label = "EEC vs TEC", pos = "EEC"),
  list(comp = "EC_vs_ART",  folder = "EC_ART",  wide = wide_ea, label = "EC vs ART",  pos = "EC")
)

## ====================================================
summary_path <- file.path(SWEEP_DIR, "topn_sweep_full_tuning_AUC.csv")
store <- list(); rows <- list()
for (cell in cells) {
  store[[cell$comp]] <- list()
  for (tn in TOP_NS) {
    res <- run_full_leak_free(cell$wide, top_n = tn,
                              label = sprintf("%s (top %d + clinical)", cell$label, tn))
    store[[cell$comp]][[as.character(tn)]] <- res
    write.csv(data.frame(SampleID = res$sample_ids, truth = res$truth, rf_prob_oof = res$pred_mean),
              file.path(SWEEP_DIR, sprintf("%s_top%d_predictions_oof.csv", cell$comp, tn)),
              row.names = FALSE)
    rows[[length(rows) + 1]] <- data.frame(
      comparison = cell$label, n_pos = res$n_pos, n_neg = res$n_neg, top_n = tn,
      n_features = tn + 2L, feat_per_positive = round((tn + 2) / res$n_pos, 2),
      full_tuning_AUC = round(res$auc, 4), ci_lower = round(res$auc_ci_lower, 4),
      ci_upper = round(res$auc_ci_upper, 4),
      ci_width = round(res$auc_ci_upper - res$auc_ci_lower, 4),
      top_feat = names(res$feat_freq)[1],
      prevalence  = round(res$prevalence, 3),
      auPR        = round(res$aupr, 4),
      pr_ci_lower = round(res$aupr_ci_lower, 4),
      pr_ci_upper = round(res$aupr_ci_upper, 4),
      stringsAsFactors = FALSE)
    write.csv(do.call(rbind, rows), summary_path, row.names = FALSE)
    cat(sprintf(">>> %s top %d: auROC=%.4f (%.3f-%.3f)  auPR=%.4f (%.3f-%.3f)\n",
                cell$label, tn, res$auc, res$auc_ci_lower, res$auc_ci_upper,
                res$aupr, res$aupr_ci_lower, res$aupr_ci_upper))
  }
}
sweep <- do.call(rbind, rows)

## ====================================================
selected_tn <- list(); imp_tables <- list()
for (cell in cells) {
  cr <- sweep[sweep$comparison == cell$label, ]
  ## use raw (unrounded) values from store so rounding can't flip the selection
  raw_aucs   <- vapply(cr$top_n, function(tn)
                       store[[cell$comp]][[as.character(tn)]]$auc, numeric(1))
  raw_widths <- vapply(cr$top_n, function(tn) {
    rr <- store[[cell$comp]][[as.character(tn)]]
    rr$auc_ci_upper - rr$auc_ci_lower
  }, numeric(1))
  best_tn <- cr$top_n[order(-raw_aucs, raw_widths)][1]
  selected_tn[[cell$comp]] <- best_tn
  res  <- store[[cell$comp]][[as.character(best_tn)]]
  fdir <- file.path(RF_DIR, cell$folder); dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
  base <- file.path(fdir, cell$comp)

  write.csv(data.frame(SampleID = res$sample_ids, truth = res$truth, rf_prob_oof = res$pred_mean),
            paste0(base, "_predictions_oof.csv"), row.names = FALSE)
  fs <- data.frame(feature = names(res$feat_freq),
                   n_folds_selected = as.integer(res$feat_freq),
                   pct_of_folds = round(as.integer(res$feat_freq) / N_FOLDS * 100, 1))
  write.csv(fs, paste0(base, "_feature_stability.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(res$best_params_log, as.data.frame)),
            paste0(base, "_best_hyperparams_per_fold.csv"), row.names = FALSE)

  shap_long <- pool_shap(cell$wide, best_tn, res$best_params_log)
  write.csv(shap_long, paste0(base, "_SHAP_long.csv"), row.names = FALSE)
  imp <- shap_importance(shap_long, cell$pos, fs)
  imp_tables[[cell$comp]] <- imp     # keep for cross-K SHAP concordance below
  write.csv(imp, paste0(base, "_SHAP_importance.csv"), row.names = FALSE)
  cat(sprintf("  %s: selected top %d (auROC %.4f, auPR %.4f); SHAP #1 = %s\n",
              cell$label, best_tn, res$auc, res$aupr, imp$feature[1]))
}

sens_rows <- list(); sens_sweep_rows <- list()
for (cell in cells) {
  if (min(sum(cell$wide$group == "pos"),
          sum(cell$wide$group == "neg")) >= OUTER_K) next
  sres <- list()
  for (tn in TOP_NS) {
    rr <- run_full_leak_free(cell$wide, top_n = tn,
                             label = sprintf("%s SENS %d-fold top %d", cell$label, SENS_K, tn),
                             k_outer = SENS_K, k_inner = SENS_K)
    sres[[as.character(tn)]] <- rr
    sens_sweep_rows[[length(sens_sweep_rows) + 1]] <- data.frame(
      comparison = cell$label, fold_K = SENS_K, top_n = tn,
      auROC = round(rr$auc, 4),
      roc_ci_lower = round(rr$auc_ci_lower, 4),
      roc_ci_upper = round(rr$auc_ci_upper, 4),
      auPR = round(rr$aupr, 4),
      pr_ci_lower = round(rr$aupr_ci_lower, 4),
      pr_ci_upper = round(rr$aupr_ci_upper, 4),
      stringsAsFactors = FALSE)
  }
  ss <- do.call(rbind, lapply(TOP_NS, function(tn) {
    rr <- sres[[as.character(tn)]]
    data.frame(top_n = tn, auc = rr$auc, ci_width = rr$auc_ci_upper - rr$auc_ci_lower)
  }))
  sens_sel <- ss$top_n[order(-ss$auc, ss$ci_width)][1]
  prim_tn  <- selected_tn[[cell$comp]]
  prim     <- store[[cell$comp]][[as.character(prim_tn)]]
  at_prim  <- sres[[as.character(prim_tn)]]
  sens_res <- sres[[as.character(sens_sel)]]
  d_roc    <- abs(prim$auc  - at_prim$auc)
  d_pr     <- abs(prim$aupr - at_prim$aupr)
  selection_ok <- (prim_tn == sens_sel)

  ## ====================================================
  shap_long_sens <- pool_shap(cell$wide, sens_sel, sens_res$best_params_log,
                              k_outer = SENS_K)
  sens_fs <- data.frame(feature = names(sens_res$feat_freq),
                        n_folds_selected = as.integer(sens_res$feat_freq),
                        pct_of_folds = round(as.integer(sens_res$feat_freq) /
                                                  (SENS_K * OUTER_R) * 100, 1))
  sens_imp <- shap_importance(shap_long_sens, cell$pos, sens_fs)
  sens_dir <- file.path(RF_DIR, sprintf("fold_sensitivity_K%d", SENS_K), cell$folder)
  dir.create(sens_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(shap_long_sens,
            file.path(sens_dir, sprintf("%s_SHAP_long.csv", cell$comp)),
            row.names = FALSE)
  write.csv(sens_imp,
            file.path(sens_dir, sprintf("%s_SHAP_importance.csv", cell$comp)),
            row.names = FALSE)

  ## ====================================================
  prim_imp  <- imp_tables[[cell$comp]]
  col_k5    <- sprintf("rank_K%d", SENS_K)
  col_k10   <- sprintf("rank_K%d", OUTER_K)
  both <- merge(sens_imp[, c("feature","rank","mean_abs_SHAP","direction")],
                prim_imp[, c("feature","rank","mean_abs_SHAP","direction")],
                by = "feature", all = TRUE,
                suffixes = c(sprintf("_K%d", SENS_K), sprintf("_K%d", OUTER_K)))
  both <- both[order(pmin(both[[col_k5]], both[[col_k10]], na.rm = TRUE)), ]
  write.csv(both, file.path(sens_dir,
            sprintf("%s_SHAP_compare_K%d_vs_K%d.csv", cell$comp, SENS_K, OUTER_K)),
            row.names = FALSE)
  common <- both[!is.na(both[[col_k5]]) & !is.na(both[[col_k10]]), ]
  spearman_rho <- if (nrow(common) >= 3)
    suppressWarnings(cor(common[[col_k5]], common[[col_k10]], method = "spearman"))
    else NA_real_
  ovl <- function(k) length(intersect(utils::head(sens_imp$feature, k),
                                       utils::head(prim_imp$feature, k)))
  ovl5 <- ovl(5); ovl10 <- ovl(10); ovl20 <- ovl(20)

  ## ====================================================
  decision <- if (d_roc >= 0.1)
                sprintf("REVIEW: |d auROC| = %.3f (>= 0.1)", d_roc)
              else if (!selection_ok)
                sprintf("performance stable but feature-count selection unstable (K%d top%d vs primary K%d top%d)",
                        SENS_K, sens_sel, OUTER_K, prim_tn)
              else
                sprintf("keep %d-fold (auROC stable and top_n agrees)", OUTER_K)

  sens_rows[[length(sens_rows) + 1]] <- data.frame(
    comparison = cell$label, primary_K = OUTER_K, sens_K = SENS_K,
    primary_selected_top_n = prim_tn,
    primary_auROC = round(prim$auc, 4), primary_auPR = round(prim$aupr, 4),
    sens_selected_top_n = sens_sel,
    sens_auROC_at_primary_topn = round(at_prim$auc, 4),
    sens_auPR_at_primary_topn  = round(at_prim$aupr, 4),
    sens_auROC_at_sens_topn    = round(sens_res$auc, 4),
    sens_auPR_at_sens_topn     = round(sens_res$aupr, 4),
    delta_auROC_matched_topn   = round(d_roc, 4),
    delta_auPR_matched_topn    = round(d_pr, 4),
    selection_agrees           = selection_ok,
    shap_spearman_rho          = if (is.na(spearman_rho)) NA_real_ else round(spearman_rho, 3),
    shap_top5_overlap          = sprintf("%d/5",  ovl5),
    shap_top10_overlap         = sprintf("%d/10", ovl10),
    shap_top20_overlap         = sprintf("%d/20", ovl20),
    decision                   = decision,
    stringsAsFactors = FALSE)

  cat(sprintf("  K%d sel top%d (auROC %.3f, auPR %.3f) | K%d sel top%d; K%d@top%d auROC %.3f (auPR %.3f) | matched |d auROC|=%.3f -> %s\n",
              OUTER_K, prim_tn, prim$auc, prim$aupr,
              SENS_K, sens_sel,
              SENS_K, prim_tn, at_prim$auc, at_prim$aupr,
              d_roc, decision))
  cat(sprintf("  SHAP concordance: Spearman rho = %s; overlap top5 %d/5, top10 %d/10, top20 %d/20\n",
              if (is.na(spearman_rho)) "NA" else sprintf("%.3f", spearman_rho),
              ovl5, ovl10, ovl20))
}
sens_df <- if (length(sens_rows)) do.call(rbind, sens_rows) else NULL
if (!is.null(sens_df))
  write.csv(sens_df, file.path(RF_DIR, "fold_sensitivity.csv"), row.names = FALSE)
if (length(sens_sweep_rows))
  write.csv(do.call(rbind, sens_sweep_rows),
            file.path(RF_DIR, "fold_sensitivity_sweep.csv"), row.names = FALSE)

## ====================================================
cfg <- c(
  seed = SEED, outer_repeats = OUTER_R, outer_K = OUTER_K, inner_K = INNER_K,
  ntree = NTREE, n_boot = N_BOOT, shap_nsim = SHAP_NSIM,
  top_ns = paste(TOP_NS, collapse = "/"),
  imputation_method = "train_median",
  strict_shap = STRICT_SHAP,
  par_backend = PAR_BACKEND, ncores = NCORES,
  sens_K = SENS_K
)
for (cell in cells)
  cfg[paste0("selected_top_n_", cell$comp)] <- selected_tn[[cell$comp]]
if (!is.null(sens_df)) for (i in seq_len(nrow(sens_df))) {
  key <- paste0("fold_sensitivity_", gsub(" ", "_", sens_df$comparison[i]))
  cfg[key] <- sprintf(
    "%d-fold(top%d) auROC %.3f vs %d-fold(top%d) auROC %.3f; matched-top_n |d|=%.3f; selection_agrees=%s; %s",
    sens_df$primary_K[i], sens_df$primary_selected_top_n[i], sens_df$primary_auROC[i],
    sens_df$sens_K[i], sens_df$sens_selected_top_n[i], sens_df$sens_auROC_at_sens_topn[i],
    sens_df$delta_auROC_matched_topn[i], sens_df$selection_agrees[i], sens_df$decision[i])
}
write.csv(data.frame(key = names(cfg), value = unname(cfg), stringsAsFactors = FALSE),
          file.path(OUT, "results", "run_config.csv"), row.names = FALSE)

if (!is.null(.CL)) parallel::stopCluster(.CL)

## ====================================================
writeLines(c(
  "status: complete",
  sprintf("completed_at: %s", format(Sys.time())),
  sprintf("seed: %d", SEED),
  sprintf("outer_repeats: %d", OUTER_R),
  sprintf("outer_K: %d", OUTER_K),
  sprintf("inner_K: %d", INNER_K),
  sprintf("top_ns: %s", paste(TOP_NS, collapse = "/")),
  sprintf("sens_K: %d", SENS_K),
  sprintf("comparisons: %s",
          paste(vapply(cells, function(c) c$comp, character(1)), collapse = ", "))
), MANIFEST)

cat("\nDone. Outputs under", RF_DIR, "\n")
cat("Run manifest:", MANIFEST, "\n")
