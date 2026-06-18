# two_group_analysis.R 

suppressPackageStartupMessages({ library(dplyr) })

# =============================================================================
# CONFIG — edit these
# =============================================================================
DEMO_PATH <- "../data/your_clinical_data.csv"
OUTPUT_CSV <- "results.csv"

# The two groups to compare.
GROUP_COL  <- "Cohort"        # column in DEMO_PATH that holds group labels
GROUP_POS  <- "GroupA"        # "positive" / first group
GROUP_NEG  <- "GroupB"        # "negative" / second group

TESTS <- list(
  c("logIntactVirus",           "mw"),
  c("CD4abs",                   "mw"),
  c("CD4CD8ratio",              "mw"),
  c("Days_since_HIV_diagnosis", "welch")
)

SUBSET_SAMPLE_IDS <- NULL
SUBSET_PARENT_GROUP <- NULL    # only used if SUBSET_SAMPLE_IDS is set

# =============================================================================
# Pipeline
# =============================================================================
d <- read.csv(DEMO_PATH, stringsAsFactors = FALSE)

for (v in vapply(TESTS, `[[`, "", 1)) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))

if (!is.null(SUBSET_SAMPLE_IDS)) {
  d <- d %>% filter(.data[[GROUP_COL]] == SUBSET_PARENT_GROUP)
  d$.group_eval <- ifelse(d$SampleID %in% SUBSET_SAMPLE_IDS, GROUP_POS, GROUP_NEG)
} else {
  d <- d %>% filter(.data[[GROUP_COL]] %in% c(GROUP_POS, GROUP_NEG))
  d$.group_eval <- d[[GROUP_COL]]
}

do_test <- function(varname, test_kind) {
  v1 <- d[[varname]][d$.group_eval == GROUP_POS]; v1 <- v1[!is.na(v1)]
  v2 <- d[[varname]][d$.group_eval == GROUP_NEG]; v2 <- v2[!is.na(v2)]
  if (test_kind == "welch") {
    res <- t.test(v1, v2, var.equal = FALSE); stat <- unname(res$statistic); sn <- "t"
    tt  <- "Welch t-test on raw values"
  } else {
    res <- suppressWarnings(wilcox.test(v1, v2, exact = FALSE)); stat <- unname(res$statistic); sn <- "W"
    tt  <- "Mann-Whitney test on raw values"
  }
  data.frame(
    variable = varname, test_type = tt, term = GROUP_COL,
    statistic = stat, statistic_name = sn, p.value = res$p.value,
    normality_min_p = NA_real_, n_total = length(v1) + length(v2),
    section = "Main Test", sample_size_warning = NA_character_,
    n_group1 = length(v1), n_group2 = length(v2),
    n_outliers = NA_integer_, outlier_strategy = "none",
    decision_reason = "Unadjusted test on raw values",
    note = sprintf("%s vs %s | group1=%d, group2=%d", GROUP_POS, GROUP_NEG, length(v1), length(v2)),
    stringsAsFactors = FALSE
  )
}

rows <- do.call(rbind, lapply(TESTS, function(t) do_test(t[1], t[2])))
write.csv(rows, OUTPUT_CSV, row.names = FALSE, quote = TRUE, na = "")
cat("Wrote", OUTPUT_CSV, "\n")
print(rows %>% select(variable, test_type, p.value, n_total, n_group1, n_group2))
