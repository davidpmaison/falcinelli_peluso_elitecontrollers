# 02_pathway_enrichment.R — (GO BP/CC/MF, KEGG)

suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi); library(dplyr)
})

# =============================================================================
# CONFIG — edit these
# =============================================================================
DATA_PATH       <- "../data/your_npx_data.csv"
ASSAY_COL       <- "Assay"
SIG_THRESHOLD   <- 0.05      # p.adjust cutoff for "significant" terms
GO_SIMPLIFY_CUT <- 0.7       # semantic similarity cutoff for clusterProfiler::simplify()

COMPARISONS <- list(
  list(ttest = "GROUPA_vs_GROUPB/GROUPA_vs_GROUPB_ttest_results.csv",
       out_prefix = "GROUPA_vs_GROUPB/GROUPA_vs_GROUPB_OUTPUT_results_")
)

# =============================================================================
# Background mapping
# =============================================================================
all_assays <- unique(read.csv(DATA_PATH, stringsAsFactors = FALSE)[[ASSAY_COL]])
sym2ent <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = all_assays,
                                 column = "ENTREZID", keytype = "SYMBOL",
                                 multiVals = "first")
sym2ent <- sym2ent[!is.na(sym2ent)]
bg_entrez <- unique(unname(sym2ent))
cat(sprintf("Background: %d / %d assays mapped to Entrez\n",
            length(bg_entrez), length(all_assays)))

# =============================================================================
# Per-comparison enrichment
# =============================================================================
run_enrichment <- function(ttest_path, out_prefix) {
  d <- read.csv(ttest_path, stringsAsFactors = FALSE)
  sig_syms <- d$Assay[d$Adjusted_pval < SIG_THRESHOLD]
  sig_ent  <- unique(unname(sym2ent[sig_syms])); sig_ent <- sig_ent[!is.na(sig_ent)]
  cat(sprintf("\n[%s]  sig=%d, mapped=%d\n", basename(ttest_path), length(sig_syms), length(sig_ent)))

  go_ora <- function(ont) {
    res <- enrichGO(gene = sig_ent, universe = bg_entrez, OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID", ont = ont,
                    pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", readable = TRUE)
    if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
    res_simp <- clusterProfiler::simplify(res, cutoff = GO_SIMPLIFY_CUT, by = "p.adjust", select_fun = min)
    d2 <- as.data.frame(res_simp); d2[d2$p.adjust < SIG_THRESHOLD, , drop = FALSE]
  }
  kegg_ora <- function() {
    res <- enrichKEGG(gene = sig_ent, universe = bg_entrez, organism = "hsa",
                      pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH")
    if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
    d2 <- as.data.frame(res); d2[d2$p.adjust < SIG_THRESHOLD, , drop = FALSE]
  }

  bp <- go_ora("BP"); cc <- go_ora("CC"); mf <- go_ora("MF"); kg <- kegg_ora()
  if (!is.null(bp)) write.csv(bp, paste0(out_prefix, "GO_BP_simplified_sig.csv"), row.names = FALSE)
  if (!is.null(cc)) write.csv(cc, paste0(out_prefix, "GO_CC_simplified_sig.csv"), row.names = FALSE)
  if (!is.null(mf)) write.csv(mf, paste0(out_prefix, "GO_MF_simplified_sig.csv"), row.names = FALSE)
  if (!is.null(kg)) write.csv(kg, paste0(out_prefix, "KEGG_sig.csv"),             row.names = FALSE)
  cat(sprintf("  BP=%d CC=%d MF=%d KEGG=%d\n",
              ifelse(is.null(bp),0,nrow(bp)), ifelse(is.null(cc),0,nrow(cc)),
              ifelse(is.null(mf),0,nrow(mf)), ifelse(is.null(kg),0,nrow(kg))))
}

for (cmp in COMPARISONS) run_enrichment(cmp$ttest, cmp$out_prefix)
