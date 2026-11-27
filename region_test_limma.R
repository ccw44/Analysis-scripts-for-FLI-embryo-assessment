## Step 05 - region-level inference: moderated t (limma eBayes) on per-sample region means
## This is the Python-implementable stand-in for dmrseq's region-level inference (dmrseq is not
## installed in this environment): candidate regions are locked from coverage only (step 04),
## every candidate is tested genome-wide (not just DSS-selected regions), and the region-level
## FDR (BH) is computed over the tested region set.
##
## Direction: coefficient 2 = FAR - Control (positive = FAR hypermethylated).
## Pre-specified variance models, both reported:
##   primary     : eBayes(trend = TRUE,  robust = TRUE)   mean-variance trend for beta values
##   sensitivity : eBayes(trend = FALSE, robust = FALSE)  constant-variance model

## Packages are expected in the default library; set R_LIBS otherwise.
suppressPackageStartupMessages({ library(data.table); library(limma) })

## ---- data locations: edit this block only ------------------------------
## Run the script from the repository root; relative paths resolve against it.
##   ref_dir  mm39 annotation (mm39.chrom.sizes)
##   out_dir  output root; intermediate tables in tables/
ref_dir <- "data/mm39"
out_dir <- "results"
TAB     <- file.path(out_dir, "tables")
if (!dir.exists(TAB)) dir.create(TAB, recursive = TRUE)


for (tissue in c("liver", "blood")) {
  cat("\n===== ", tissue, " =====\n")
  summary_rows <- data.table()
  mm <- fread(file.path(TAB, paste0("region_mean_matrix_", tissue, ".tsv.gz")))
  sids <- setdiff(names(mm), "region_id")
  ctrl <- grep("^Control-", sids, value = TRUE)
  far  <- grep("^FAR-", sids, value = TRUE)
  M <- as.matrix(mm[, ..sids]); rownames(M) <- as.character(mm$region_id)
  group <- factor(c(ifelse(sids %in% far, "FAR", "Control")), levels = c("Control", "FAR"))
  design <- model.matrix(~group)
  cat("design:\n"); print(design)
  cat("regions in matrix:", nrow(M),
      "| complete cases (all 6):", sum(complete.cases(M)), "\n")

  fit <- lmFit(M, design)
  for (variant in c("primary", "sensitivity")) {
    eb <- if (variant == "primary") eBayes(fit, trend = TRUE,  robust = TRUE)
          else                      eBayes(fit, trend = FALSE, robust = FALSE)
    tt <- topTable(eb, coef = 2, number = Inf, sort.by = "none")
    rid <- as.integer(rownames(tt))
    fit_id <- as.integer(rownames(eb$coefficients))
    idx <- match(rid, fit_id)
    dt <- data.table(region_id = rid,
                     logFC = tt$logFC, t = tt$t, p_value = tt$P.Value,
                     q_value = tt$adj.P.Val, s2_post = eb$s2.post[idx],
                     df_residual = eb$df.residual[idx], df_prior = eb$df.prior,
                     AveExpr = tt$AveExpr)
    setnames(dt, c("logFC", "t", "p_value", "q_value", "s2_post", "AveExpr"),
             paste0(c("logFC", "t", "p", "q", "s2_post", "AveExpr"), "_", variant))
    setnames(dt, c("df_residual", "df_prior"), paste0(c("df_residual", "df_prior"), "_", variant))
    if (variant == "primary") out <- dt else out <- merge(out, dt, by = "region_id", all = TRUE)
    qv <- dt[[paste0("q_", variant)]]
    dpr <- eb$df.prior
    dpr_txt <- if (length(unique(dpr)) == 1) format(dpr[1], digits = 4)
               else paste0(format(min(dpr), digits = 3), "..", format(max(dpr), digits = 3))
    row <- data.table(tissue = tissue, variant = variant, df_prior = dpr_txt,
                      n_tested = sum(!is.na(qv)),
                      n_q_lt_005 = sum(qv < 0.05, na.rm = TRUE),
                      n_q_lt_010 = sum(qv < 0.10, na.rm = TRUE),
                      n_q_lt_025 = sum(qv < 0.25, na.rm = TRUE),
                      min_q = min(qv, na.rm = TRUE),
                      n_q005_and_absdiff_gt_020 =
                        sum(qv < 0.05 & abs(dt[[paste0("logFC_", variant)]]) > 0.20, na.rm = TRUE))
    summary_rows <- rbind(summary_rows, row)
    cat(sprintf("  %-12s df.prior=%s tested=%d q<0.05=%d q<0.10=%d q<0.25=%d min_q=%.4g |q<0.05 & |logFC|>0.20|=%d\n",
                variant, dpr_txt, row$n_tested, row$n_q_lt_005, row$n_q_lt_010,
                row$n_q_lt_025, row$min_q, row$n_q005_and_absdiff_gt_020))
  }
  fwrite(out, file.path(TAB, paste0("region_tests_", tissue, ".tsv.gz")), sep = "\t")
  cat("wrote region_tests_", tissue, ".tsv.gz\n", sep = "")

  ## consistency check: limma logFC must equal the contract effect size (FAR - Control)
  eff <- fread(file.path(TAB, paste0("regions_all_", tissue, ".tsv.gz")),
               select = c("region_id", "diff", "n_valid_control", "n_valid_far"))
  chk <- merge(out[, .(region_id, logFC_primary)], eff, by = "region_id")
  chk <- chk[is.finite(diff) & is.finite(logFC_primary)]
  cat(sprintf("  effect-size vs limma logFC: n=%d  max|diff|=%g  cor=%.6f\n",
              nrow(chk), max(abs(chk$diff - chk$logFC_primary)), cor(chk$diff, chk$logFC_primary)))
  fwrite(data.table(tissue = tissue, n_checked = nrow(chk),
                    max_abs_diff = max(abs(chk$diff - chk$logFC_primary)),
                    cor = cor(chk$diff, chk$logFC_primary)),
         file.path(TAB, paste0("effectsize_vs_logFC_check_", tissue, ".csv")))
  fwrite(summary_rows, file.path(TAB, paste0("region_test_summary_", tissue, ".csv")))
}
cat("\nsessionInfo:\n"); print(sessionInfo())
