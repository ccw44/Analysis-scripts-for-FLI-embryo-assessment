## DSS: DMLtest + callDMR on the locked comparable CpG set.
## Direction: group1 = FAR, group2 = Control, so DML$diff > 0 means FAR higher.

suppressPackageStartupMessages({
  library(data.table); library(bsseq); library(DSS); library(GenomicRanges)
  library(GenomeInfoDb); library(IRanges)
})
## Data locations - edit this block only; run from the repository root.
##   ref_dir  mm39 annotation (mm39.chrom.sizes)
##   out_dir  output root; tables/ and figures/
ref_dir <- "data/mm39"
out_dir <- "results"
TAB     <- file.path(out_dir, "tables")
if (!dir.exists(TAB)) dir.create(TAB, recursive = TRUE)


sizes <- fread(file.path(ref_dir, "mm39.chrom.sizes"), header = FALSE,
               col.names = c("chrom", "size"))
sizes[chrom == "chrM", size := size]
seq_order <- c(paste0("chr", 1:19), "chrX", "chrY", "chrM")
si <- Seqinfo(seqnames = seq_order,
              seqlengths = sizes$size[match(seq_order, sizes$chrom)],
              genome = "mm39")

RUN <- list(
  smoothing      = FALSE,   # WGBS setting; 500-bp moving-average window (historical parameter)
  smoothing.span = 500,
  equal.disp     = FALSE,
  delta          = 0,       # callDMR default (not part of the locked parameter list)
  p.threshold    = 1e-5,
  minlen         = 50,
  minCG          = 3,
  dis.merge      = 100,
  pct.sig        = 0.50,
  ncores         = 4        # resource cap: <= 8 cores per task on this server
)
print(RUN)

for (tissue in c("liver", "blood")) {
  cat("\n================ ", tissue, " ================\n")
  t0 <- Sys.time()
  dt <- fread(file.path(TAB, paste0("comparable_", tissue, ".tsv.gz")))
  ctrl <- grep("^Control-", grep("_M$", names(dt), value = TRUE), value = TRUE)
  far  <- grep("^FAR-",     grep("_M$", names(dt), value = TRUE), value = TRUE)
  ctrl <- sub("_M$", "", ctrl); far <- sub("_M$", "", far)
  sids <- c(far, ctrl)
  cat("samples:", paste(sids, collapse = ", "), "\n")
  cat("comparable CpGs in table:", nrow(dt), "\n")

  M   <- as.matrix(dt[, paste0(sids, "_M"),   with = FALSE])
  Cov <- as.matrix(dt[, paste0(sids, "_cov"), with = FALSE])
  colnames(M) <- sids; colnames(Cov) <- sids   # BSseq requires assay colnames == sampleNames
  gr  <- GRanges(seqnames = dt$chrom, ranges = IRanges(start = dt$start + 1L, end = dt$end))
  seqinfo(gr) <- si
  BSobj <- BSseq(gr = gr, M = M, Cov = Cov, sampleNames = sids)
  BSobj <- sort(BSobj)

  cov_mat <- getCoverage(BSobj, type = "Cov")
  keep <- rowSums(cov_mat > 0) >= 2          # DSS needs >=2 samples with coverage at a locus
  cat("loci with coverage in >=2 samples:", sum(keep), "of", nrow(cov_mat),
      "-> dropped", sum(!keep), "\n")
  BSobj <- BSobj[keep, ]

  dml <- DMLtest(BSobj, group1 = far, group2 = ctrl,
                 smoothing = RUN$smoothing, smoothing.span = RUN$smoothing.span,
                 equal.disp = RUN$equal.disp, ncores = RUN$ncores)
  dt_dml <- as.data.table(dml)
  dt_dml[, direction := fifelse(diff > 0, "FAR_high", fifelse(diff < 0, "FAR_low", "none"))]
  fwrite(dt_dml, file.path(TAB, paste0("DML_all_", tissue, ".tsv.gz")), sep = "\t")
  cat("DML rows written:", nrow(dt_dml),
      "| FDR<0.05:", sum(dt_dml$fdr < 0.05, na.rm = TRUE),
      "| p<1e-5:", sum(dt_dml$pval < 1e-5, na.rm = TRUE), "\n")

  dmr <- callDMR(dml, delta = RUN$delta, p.threshold = RUN$p.threshold,
                 minlen = RUN$minlen, minCG = RUN$minCG, dis.merge = RUN$dis.merge,
                 pct.sig = RUN$pct.sig)
  dmr <- as.data.table(dmr)
  if (nrow(dmr) > 0) {
    dmr[, direction := fifelse(diff.Methy > 0, "FAR_high", fifelse(diff.Methy < 0, "FAR_low", "none"))]
    dmr[, tissue := tissue]
  } else {
    dmr <- data.table(chr = character(0), start = integer(0), end = integer(0),
                      length = integer(0), nCG = integer(0), meanMethy1 = numeric(0),
                      meanMethy2 = numeric(0), diff.Methy = numeric(0), areaStat = numeric(0),
                      direction = character(0), tissue = character(0))
  }
  fwrite(dmr, file.path(TAB, paste0("DSS_regions_", tissue, ".tsv")), sep = "\t")
  cat("DSS regions:", nrow(dmr), "\n")

  ## coverage-aware statistics for the record
  qc <- data.table(
    tissue = tissue,
    n_comparable_cpg = nrow(dt),
    n_cpg_used_by_DSS = sum(keep),
    n_cpg_dropped_cov_lt2_samples = sum(!keep),
    n_dml = nrow(dt_dml),
    n_dml_fdr05 = sum(dt_dml$fdr < 0.05, na.rm = TRUE),
    n_dml_p1e5 = sum(dt_dml$pval < 1e-5, na.rm = TRUE),
    n_dss_regions = nrow(dmr),
    n_dss_regions_FARhigh = sum(dmr$diff.Methy > 0),
    n_dss_regions_FARlow = sum(dmr$diff.Methy < 0),
    dss_params = paste(names(RUN), unlist(RUN), sep = "=", collapse = ";"),
    elapsed_min = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
  )
  fwrite(qc, file.path(TAB, paste0("dss_run_qc_", tissue, ".csv")))
  print(qc)
  rm(dt, M, Cov, gr, BSobj, cov_mat, keep, dml, dt_dml, dmr); gc()
}

cat("\n=== sessionInfo ===\n")
print(sessionInfo())
