#!/usr/bin/env Rscript
#
# Three DAR tables with BH-FDR, all derived from the SAME unfiltered MAST fit so
# that p-values are directly comparable across tables (nested subsets).
#
#   A  filtered    det_frac > 0.10
#   B  compromise  det_frac > 0.01
#   C  unfiltered  every tested DMR
#
# Two FDR columns:
#   padjust_BH      BH within that table's own test family (what you would report
#                   if this filter were the analysis)
#   padjust_BH_all  BH over the full 74k-feature family (common scale, so the
#                   effect of the filter can be compared directly)

suppressPackageStartupMessages(library(data.table))

DAR_DIR <- "/data2st1/junyi/output/atac1112/dar"
OUT_DIR <- file.path(DAR_DIR, "fdr_tables")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

det <- fread(file.path(OUT_DIR, "region_l1_detection.csv"))
unf <- fread(file.path(DAR_DIR, "dmr_region_l1_cums_unfiltered/MASTNG_unfiltered_all_comparisons.csv"))

C <- merge(unf[tested == TRUE], det, by = c("region_l1", "primerid"), all.x = TRUE)
C[, padjust_BH_all := p.adjust(`Pr(>Chisq)`, "BH"), by = region_l1]

build <- function(dt, cut) {
  if (!is.null(cut)) dt <- dt[!is.na(det_frac) & det_frac > cut]
  dt <- copy(dt)
  dt[, padjust_BH := p.adjust(`Pr(>Chisq)`, "BH"), by = region_l1]
  setorder(dt, region_l1, `Pr(>Chisq)`)
  dt
}

A <- build(C, 0.10)
B <- build(C, 0.01)
C <- build(C, NULL)

fwrite(A, file.path(OUT_DIR, "tableA_filtered_freq0.10_BH.csv"))
fwrite(B, file.path(OUT_DIR, "tableB_freq0.01_BH.csv"))
fwrite(C, file.path(OUT_DIR, "tableC_unfiltered_BH.csv"))

summ <- function(dt, label) {
  dt[, .(table = label, n_tested = .N,
         sig_BH = sum(padjust_BH < 0.05),
         sig_BH_all = sum(padjust_BH_all < 0.05),
         up_BH = sum(padjust_BH < 0.05 & !is.na(coef) & coef > 0),
         down_BH = sum(padjust_BH < 0.05 & !is.na(coef) & coef < 0)), by = region_l1]
}
S <- rbindlist(list(summ(A, "A_freq0.10"), summ(B, "B_freq0.01"), summ(C, "C_unfiltered")))
setorder(S, region_l1, table)
fwrite(S, file.path(OUT_DIR, "summary_BH.csv"))

cat("=== totals ===\n")
print(S[, .(n_tested = sum(n_tested), sig_BH = sum(sig_BH), sig_BH_all = sum(sig_BH_all),
            up = sum(up_BH), down = sum(down_BH)), by = table])

cat("\n=== significant per group (BH within table / BH on full family) ===\n")
w <- dcast(S, region_l1 ~ table, value.var = c("sig_BH"), fill = 0)
print(w[order(-C_unfiltered)])
cat("\n")
w2 <- dcast(S, region_l1 ~ table, value.var = c("sig_BH_all"), fill = 0)
print(w2[order(-C_unfiltered)])
