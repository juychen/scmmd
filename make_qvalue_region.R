#!/usr/bin/env Rscript
#
# Region-level DAR (20k cells, freq > 0.01): add FDR columns.
#
#   padjust_BH : Benjamini-Hochberg (what was used before)
#   qvalue     : Storey q-value -- the "FDR" estimate requested here
#
# Both control FDR; Storey usually reports more hits because it estimates the
# proportion of true nulls (pi0) instead of assuming pi0 = 1 like BH does.

suppressPackageStartupMessages({
  library(data.table)
  library(qvalue)
})

OUT_DIR <- "/data2st1/junyi/output/atac1112/dar/dmr_region_freq0.01_20k"
in_csv  <- file.path(OUT_DIR, "MASTNG_region20k_f01_all_comparisons.csv")

d <- fread(in_csv)
d[, padjust_BH := p.adjust(`Pr(>Chisq)`, "BH"), by = Region]

d[, qvalue := {
  p <- `Pr(>Chisq)`
  out <- rep(NA_real_, length(p))
  ok  <- !is.na(p)
  if (sum(ok) > 100) {
    q <- tryCatch(qvalue(p[ok])$qvalue, error = function(e) NULL)
    if (is.null(q)) out[ok] <- p.adjust(p[ok], "BH") else out[ok] <- q
  } else {
    out[ok] <- p.adjust(p[ok], "BH")
  }
  out
}, by = Region]

setorder(d, Region, `Pr(>Chisq)`)
fwrite(d, file.path(OUT_DIR, "MASTNG_region20k_f01_FDR.csv"))

s <- d[, .(n_tested = .N,
           sig_q_0.05  = sum(qvalue < 0.05, na.rm = TRUE),
           sig_q_0.10  = sum(qvalue < 0.10, na.rm = TRUE),
           up_q        = sum(qvalue < 0.05 & coef > 0, na.rm = TRUE),
           down_q      = sum(qvalue < 0.05 & coef < 0, na.rm = TRUE),
           sig_BH_0.05 = sum(padjust_BH < 0.05, na.rm = TRUE),
           pi0         = 1 - sum(qvalue < 0.05, na.rm = TRUE) / .N),
       by = Region]
fwrite(s, file.path(OUT_DIR, "summary_FDR_region.csv"))

cat("=== region-level DAR, 20k cells, freq > 0.01 ===\n")
print(s[, .(Region, n_tested, sig_q_0.05, up_q, down_q, sig_BH_0.05)])
cat("\nmin p per region:\n")
print(d[, .(min_p = min(`Pr(>Chisq)`, na.rm = TRUE)), by = Region])
