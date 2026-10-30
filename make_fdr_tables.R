#!/usr/bin/env Rscript
#
# Rebuild the three DAR tables with BH-FDR instead of Bonferroni, and attach
# the per-group detection rate of every DMR so the tables can be re-filtered
# without re-running MAST.
#
#   A  filtered   : freq > 0.10  (the original run)          -> BH
#   B  compromise : freq > 0.01  (subset of the unfiltered)  -> BH
#   C  unfiltered : every tested DMR                         -> BH
#
# Detection rate = fraction of cells of the group (all cells, not the 10k
# subsample) with count > 0; this is the same quantity MAST's freq() uses.

suppressPackageStartupMessages({
  library(Seurat); library(schard); library(Matrix); library(data.table)
})

DAR_DIR  <- "/data2st1/junyi/output/atac1112/dar"
IN_H5AD  <- "/data2st1/junyi/output/atac1112/3REGIONS_dmr_filtered.h5ads"
OUT_DIR  <- file.path(DAR_DIR, "fdr_tables")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------------- detection --
cat("loading h5ad ...\n")
seo <- schard::h5ad2seurat(IN_H5AD)
md  <- seo@meta.data

# composite grouping, identical to the DAR script
l1  <- as.character(md$celltype.L1)
rg  <- as.character(md$Region)
nt  <- as.character(md$celltype.L1_ct)
ntc <- as.character(md$Neurotransmitter_celltype)
fb  <- ifelse(ntc == "GABAergic", "GABA", ifelse(ntc == "Glutamatergic", "Glut", NA_character_))
isn <- !is.na(l1) & l1 == "Neuron"
bad <- isn & !(nt %in% c("GABA", "Glut"))
nt[bad] <- fb[bad]
key <- ifelse(isn, nt, l1)
grp <- ifelse(is.na(key) | is.na(rg), NA_character_, paste0(rg, "_", key))

cat("computing detection rates ...\n")
X <- GetAssayData(seo, slot = "counts")
B <- X > 0                                    # logical sparse, computed once
det_list <- list()
for (g in sort(unique(na.omit(grp)))) {
  cells <- colnames(seo)[!is.na(grp) & grp == g]
  d <- Matrix::rowSums(B[, cells, drop = FALSE])
  det_list[[g]] <- data.table(region_l1 = g, primerid = rownames(X),
                              n_cells = length(cells), n_detected = d,
                              det_frac = d / length(cells))
  cat(sprintf("  %-18s %6d cells\n", g, length(cells)))
}
det <- rbindlist(det_list)
fwrite(det, file.path(OUT_DIR, "region_l1_detection.csv"))
cat("detection table:", nrow(det), "rows\n\n")

## ------------------------------------------------------------------ tables --
filt <- fread(file.path(DAR_DIR, "dmr_region_l1_cums/MASTNG_regionL1_cums_all_comparisons.csv"))
unf  <- fread(file.path(DAR_DIR, "dmr_region_l1_cums_unfiltered/MASTNG_unfiltered_all_comparisons.csv"))

add_det <- function(dt) merge(dt, det, by = c("region_l1", "primerid"), all.x = TRUE, sort = FALSE)
add_fdr <- function(dt) {
  dt <- dt[!is.na(`Pr(>Chisq)`)]
  dt[, padjust_BH := p.adjust(`Pr(>Chisq)`, "BH"), by = region_l1]
  dt[order(region_l1, `Pr(>Chisq)`)]
}

# A: original freq > 0.10
A <- add_fdr(add_det(filt))
fwrite(A, file.path(OUT_DIR, "tableA_filtered_freq0.10_BH.csv"))

# B + C: from the unfiltered run
U <- add_det(unf)
C <- add_fdr(U[tested == TRUE])
fwrite(C, file.path(OUT_DIR, "tableC_unfiltered_BH.csv"))

B <- add_fdr(U[tested == TRUE & !is.na(det_frac) & det_frac > 0.01])
fwrite(B, file.path(OUT_DIR, "tableB_freq0.01_BH.csv"))

## ---------------------------------------------------------------- summary --
summ <- function(dt, label) {
  dt[, .(table = label,
         n_tested = .N,
         sig_BH_0.05 = sum(padjust_BH < 0.05, na.rm = TRUE),
         sig_BH_0.10 = sum(padjust_BH < 0.10, na.rm = TRUE),
         up_BH  = sum(padjust_BH < 0.05 & coef > 0, na.rm = TRUE),
         down_BH = sum(padjust_BH < 0.05 & coef < 0, na.rm = TRUE)), by = region_l1]
}
S <- rbindlist(list(summ(A, "A_freq0.10"), summ(B, "B_freq0.01"), summ(C, "C_unfiltered")))
S[, sig_BH_0.05 := as.integer(sig_BH_0.05)]
setorder(S, region_l1, table)
fwrite(S, file.path(OUT_DIR, "summary_BH.csv"))

cat("\n=== per-table totals ===\n")
print(S[, .(n_tested = sum(n_tested), sig_BH_0.05 = sum(sig_BH_0.05)), by = table])

cat("\n=== significant DMRs (BH < 0.05) per group ===\n")
w <- dcast(S, region_l1 ~ table, value.var = "sig_BH_0.05", fill = 0)
print(w[order(-`C_unfiltered`)])

cat("\ndone ->", OUT_DIR, "\n")
