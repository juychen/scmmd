#!/usr/bin/env Rscript
#
# DAR on DMR regions -- MAST + ngeneson (NG) covariate, arbitrary grouping/contrast.
# Adapted from 11darmast_new_dmr.R (the 2026-07-15 celltype.L2 run).
#
# !!! Design notes -- do not "simplify" these away !!!
#   1. The h5ad count matrix X is the DMR SUBSET, but normalisation uses
#      `total_counts` from obs, which was computed on the FULL peak matrix by
#      10mkdmrmatrix.ipynb  ->  CPM against the real library depth.
#   2. `ngeneson` is read from obs (z-scored `npeaks_on` on the FULL peak
#      matrix). It is NOT recomputed from the DMR subset; recomputation only
#      happens as a fallback when the column is missing / all NA.
#   3. Model: zlm(~ compare_group + ngeneson), optionally + batch.
#      logFC = group2 - group1 (default: CUMS - CONTROL).
#
# Usage:
#   Rscript 11darmast_dmr_cums.R --input 3REGIONS_dmr_filtered.h5ads \
#       --output <outdir> --group_by Region --group1 MW --group2 MC

suppressPackageStartupMessages({
  suppressWarnings({
    library(Seurat)
    library(MAST)
    library(dplyr)
    library(data.table)
    library(argparse)
    library(schard)
  })
})


parser <- ArgumentParser(description = "DAR on DMR regions with MAST + ngeneson covariate")

parser$add_argument("--input", type = "character", required = TRUE,
                    help = "DMR h5ad (e.g. 3REGIONS_dmr_filtered.h5ads)")
parser$add_argument("--output", type = "character", required = TRUE,
                    help = "output directory (created if missing)")
parser$add_argument("--group_by", type = "character", default = "celltype.L2",
                    help = "obs column to loop over, e.g. Region / celltype.L1 / celltype.L2")
parser$add_argument("--compare_by", type = "character", default = "Condition",
                    help = "obs column holding the two conditions (Condition == expriment here)")
parser$add_argument("--group1", type = "character", default = "MW",
                    help = "reference / CONTROL level (MW = CON in this dataset)")
parser$add_argument("--group2", type = "character", default = "MC",
                    help = "case level (MC = CUMS/SUS); logFC = group2 - group1")
parser$add_argument("--batch_by", type = "character", default = "date",
                    help = "batch column, comes from --qc_csv merge")
parser$add_argument("--use_batch", type = "character", default = "no",
                    choices = c("no", "yes", "auto"),
                    help = "no = reproduce the 2026-07-15 model; auto = keep batch only when every batch has both conditions")
parser$add_argument("--exclude", type = "character", default = "",
                    help = "comma separated levels of --group_by to drop, e.g. 'OPC' (default: none)")
parser$add_argument("--prefix", type = "character", default = "",
                    help = "prefix of output file names")
parser$add_argument("--qc_csv", type = "character", default = "/data2st1/junyi/output/atac0627/frac_qc.csv",
                    help = "per-sample QC table (sample, frac_peak, date); skipped if missing")
parser$add_argument("--max_cells", type = "integer", default = 10000,
                    help = "cells kept per group when subsampling (0 = never subsample)")
parser$add_argument("--subsample_threshold", type = "integer", default = 25000,
                    help = "subsample only when a group has more cells than this")
parser$add_argument("--min_cells_peak", type = "integer", default = 10,
                    help = "a DMR must be detected in at least this many cells of the group")
parser$add_argument("--freq_expressed", type = "double", default = 0.1,
                    help = "MAST freq() filter on the compared cells (0 = keep every detected DMR)")
parser$add_argument("--chunk_size", type = "integer", default = 4000,
                    help = "DMRs per zlm fit; smaller = less memory, same results")
parser$add_argument("--ncores", type = "integer", default = 1,
                    help = "parallel workers for the chunked fits (Linux fork)")
parser$add_argument("--keep_all", action = "store_true", default = FALSE,
                    help = "write every DMR of the group, untested ones as NA")
parser$add_argument("--scale_factor", type = "double", default = 100000,
                    help = "CPM scale factor (1e5, same as the original run)")

args <- parser$parse_args()

input      <- args$input
output     <- args$output
group_by   <- args$group_by
compare_by <- args$compare_by
group1     <- args$group1
group2     <- args$group2
batch_by   <- args$batch_by
use_batch  <- args$use_batch
prefix     <- args$prefix
qc_csv     <- args$qc_csv

exclude_vals <- if (nzchar(args$exclude)) {
  trimws(strsplit(args$exclude, ",")[[1]])
} else {
  character(0)
}

dir.create(output, showWarnings = FALSE, recursive = TRUE)
setwd(output)

cat("=================================================\n")
cat("input      :", input, "\n")
cat("output     :", output, "\n")
cat("group_by   :", group_by, "\n")
cat("compare_by :", compare_by, " (", group1, "vs", group2, ")\n")
cat("use_batch  :", use_batch, "on", batch_by, "\n")
cat("=================================================\n")


## ---------------------------------------------------------------- loading ----
cat("Reading h5ad ...\n")
seo <- schard::h5ad2seurat(input)
cat("Loaded:", nrow(seo), "features x", ncol(seo), "cells\n")

for (col in setdiff(c(group_by, compare_by), "region_l1")) {
  if (!col %in% colnames(seo@meta.data)) {
    stop(sprintf("column '%s' not found in obs; available: %s",
                 col, paste(colnames(seo@meta.data), collapse = ", ")))
  }
}


## ------------------------------------------------------ normalisation ----
# NOTE: CPM + log1p against the FULL peak matrix depth (obs$total_counts) is a
# purely per-cell operation, so it is applied per group inside the loop below
# (identical values, but avoids materialising a dense 74k x 152k matrix).
if (!"total_counts" %in% colnames(seo@meta.data)) {
  cat("!! total_counts not in meta.data -- will fall back to LogNormalize per group\n")
}


## ------------------------------------------------------------ sample QC ----
if (file.exists(qc_csv) && "sample" %in% colnames(seo@meta.data)) {
  cat("Merging per-sample QC from", qc_csv, "\n")
  df_qc <- read.csv(qc_csv, row.names = 1)
  df_qc$sample <- rownames(df_qc)

  metadata <- seo@meta.data
  metadata$cell_bc <- rownames(metadata)
  meta2 <- merge(metadata, df_qc, by = "sample", all.x = TRUE)
  rownames(meta2) <- meta2$cell_bc
  seo@meta.data <- meta2

  frac_col <- grep("^Fraction.of.high.quality", colnames(seo@meta.data), value = TRUE)
  if (length(frac_col) == 1) {
    colnames(seo@meta.data)[colnames(seo@meta.data) == frac_col] <- "frac_peak"
  }
  if (!batch_by %in% colnames(seo@meta.data)) {
    cat("!! batch column", batch_by, "not available after QC merge\n")
  }
} else {
  cat("!! QC csv not found or no 'sample' column -- skipping QC merge\n")
}


## ------------------------------------------- optional composite grouping ----
# `--group_by region_l1` builds Region x (L1, neurons split into GABA / Glut):
#   AMY_GABA, AMY_Glut, AMY_Astro-Epen, AMY_Immune, AMY_OPC-Oligo, AMY_Vascular, ...
if (identical(group_by, "region_l1")) {
  need <- c("Region", "celltype.L1", "celltype.L1_ct", "Neurotransmitter_celltype")
  miss <- setdiff(need, colnames(seo@meta.data))
  if (length(miss) > 0) {
    stop("group_by=region_l1 requires obs columns: ", paste(miss, collapse = ", "))
  }

  l1  <- as.character(seo@meta.data$celltype.L1)
  rg  <- as.character(seo@meta.data$Region)
  nt  <- as.character(seo@meta.data$celltype.L1_ct)          # GABA / Glut for neurons
  ntc <- as.character(seo@meta.data$Neurotransmitter_celltype)
  fb  <- ifelse(ntc == "GABAergic", "GABA",
                ifelse(ntc == "Glutamatergic", "Glut", NA_character_))

  is_neuron <- !is.na(l1) & l1 == "Neuron"
  bad_nt    <- is_neuron & !(nt %in% c("GABA", "Glut"))
  if (any(bad_nt)) {
    cat("  NT fallback applied to", sum(bad_nt), "Neuron cells\n")
    nt[bad_nt] <- fb[bad_nt]
  }

  key <- ifelse(is_neuron, nt, l1)
  seo@meta.data$region_l1 <- ifelse(is.na(key) | is.na(rg), NA_character_,
                                    paste0(rg, "_", key))
  group_by <- "region_l1"

  cat("Composite grouping 'region_l1' created:", 
      length(unique(na.omit(seo@meta.data$region_l1))), "groups\n")
  tb <- as.data.frame(table(seo@meta.data$region_l1, useNA = "ifany"))
  colnames(tb) <- c("group", "n_cells")
  print(tb[tb$n_cells > 0, ], row.names = FALSE)
}


## ------------------------------------------------------- MAST analysis ----
perform_mast_analysis <- function(seurat_obj,
                                  group.by,
                                  compare.by,
                                  group1,
                                  group2,
                                  batch.by,
                                  use_batch_mode = "no",
                                  freq_expressed = 0.1,
                                  min_cells_peak = 10,
                                  max_cells = 10000,
                                  subsample_threshold = 25000,
                                  scale_factor = 100000,
                                  chunk_size = 4000,
                                  ncores = 1,
                                  keep_all = FALSE,
                                  save_as_tmp = TRUE) {
  library(MAST)
  library(data.table)

  metadata     <- seurat_obj@meta.data
  group_levels <- unique(as.character(metadata[[group.by]]))
  group_levels <- group_levels[!is.na(group_levels)]
  group_levels <- setdiff(group_levels, exclude_vals)

  all_results <- data.frame()

  for (cell_group in group_levels) {
    tryCatch({
      cat("====================================\n")
      cat("Analyzing", group.by, ":", cell_group, "\n")

      cells_in_group <- rownames(metadata)[
        !is.na(metadata[[group.by]]) & as.character(metadata[[group.by]]) == cell_group
      ]
      if (length(cells_in_group) < 2) {
        cat("  Skipped: too few cells\n")
        next
      }
      seo_subset <- subset(seurat_obj, cells = cells_in_group)

      n_before <- length(cells_in_group)
      if (max_cells > 0 && n_before > subsample_threshold) {
        set.seed(123)
        cells_in_group <- sample(cells_in_group, max_cells)
        cat(sprintf("  subsampled %d -> %d cells (seed 123)\n", n_before, max_cells))
        seo_subset <- subset(seurat_obj, cells = cells_in_group)
      }

      # a DMR has to be seen in a minimum number of cells of this group
      exprs <- GetAssayData(seo_subset, slot = "counts")
      peaks_to_keep  <- Matrix::rowSums(exprs > 0) >= min_cells_peak
      seo_subset2    <- seo_subset[peaks_to_keep, ]
      cat(sprintf("  %d / %d DMRs kept (>=%d cells)\n",
                  sum(peaks_to_keep), length(peaks_to_keep), min_cells_peak))

      # CPM + log1p against the FULL peak matrix depth (obs$total_counts).
      # per-cell operation -> same values as doing it on the whole object.
      if ("total_counts" %in% colnames(seo_subset2@meta.data)) {
        tc <- seo_subset2@meta.data$total_counts
        tc[is.na(tc) | tc == 0] <- 1
        counts_sub <- GetAssayData(seo_subset2, slot = "counts")
        # column scaling via a sparse product -- sweep() would densify the matrix
        D <- Matrix::Diagonal(x = scale_factor / tc)
        norm_mat <- log1p(counts_sub %*% D)
        dimnames(norm_mat) <- dimnames(counts_sub)
        seo_subset2[["RNA"]]@data <- norm_mat
      } else {
        cat("  !! total_counts missing -- falling back to LogNormalize\n")
        seo_subset2 <- NormalizeData(seo_subset2,
                                     normalization.method = "LogNormalize",
                                     scale.factor = scale_factor)
      }

      md        <- seo_subset2@meta.data
      cells_g1  <- rownames(md)[as.character(md[[compare.by]]) == group1]
      cells_g2  <- rownames(md)[as.character(md[[compare.by]]) == group2]
      if (length(cells_g1) == 0 || length(cells_g2) == 0) {
        cat("  Skipped: one of the conditions has zero cells\n")
        next
      }
      cat(sprintf("  cells: %s=%d  %s=%d\n", group1, length(cells_g1), group2, length(cells_g2)))

      selected_cells <- c(cells_g1, cells_g2)
      anno.tmp <- md[selected_cells, ]

      # batch usability
      use_batch_here <- FALSE
      if (use_batch_mode != "no" && batch.by %in% colnames(anno.tmp)) {
        bt <- table(anno.tmp[[batch.by]], anno.tmp[[compare.by]])
        cat("  batch x condition:\n"); print(bt)
        both_conditions <- all(rowSums(bt > 0) == ncol(bt))
        if (use_batch_mode == "yes") {
          use_batch_here <- TRUE
        } else if (use_batch_mode == "auto") {
          use_batch_here <- nrow(bt) > 1 && both_conditions
        }
        cat("  batch in model:", use_batch_here, "\n")
      }

      contrast_name <- paste0("compare_group", group2)

      # ---- one MAST fit on a chunk of DMRs --------------------------------
      fit_one <- function(idx) {
        dat <- as.matrix(expr_matrix[idx, selected_cells, drop = FALSE])
        sca <- suppressMessages(FromMatrix(dat, cData = anno.tmp))

        sg <- freq(sca) > freq_expressed
        if (!any(sg)) return(NULL)
        sca <- sca[sg, ]

        colData(sca)$compare_group <- factor(colData(sca)[[compare.by]],
                                             levels = c(group1, group2))

        # NG: reuse what is in obs (full peak matrix); compute only if absent
        if (!"ngeneson" %in% colnames(colData(sca)) ||
            all(is.na(colData(sca)$ngeneson))) {
          colData(sca)$ngeneson <- scale(colSums(assay(sca) > 0))
        } else {
          colData(sca)$ngeneson <- as.numeric(colData(sca)$ngeneson)
        }

        if (use_batch_here) {
          colData(sca)$batch <- factor(colData(sca)[[batch.by]])
          fit <- zlm(~ compare_group + ngeneson + batch, sca)
        } else {
          fit <- zlm(~ compare_group + ngeneson, sca)
        }

        s <- summary(fit, doLRT = contrast_name)$datatable
        fc <- merge(
          s[contrast == contrast_name & component == 'H',
            .(primerid, `Pr(>Chisq)`)],
          s[contrast == contrast_name & component == 'logFC',
            .(primerid, coef, ci.hi, ci.lo)],
          by = 'primerid'
        )
        if (nrow(fc) == 0) return(NULL)
        as.data.frame(fc)
      }

      # pre-compute the MAST freq() criterion on the sparse counts (identical to
      # freq(sca) > freq_expressed) so chunks never contain untestable DMRs
      counts_sel <- GetAssayData(seo_subset2, slot = "counts")[,
                        selected_cells, drop = FALSE]
      det_frac   <- Matrix::rowSums(counts_sel > 0) / length(selected_cells)
      keep_feat  <- det_frac > freq_expressed
      if (!any(keep_feat)) {
        cat("  Warning: no DMR passed the freq filter\n")
        next
      }
      expr_matrix <- GetAssayData(seo_subset2, layer = "data")[keep_feat, ,
                                                                drop = FALSE]

      feat_rows <- seq_len(nrow(expr_matrix))
      chunks    <- split(feat_rows, ceiling(seq_along(feat_rows) / chunk_size))
      cat(sprintf("  fitting %d DMRs in %d chunk(s) with %d core(s), freq > %.3f\n",
                  length(feat_rows), length(chunks), ncores, freq_expressed))

      res_list <- if (ncores > 1) {
        parallel::mclapply(chunks, fit_one, mc.cores = ncores,
                           mc.preschedule = FALSE, mc.silent = TRUE)
      } else {
        lapply(chunks, fit_one)
      }
      res_list <- Filter(Negate(is.null), res_list)

      if (length(res_list) == 0) {
        cat("  Warning: no DMR passed filtering\n")
        next
      }
      fcHurdle <- as.data.frame(data.table::rbindlist(res_list, fill = TRUE))
      cat(sprintf("  %d DMRs tested after MAST freq > %.3f\n",
                  nrow(fcHurdle), freq_expressed))

      fcHurdle$padjust         <- p.adjust(fcHurdle[["Pr(>Chisq)"]], "bonferroni")
      fcHurdle[[group.by]]     <- cell_group
      fcHurdle$n_cells_group1  <- length(cells_g1)
      fcHurdle$n_cells_group2  <- length(cells_g2)

      # keep every DMR of this group in the table: untested ones get NA
      if (keep_all) {
        all_ids <- rownames(seo_subset)
        mm <- match(all_ids, fcHurdle$primerid)
        out <- fcHurdle[mm, , drop = FALSE]
        out$primerid <- all_ids
        out[[group.by]] <- cell_group
        out$n_cells_group1 <- length(cells_g1)
        out$n_cells_group2 <- length(cells_g2)
        out$tested <- !is.na(out[["Pr(>Chisq)"]])
        fcHurdle <- out
      }

      all_results <- rbind(all_results, fcHurdle)

      if (save_as_tmp) {
        safe_prefix <- if (nzchar(prefix)) paste0(prefix, "_") else ""
        safe_group  <- gsub(" ", "_", gsub("/", "-", cell_group, fixed = TRUE))
        fname       <- paste0("tmp.degene.result.", safe_prefix, safe_group)
        save(fcHurdle, file = paste0(fname, ".rda"))
        write.csv(fcHurdle, file = paste0(fname, ".csv"), row.names = FALSE)
      }
    }, error = function(e) {
      cat("  Error in", cell_group, ":", e$message, "\n")
    })
  }

  if (nrow(all_results) == 0) {
    warning("No differential results were generated.")
  }
  all_results
}


res <- perform_mast_analysis(
  seo,
  group.by             = group_by,
  compare.by           = compare_by,
  group1               = group1,
  group2               = group2,
  batch.by             = batch_by,
  use_batch_mode       = use_batch,
  freq_expressed       = args$freq_expressed,
  min_cells_peak       = args$min_cells_peak,
  max_cells            = args$max_cells,
  subsample_threshold  = args$subsample_threshold,
  scale_factor         = args$scale_factor,
  chunk_size           = args$chunk_size,
  ncores               = args$ncores,
  keep_all             = args$keep_all
)


## ---------------------------------------------------------------- output ----
if (nrow(res) > 0) {
  safe_prefix <- if (nzchar(prefix)) paste0(prefix, "_") else ""
  out_csv <- file.path(output, paste0("MASTNG_", safe_prefix, "all_comparisons.csv"))
  write.csv(res, out_csv, row.names = FALSE)
  cat("\nSaved combined table:", out_csv, "\n")
  cat("rows:", nrow(res), "\n")

  sig <- as.data.table(res)[padjust < 0.05]
  cat("significant (bonferroni padjust < 0.05):", nrow(sig), "\n")
  if (nrow(sig) > 0) {
    summ <- sig[, .(n_sig = .N,
                    n_up   = sum(coef > 0),
                    n_down = sum(coef < 0)),
                by = group_by]
    print(summ)
    write.csv(summ, file.path(output, paste0("MASTNG_", safe_prefix, "sig_summary.csv")),
              row.names = FALSE)
  }
} else {
  cat("\nNo results to write.\n")
}

cat("\nDone.\n")
