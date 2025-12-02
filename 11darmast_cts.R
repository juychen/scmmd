suppressPackageStartupMessages({
  suppressWarnings({
    library(Seurat)
    library(MAST)
    library(dplyr)
    library(data.table)
    library(ggplot2)
    library(schard)
    library(argparse)
  })
})


parser <- ArgumentParser(description = "Description of your script")
parser$add_argument("--input", type = "character", required = TRUE, help = "Path to the input file")
parser$add_argument("--output", type = "character", required = TRUE, help = "Path to the output file")
parser$add_argument("--chunk", type = "integer", default = 30000, help = "Number of regions to process in each chunk")
parser$add_argument("--ncores", type = "integer", default = 8, help = "Number of cores to use for parallel processing")
parser$add_argument("--region", type = "character", default = "PFC", help = "Column name for region counts in the input file")

args <- parser$parse_args()
input <- args$input
output <- args$output
n_cores <- args$ncores
region <- args$region


# seo = schard::h5ad2seurat(paste0(folder,'/',brainregion, ".h5ad"))

# seo <- NormalizeData(seo, normalization.method = "LogNormalize", scale.factor = 10000)


# region = brainregion

seo = schard::h5ad2seurat(input)
#seo <- subset(seo, subset = !(sample %in% c("MW26A_PFC", "MC25A_PFC")))
df_qc <- read.csv('/data2st1/junyi/output/atac0627/frac_qc.csv', row.names = 1)
df_qc$sample <- rownames(df_qc)
metadata <- seo@meta.data
cell_groups <- unique(metadata[["celltype.L2"]])
metadata$cell_bc <- rownames(metadata)
meta2 <- merge(metadata, df_qc, by = "sample", all.x = TRUE)
rownames(meta2) <- meta2$cell_bc
seo@meta.data <- meta2
colnames(seo@meta.data)[colnames(seo@meta.data) == "Fraction.of.high.quality.fragments.overlapping.peaks"] <- "frac_peak"
setwd(output)

if (ncol(seo) > 10000) {
  # Randomly sample 10000 cells
  set.seed(123)
  cells_keep <- sample(colnames(seo), 10000)
  seo_subset <- subset(seo, cells = cells_keep)
}



perform_mast_celltype_specific <- function(
    seurat_obj,
    group.by = "celltype.L2",    # cell type column
    batch.by = NULL,
    freq_expressed = 0.1,
    save.as.tmp = TRUE
){
  library(MAST)
  library(data.table)

  metadata <- seurat_obj@meta.data
  celltypes <- unique(metadata[[group.by]])

  all_results <- data.frame()

  for (ct in celltypes) {

    cat("====================================\n")
    cat("Celltype-specific DAR for:", ct, "\n")
    if (region == "") {
        region_suffix <- ""
      } else {
        region_suffix <- paste0(region,"_")
    }
    if (file.exists(paste0("DAR_celltype_specific_", region_suffix, make.names(ct), ".csv"))){
      cat("  Skipping existing result for:", ct, "\n")
      next
    }
    tryCatch({

    # Define labels
    metadata$compare_group <- ifelse(metadata[[group.by]] == ct, ct, "others")
    metadata$compare_group <- factor(metadata$compare_group, levels=c(ct,"others"))
    seurat_obj$compare_group <- metadata$compare_group

    # Subset only relevant cells (to reduce memory)
    # keep_cells <- rownames(metadata)
    # seo_subset <- subset(seurat_obj, cells = keep_cells)

    # Filter peaks expressed in >=10 cells
    exprs <- GetAssayData(seurat_obj, slot = "counts")
    keep <- rowSums(exprs > 0) >= 10
    seo_subset <- seurat_obj[keep, ]

    expr_matrix <- GetAssayData(seo_subset, slot = "data")
    anno <- seo_subset@meta.data

    # Build MAST object
    sca <- FromMatrix(as.matrix(expr_matrix), cData = anno)

    # Frequency filter
    select.genes <- freq(sca) > freq_expressed
    sca <- sca[select.genes, ]

    # Add CDR2 (ngenes)
    cdr2 <- colSums(assay(sca) > 0)
    colData(sca)$ngeneson <- scale(cdr2)

    # Add batch effect if exists
    # if (!is.null(batch.by)) {
    #   colData(sca)$batch <- factor(colData(sca)[[batch.by]])
    #   zlm_model <- zlm(~ compare_group + ngeneson + batch, sca)
    # } else {
    zlm_model <- zlm(~ compare_group + ngeneson, sca)
    # }

    # MAST contrast test
    contrast <- paste0("compare_groupothers")  
    summary_result <- summary(zlm_model, doLRT = contrast)
    dt <- summary_result$datatable

    fcHurdle <- merge(
      dt[contrast == contrast & component == "H", .(primerid, `Pr(>Chisq)`)],
      dt[contrast == contrast & component == "logFC",
         .(primerid, coef, ci.hi, ci.lo)],
      by="primerid"
    )

    if (nrow(fcHurdle)==0) {
      cat("No results for:", ct, "\n")
      next
    }

    fcHurdle[, padjust := p.adjust(`Pr(>Chisq)`, "bonferroni")]
    fcHurdle$celltype <- ct

    all_results <- rbind(all_results, fcHurdle)

    if (save.as.tmp) {
      if (region == "") {
        region_suffix <- ""
      } else {
        region_suffix <- paste0(region,"_")
      }
      write.csv(fcHurdle,
                paste0("DAR_celltype_specific_", region_suffix, make.names(ct), ".csv"))
    }
    }, error = function(e){
      cat("  Error in cell type ", ct, " : ", e$message, "\n")
    })
  }

  return(all_results)
}
seo <- subset(seo, subset = celltype.L1_ct != "OPC")
r.dar_celltype <- perform_mast_celltype_specific(
    seo,
    group.by = "celltype.L2",
    batch.by = "date"
)


# colnames(r.deg_M) = c("gene","Pr","avg_log2FC","ci.hi","ci.lo","p_val_adj","celltype")
# write.csv(r.deg_M,"MAST_DEG_M.csv")

# filtered_r.deg_M = r.deg_M %>%
#   filter(p_val_adj < 0.05)

# write.csv(filtered_r.deg_M,"MAST_DEG_M_filtered.csv")