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
seo <- subset(seo, subset = !(sample %in% c("MW26A_PFC", "MC25A_PFC")))
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


perform_mast_analysis <- function(seurat_obj,
                                  group.by, 
                                  compare.by,
                                  group1, 
                                  group2, 
                                  batch.by,
                                  freq_expressed = 0.1,
                                  save.as.tmp = TRUE) {
  library(MAST)
  library(data.table)

  
  #expr_matrix <- GetAssayData(seurat_obj, layer = "data")
      
  #filter_genes <- !grepl("^Rp", rownames(expr_matrix)) & !grepl("^mt-", rownames(expr_matrix))
  #expr_matrix <- expr_matrix[filter_genes, ]
    
  metadata <- seurat_obj@meta.data
  cell_groups <- unique(metadata[[group.by]])
  
  all_results <- data.frame()
  
  for (cell_group in cell_groups) {
    tryCatch({
      cat("====================================\n")
      cat("Analyzing cell group:", cell_group, "\n")
      

      seo_subset = subset(seurat_obj, subset = celltype.L2 == cell_group)
      
      if (ncol(seo_subset) > 25000) {
        # Randomly sample 10000 cells
        set.seed(123)
        cells_keep <- sample(colnames(seo_subset), 10000)
        seo_subset <- subset(seo_subset, cells = cells_keep)
      }

      # Filter genes expressed in at least 5 cells
      exprs <- GetAssayData(seo_subset, slot = "counts")
      peaks_to_keep <- rowSums(exprs > 0) >= 10
      seo_subset2 <- seo_subset[peaks_to_keep, ]
      expr_matrix <- GetAssayData(seo_subset2, layer = "data")

      metadata_subset <- seo_subset2@meta.data
      cells_in_group <- rownames(metadata_subset[metadata_subset[[group.by]] == cell_group, ])
      cells_group1 <- cells_in_group[metadata_subset[cells_in_group, compare.by] == group1]
      cells_group2 <- cells_in_group[metadata_subset[cells_in_group, compare.by] == group2]
      
      if(length(cells_group1) == 0 || length(cells_group2) == 0) {
        cat("  Skipped: one of the treatment groups has zero cells\n")
        next
      }


      selected_cells <- c(cells_group1, cells_group2)
      dat.tmp <- expr_matrix[, selected_cells]
      anno.tmp <- metadata[selected_cells, ]
      
      # 打印batch和treatment分布表
      batch_treatment_table <- table(anno.tmp[[batch.by]], anno.tmp[[compare.by]])
      cat("  Batch x Treatment distribution:\n")
      print(batch_treatment_table)
      
      # 判断每个batch是否包含两个treatment组
      batches_all_have_both_treatments <- all(rowSums(batch_treatment_table > 0) == ncol(batch_treatment_table))
      
      multiple_batches_exist <- nrow(batch_treatment_table) > 1
      
      # 根据分布决定是否用batch效应
      if (multiple_batches_exist && batches_all_have_both_treatments) {
        use_batch_effect <- TRUE
        cat("  Will use batch effect in model\n")
      } else {
        use_batch_effect <- FALSE
        if(!multiple_batches_exist) {
          cat("  Only one batch found - batch effect not used\n")
        } else if (!batches_all_have_both_treatments) {
          cat("  Not all batches contain both treatments - batch effect not used\n")
        }
      }
      use_batch_effect <- FALSE
      # 创建MAST对象
      sca <- FromMatrix(as.matrix(dat.tmp), cData = anno.tmp)
      
      # 基因表达频率过滤
      select.genes <- freq(sca) > freq_expressed
      sca <- sca[select.genes, ]
      
      # 设定compare_group factor
      colData(sca)$compare_group <- factor(colData(sca)[[compare.by]], levels = c(group1, group2))
      
      # add ngenes
      cdr2 <- colSums(assay(sca) > 0)
      colData(sca)$ngeneson <- scale(cdr2)

      if (use_batch_effect) {
        colData(sca)$batch <- factor(colData(sca)[[batch.by]])
        zlm_model <- zlm(~ compare_group + ngeneson + batch, sca)
        #        zlm_model <- zlm(~ compare_group + ngeneson + batch + frac_peak, sca)

      } else {
        zlm_model <- zlm(~ compare_group + ngeneson, sca)
      }
      
      # 差异分析
      contrast_name <- paste0("compare_group", group2)
      summary_result <- summary(zlm_model, doLRT = contrast_name)
      summaryDt <- summary_result$datatable
      
      fcHurdle <- merge(
        summaryDt[contrast == contrast_name & component == 'H', .(primerid, `Pr(>Chisq)`)],
        summaryDt[contrast == contrast_name & component == 'logFC', .(primerid, coef, ci.hi, ci.lo)],
        by = 'primerid'
      )
      if(nrow(fcHurdle) == 0){
        cat("  Warning: no genes passed filtering\n")
        next
      }
      
      fcHurdle[, padjust := p.adjust(`Pr(>Chisq)`, 'bonferroni')]
      fcHurdle[[group.by]] <- cell_group
      
      all_results <- rbind(all_results, fcHurdle)
      
      if(save.as.tmp){
        if (region == "") {
          region_suffix <- ""
        } else {
          region_suffix <- paste0(region,"_")
        }
        save(fcHurdle, file = paste0("tmp.degene.result.",region_suffix, make.names(cell_group), ".rda"))
        write.csv(fcHurdle, file = paste0("tmp.degene.result.",region_suffix, make.names(cell_group), ".csv"))
      }
    }, error = function(e){
      cat("  Error in cell group ", cell_group, " : ", e$message, "\n")
    })
  }
  
  if (nrow(all_results) == 0) {
    warning("No differential expression results were generated.")
  }
  
  return(all_results)
}
seo <- subset(seo, subset = celltype.L1_ct != "OPC")
r.deg_M <- perform_mast_analysis(seo,
                                group.by = "celltype.L2",
                                compare.by = "expriment",
                                group1 = "MW",
                                group2 = "MC",
                            batch.by = "date")

# colnames(r.deg_M) = c("gene","Pr","avg_log2FC","ci.hi","ci.lo","p_val_adj","celltype")
# write.csv(r.deg_M,"MAST_DEG_M.csv")

# filtered_r.deg_M = r.deg_M %>%
#   filter(p_val_adj < 0.05)

# write.csv(filtered_r.deg_M,"MAST_DEG_M_filtered.csv")