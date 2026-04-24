# %%
#library('clusterProfiler')
#library(org.Mm.eg.db)
library(ggplot2)
library(gprofiler2)
library(optparse)
library(dplyr)


option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "/data7/mark/STG/dataset/snRNA/merge_SCH/SUS_500_1000gene/merged_mast_wilcox_degs_all_log2fc0.1/merged_mast_wilcox_degs_all_log2fc0.1.csv",
              help = "输入数据文件路径", metavar = "character"),
  make_option(c("-o", "--outdir"), type = "character", default = "/data7/mark/STG/dataset/snRNA/merge_SCH/SUS_500_1000gene/merged_mast_wilcox_degs_all_log2fc0.1/GO_enrichment_results2",
              help = "输出文件名称 [默认值: %default]", metavar = "character"),
  make_option(c("-r", "--region"), type = "character", default = "AMY", 
              help = "输出文件名称 [默认值: %default]", metavar = "character"),
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
              help = "开启详细输出模式")
)
parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

# Pass a parameter to the gost funct
#args <- commandArgs(trailingOnly = TRUE)
# cat("First argument:", args[1], "\n")
# brainregion <- args[1]  # e.g., "TH" or "AMY"
# %%
indir <- opt$input
outdir <- opt$outdir
region <- opt$region

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}
# df_important_TF <- read.csv(paste0(indir, '/TFregulon_analysis_important_TF_gender_new.csv'))
# df_grns <- list.files('/data2st1/junyi/output/atac0627/snregulation/ctx/', pattern = 'ctx_.*\\.tsv', full.names = TRUE)
overwrite <- FALSE
df_deg <- read.csv(indir)

# add a column named ctstatus, ct+status
# replace the space to "_",replace "/" to "-"

if(!'region' %in% colnames(df_deg)){
  if ('Region' %in% colnames(df_deg)){
    df_deg$region <- df_deg$Region

  }else{
    print("No region column found!")
    stop()
  }
}

if (region %in% unique(df_deg$region)) {
  df_deg <- df_deg[df_deg$region == region, ]
} else {
  print(paste("Region", region, "not found in df_deg$region!"))
  stop()
}

if(!'status' %in% colnames(df_deg)){
  if ('Direction' %in% colnames(df_deg)){
    df_deg$status <- df_deg$Direction
  }else{
    print("No status column found!")
    stop()
  }
} 
if(!'gene' %in% colnames(df_deg)){
  if ('Gene' %in% colnames(df_deg)){
    df_deg$gene <- df_deg$Gene
  }else{
    print("No status column found!")
    stop()
  }
  if ('gene_name' %in% colnames(df_deg)){
    df_deg$gene <- df_deg$gene_name
  }else{
    print("No status column found!")
    stop()
  }

} 

if(!'sex' %in% colnames(df_deg)){
  if ('Sex' %in% colnames(df_deg)){
    df_deg$sex <- df_deg$Sex
  }else{
    print("No status column found!")
    stop()
  }
}

if (!'Region_subclass' %in% colnames(df_deg)) {
  if ('Region Subclass' %in% colnames(df_deg)) {
    df_deg$Region_subclass <- df_deg$`Region Subclass`
  } else {
    if (all(c("region", "Subclass", "Neurotransmitter") %in% colnames(df_deg))) {
      df_deg$Region_subclass <- ifelse(
        df_deg$Neurotransmitter == "NN",
        paste(df_deg$region, df_deg$Subclass, sep = " "),
        df_deg$Subclass
      )
    } else {
      print("No Region_subclass column found and cannot construct from region, Subclass, Neurotransmitter!")
      stop()
    }
  }
}


df_deg$ctstatus <- paste0(df_deg$Region_subclass,"_",df_deg$sex,"_", df_deg$status)
df_deg$ctstatus <- gsub(" ", "_", df_deg$ctstatus)
df_deg$ctstatus <- gsub("/", "-", df_deg$ctstatus)

# %%
all_results <- list()
for (ctstatus in unique(df_deg$ctstatus)) {
  print(ctstatus)

  # if file already exists and overwrite is FALSE, skip
  #df_grn <- read.table(file, header = TRUE, sep = "\t")
  if (file.exists(paste(outdir, "GO_enrichment",ctstatus,".csv",sep = "_")) && !overwrite) {
    print(paste("File already exists:", paste(outdir, "GO_enrichment",ctstatus,".csv",sep = "_")))
    next
  }
    tryCatch(
    {
      df_deg_subset <- df_deg[df_deg$ctstatus == ctstatus, ]
      gene_type_label <- "gene"
      gene_list <- unique(df_deg_subset$gene)
      Direction <- unique(df_deg_subset$status)
      Neurotransmitter <- unique(df_deg_subset$Neurotransmitter)
      Region_subclass <- unique(df_deg_subset$Region_subclass)
      Subclass <- unique(df_deg_subset$Subclass)
      Sex <- unique(df_deg_subset$Sex)
      Region <- unique(df_deg_subset$region)

      if(length(gene_list) > 0) {
  
        go_results <- gost(query = gene_list,
                        organism = "mmusculus",
                        sources = c("GO:MF", "GO:CC", "GO:BP"),
                        ordered_query = TRUE,
                        significant = TRUE,
                        user_threshold = 0.05,
                        correction_method = "g_SCS",
                        evcodes = TRUE)
        
        if(!is.null(go_results$result) && nrow(go_results$result) > 0) {
          res <- go_results$result
          res$Sex <- Sex
          res$Subclass <- Subclass
          res$Region <- Region
          res$`Region subclass` <- Region_subclass
          res$Neurotransmitter <- Neurotransmitter
          res$Direction <- Direction
          res$Region_subclass_sex_direction <- ctstatus
          res$gene_type_group <- gene_type_label
          res = subset(res, p_value<=0.05)
          res$nlog10p <- -log10(res$p_value)
          res$nlog10p <- ifelse(res$Direction == "Down", -res$nlog10p, res$nlog10p)

          for (col_name in names(res)) {
            if (class(res[[col_name]]) == "list") {
              # Convert list elements to a single string, e.g., comma-separated
              res[[col_name]] <- sapply(res[[col_name]], function(x) paste(x, collapse = ","))
            }
          }

          write.csv(res, file = file.path(outdir,paste0("GO_enrichment_",ctstatus,".csv")), row.names = FALSE)

          all_results[[ctstatus]] <- res

        } else {
          cat("No enrichment results for gene_type =", gene_type_label, "celltype =", ctstatus, "\n")
        }
      } else {
        cat("No genes for gene_type =", gene_type_label, "celltype =", ctstatus, "\n")
      }
    }, 
    error = function(e) {
      print(paste("Error in subsetting data for", ctstatus, Sex))
  })
}

if (length(all_results) > 0) {
  full_results <- bind_rows(all_results)
  write.csv(full_results, file = file.path(outdir,paste0("GO_enrichment_",region,"_all.csv")), row.names = FALSE)
  print(paste("Full results saved to", file.path(outdir,paste0("GO_enrichment_",region,"_all.csv"))))
} else {
  print("No results to save in full CSV.")
}