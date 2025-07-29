# %%
library('clusterProfiler')
library(org.Mm.eg.db)
library(ggplot2)

# %%
indir <- '/data2st1/junyi/output/atac0627/snregulation'
outdir <- '/data2st1/junyi/output/atac0627/snregulation/enrich/'
df_important_TF <- read.csv(paste0(indir, '/TFtarget_analysis_important_TF.csv'))
df_grns <- list.files('/data1st1/yejun/pyscenic/grn/', pattern = 'adj_.*\\.tsv', full.names = TRUE)
overwrite <- FALSE

# %%
df_important_TF.unique <- unique(df_important_TF$TF)

# %%
# for file in df_grns[:5]:
#     region = file.split("adj_")[1].split("_")[0]
#     gender = file.split(".")[0].split("_")[-1]
#     celltype = file.split("adj_")[1].replace(f'_{gender}.tsv','')

#     # region = 'TH'
#     # gender = 'M'
#     # celltype = 'Tnc_Adgrf5_Glut'
#     df_grn = pd.read_csv(file, sep='\t')
# for loop to read GRN files
# subset 5 files for testing
for (file in df_grns) {
  celltype <- gsub("^.*adj_([A-Za-z0-9_]+)_.*$", "\\1", file) 
  gender <-gsub(".*_([A-Za-z]+)\\.tsv", "\\1", file)
  region <- gsub("^.*adj_([A-Za-z0-9]+)_([A-Za-z0-9_]+).*$", "\\1", file) 

  # if file already exists and overwrite is FALSE, skip
  df_grn <- read.table(file, header = TRUE, sep = "\t")

  # filter the GRN dataframe to only include important TFs
  df_grn <- df_grn[df_grn$TF %in% df_important_TF.unique, ]

  # for each tf in df_grn.TF.unique:
  for (tf in unique(df_grn$TF)) {

      if (file.exists(paste(outdir, "GO_enrichment_memento",celltype,gender,tf,".csv",sep = "_")) && !overwrite) {
        next
      }


        tryCatch(
        {
          print(paste("Processing:", celltype, gender))
          subset_df <- df_grn[df_grn$TF == tf, ]

          print(paste("Number of genes:", nrow(subset_df)))

          gene_entrez <- bitr(unique(subset_df$target), 
                              fromType = "SYMBOL", 
                              toType = "ENTREZID", 
                              OrgDb = "org.Mm.eg.db")
                              
          entrez_ids <- gene_entrez$ENTREZID
          ego <- enrichGO(gene = entrez_ids,
                          OrgDb = org.Mm.eg.db,
                          keyType = "ENTREZID",
                          ont = "ALL",  # "BP", "MF", "CC" 或 "ALL"
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05,
                          readable = TRUE)  # 结果中显示基因符号

          # 查看结果
          write.csv(ego, file = paste(outdir, "GO_enrichment_memento",celltype,gender,tf,".csv",sep = "_"), row.names = FALSE)
          p <- dotplot(ego, showCategory = 30)
          ggsave(paste(outdir, "GO_enrichment_memento_dot",celltype,gender,tf,".pdf",sep = "_"), plot = p, width = 12, height = 15)  # 宽度12英寸，高度15英寸

          p<-barplot(ego, showCategory = 30)  # 条形图
          ggsave(paste(outdir, "GO_enrichment_memento_bar",celltype,gender,tf,".pdf",sep = "_"), plot = p, width = 12, height = 15)  # 宽度12英寸，高度15英寸

          kk <- enrichKEGG(gene = entrez_ids,
                          organism = "mmu",  # 小鼠KEGG代码
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)

          # 可视化
          p<-dotplot(kk, showCategory = 30)
          write.csv(kk, file = paste(outdir, "KEGG_enrichment_memento",celltype,gender,tf,".csv",sep = "_"), row.names = FALSE)
          ggsave(paste(outdir, "kegg_enrichment_memento_dot",celltype,gender,tf,".pdf",sep = "_"), plot = p, width = 12, height = 15)  # 宽度12英寸，高度15英寸
        }, 
        error = function(e) {
          print(paste("Error in subsetting data for", celltype, gender))
      })
  }
  
}



