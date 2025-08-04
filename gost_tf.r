# %%
library('clusterProfiler')
library(org.Mm.eg.db)
library(ggplot2)
library(gprofiler2)

# Pass a parameter to the gost funct

# %%
indir <- '/data2st1/junyi/output/atac0627/snregulation'
outdir <- '/data2st1/junyi/output/atac0627/snregulation/gost/'
df_important_TF <- read.csv(paste0(indir, '/TFtarget_analysis_important_TF.csv'))
df_grns <- list.files('/data1st1/yejun/pyscenic/grn/', pattern = 'adj_.*\\.tsv', full.names = TRUE)
overwrite <- FALSE

# %%
df_important_TF.unique <- unique(df_important_TF$TF)

# %%

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
#df_grns <- df_grns[1:3]  # Adjust this line to read all files if needed
for (file in df_grns) {
  celltype <- gsub("^.*adj_([A-Za-z0-9_]+)_.*$", "\\1", file) 
  gender <-gsub(".*_([A-Za-z]+)\\.tsv", "\\1", file)
  region <- gsub("^.*adj_([A-Za-z0-9]+)_([A-Za-z0-9_]+).*$", "\\1", file) 

  # if file already exists and overwrite is FALSE, skip
  df_grn <- read.table(file, header = TRUE, sep = "\t")

  # filter the GRN dataframe to only include important TFs
  df_grn <- df_grn[df_grn$TF %in% df_important_TF.unique, ]

  for (tf in unique(df_grn$TF)) {

      if (file.exists(paste(outdir, "GO_enrichment",celltype,gender,tf,".csv",sep = "_")) && !overwrite) {
        next
      }
        tryCatch(
        {
          subset_df <- df_grn[df_grn$TF == tf, ]
          gene_type_label = "gene"
          gene_list <- unique(subset_df$target)
          if(length(gene_list) > 0) {
      
            gostres <- gost(query = gene_list,
                            organism = "mmusculus",
                            sources = c("GO:MF", "GO:CC", "GO:BP"),
                            ordered_query = TRUE,
                            significant = TRUE,
                            user_threshold = 0.05,
                            correction_method = "g_SCS",
                            evcodes = TRUE)
            
            if(!is.null(gostres$result) && nrow(gostres$result) > 0) {
              res <- gostres$result
              res$celltype <- celltype
              res$gene_type_group <- gene_type_label


              for (col_name in names(res)) {
                if (class(res[[col_name]]) == "list") {
                  # Convert list elements to a single string, e.g., comma-separated
                  res[[col_name]] <- sapply(res[[col_name]], function(x) paste(x, collapse = "/"))
              }}
              write.csv(res, file = paste(outdir, "GO_enrichment",celltype,gender,tf,".csv",sep = "_"), row.names = FALSE)

              #all_results[[celltype]] <- res
            } else {
              cat("No enrichment results for gene_type =", gene_type_label, "celltype =", celltype, "\n")
            }
          } else {
            cat("No genes for gene_type =", gene_type_label, "celltype =", celltype, "\n")
          }
        }, 
        error = function(e) {
          print(paste("Error in subsetting data for", celltype, gender))
      })
  }
  
}