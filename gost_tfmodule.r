# %%
library('clusterProfiler')
library(org.Mm.eg.db)
library(ggplot2)
library(gprofiler2)

# Pass a parameter to the gost funct
args <- commandArgs(trailingOnly = TRUE)
cat("First argument:", args[1], "\n")
brainregion <- args[1]  # e.g., "TH" or "AMY"

# %%
indir <- '/data2st1/junyi/output/atac0627/snregulation'
outdir <- '/data2st1/junyi/output/atac0627/snregulation/gostmodule/'
df_important_TF <- read.csv(paste0(indir, '/genemodule.csv'))
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
  # if file name does not contain brainregion, skip
  # If it is non neuron cells
  if (brainregion=='NN') {
    regions <- c('HPF','Isocortex', 'AMY','PFC', 'TH', 'STR',  'HY', 'MB')
    # if the file name contains any of the regions, skip
    if (any(grepl(paste(regions, collapse = "|"), file))) {
      #print(paste("Skipping file:", file, "because it is not NN name."))
      next
    }
  } else if (!grepl(brainregion, file)) {
    next
  }
  # celltype <- gsub("^.*adj_([A-Za-z0-9_]+)_.*$", "\\1", file) 
  # gender <-gsub(".*_([A-Za-z]+)\\.tsv", "\\1", file)
  # region <- gsub("^.*adj_([A-Za-z0-9]+)_([A-Za-z0-9_]+).*$", "\\1", file) 

  region <- strsplit(strsplit(file, "adj_")[[1]][2], "_")[[1]][1]
  gender <- strsplit(strsplit(file, "\\.")[[1]][1], "_")[[1]][length(strsplit(strsplit(file, "\\.")[[1]][1], "_")[[1]])]
  celltype <- gsub(paste0("_", gender, "\\.tsv"), "", strsplit(file, "adj_")[[1]][2])


  # if file already exists and overwrite is FALSE, skip
  df_grn <- read.table(file, header = TRUE, sep = "\t")

  # filter the GRN dataframe to only include important TFs
  #df_grn <- df_grn[df_grn$TF %in% df_important_TF.unique, ]

for (module in unique(df_important_TF$module)) {
      tf <- paste0("TFmoudule",module)

      if (file.exists(paste(outdir, "GO_enrichment",celltype,gender,tf,".csv",sep = "_")) && !overwrite) {
        print(paste("File already exists:", paste(outdir, "GO_enrichment",celltype,gender,tf,".csv",sep = "_")))
        next
      }
        tryCatch(
        {

          tf_moduel <- df_important_TF[df_important_TF$module == module, ]
          subset_df <- df_grn[df_grn$TF %in% tf_moduel$TF, ]
          write.csv(subset_df, file = paste(outdir, "GRN_",celltype,gender,tf,".csv",sep = "_"), row.names = FALSE)
          gene_type_label <- "gene"
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