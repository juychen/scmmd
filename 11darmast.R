library(zellkonverter)
library(MAST)
library(SingleCellExperiment)
library(argparse)
library(parallel)
library(foreach)
library(doParallel)

parser <- ArgumentParser(description = "Description of your script")
parser$add_argument("--input", type = "character", required = TRUE, help = "Path to the input file")
parser$add_argument("--output", type = "character", required = TRUE, help = "Path to the output file")
parser$add_argument("--chunk", type = "integer", default = 30000, help = "Number of regions to process in each chunk")
parser$add_argument("--ncores", type = "integer", default = 8, help = "Number of cores to use for parallel processing")
parser$add_argument("--celltype_name", type = "character", default = "PFC_PFC_Glut", help = "Column name for region counts in the input file")

args <- parser$parse_args()
input <- args$input
output <- args$output
n_cores <- args$ncores
region_ct <- args$celltype_name


sce <- readH5AD(input)

n_genes <- nrow(sce)
batch_size <- min(args$chunk, n_genes)
batches <- split(1:n_genes, ceiling(seq_along(1:n_genes) / batch_size))

sce_subsets <- lapply(batches, function(idx) sce[idx, ])


n_cores <- min(n_cores, detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterExport(cl, c("assay", "FromMatrix", "zlm", "colData", "rowData","summary",'merge'))
results <- foreach(i = seq_along(sce_subsets), .combine = c) %dopar% {
  sce_subset <- sce_subsets[[i]]
  counts_matrix <- assay(sce_subset, "X")
  cell_metadata <- colData(sce_subset)
  gene_metadata <- rowData(sce_subset)

  print(paste("Processing subset", i, "of", length(sce_subsets)))
  
  sca <- FromMatrix(
    exprsArray = as.matrix(counts_matrix),
    cData = cell_metadata,
    fData = gene_metadata
  )

  print(paste("Creating SingleCellAssay for subset", i))
  
  zlmCond <- zlm(~expriment + fracHQP, sca)
  summaryCond <- summary(zlmCond, doLRT = "exprimentMW")
  summaryDt <- summaryCond$datatable

  print(paste("Summarizing results for subset", i))

  write.csv(
    summaryDt,
    file = paste0(output, "/",region_ct,"_part", i, ".csv")  )
  
  # fcHurdle <- merge(
  #   summaryDt[contrast == 'exprimentMW' & component == 'H', .(primerid, `Pr(>Chisq)`)],
  #   summaryDt[contrast == 'exprimentMW' & component == 'logFC', .(primerid, coef, ci.hi, ci.lo)],
  #   by = 'primerid'
  # )
  list(summaryCond)  # Return as list to combine later
}

# Clean up
stopCluster(cl)