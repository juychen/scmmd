library(data.table)
library(limma)

## 1. 导入元数据 --------------------------------------------------------------
meta <- fread("/data1st1/yejun/metadata.csv")
meta <- as.data.frame(meta)
rownames(meta) <- meta$V1
meta <- meta[, -1]

## 2. 遍历 AUCell 结果 --------------------------------------------------------
auc_dir  <- "/data1st2/yejun/pyscenic/aucell"
auc_files <- list.files(auc_dir, pattern = "\\.tsv$", full.names = TRUE)

t_list <- list()   # 用来存每个文件的 t 值向量

for (f in auc_files) {
  ## 2.1 读取单个 AUCell 结果
  df <- fread(f)
  df <- as.data.frame(df)
  rownames(df) <- df$Cell          # 细胞在行
  df <- df[, -1]                   # 去掉 Cell 列
  head(row.names(meta))
  head(row.names(df))
  ## 2.2 匹配元数据
  mt <- meta[match(rownames(df), rownames(meta)), , drop = FALSE]
  
  ## 2.3 limma 差异分析
  expr   <- as.matrix(df)          # 行 = 细胞，列 = regulon
  group  <- factor(mt$status, levels = c("C", "W"))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- c("C", "W")
  
  fit  <- lmFit(t(expr), design)   # 注意转置：行 = regulon
  fit  <- contrasts.fit(fit, makeContrasts(WvsC = C - W, levels = design))
  fit  <- eBayes(fit)
  
  res <- topTable(fit, coef = "WvsC", number = Inf, sort.by = "none")
  
  ## 2.4 保存 t 值
  celltype <- tools::file_path_sans_ext(basename(f))
  t_list[[celltype]] <- setNames(res$t, rownames(res)) 
}

## 3. 合并所有 t 值 -----------------------------------------------------------
all_genes <- unique(unlist(lapply(t_list, names)))
merged_mat <- matrix(
  NA_real_,
  nrow = length(all_genes),
  ncol = length(t_list),
  dimnames = list(all_genes, names(t_list))
)

for (ct in names(t_list)) {
  merged_mat[names(t_list[[ct]]), ct] <- t_list[[ct]]
}

## ---------- 下面是新增部分 ---------------------------------------------------
merged_mat[is.na(merged_mat)] <- 0              # ① NA ➜ 0
merged_df <- as.data.frame(merged_mat)
colnames(merged_df) <- sub("^auc_", "",        # ③ 去掉 "auc_" 前缀
                           colnames(merged_df))
write.xlsx(merged_df, "/data1st2/yejun/pyscenic/aucell_p_values_merged_filtered.xlsx", rowNames = TRUE)
keep <- apply(abs(merged_df) >= 2, 1, any)      # ② 至少一列 |t| ≥ 2
filtered_df <- merged_df[keep, , drop = FALSE]

colnames(filtered_df) <- sub("^auc_", "",        # ③ 去掉 "auc_" 前缀
                             colnames(filtered_df))

heatmap(as.matrix(filtered_df))
colnames(filtered_df) 


## 4. 导出 --------------------------------------------------------------------
out_file <- file.path(auc_dir,
                      "aucell_t_values_merged_filtered.csv")
library(openxlsx)
write.csv(filtered_df, "/data1st2/yejun/pyscenic/aucell_t_values_merged_filtered_new.csv")
# Write the data to an Excel file with row names
write.xlsx(filtered_df, "/data1st2/yejun/pyscenic/aucell_t_values_merged_filtered.xlsx", rowNames = TRUE)