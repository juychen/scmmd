# ============================================================================
# DNA Methylation Gene DEG Heatmap
# ============================================================================
# 读取 DEG Excel 文件，筛选 DNA 甲基化相关基因，用 ComplexHeatmap 画热图
# 行 = 基因, 列 = 细胞类型, 值 = log2FC 或 nlogp
# ============================================================================

set.seed(123)

library(tidyr)
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(readxl)
library(grDevices)
library(grid)

# ============================================================================
# 0. 配置参数
# ============================================================================

# 输入文件
DEG_FILE <- "/data2st2/junyi/output/1-1 six_dataset_DEG_list_fdr_FC01_filtered_7region_rmMB 20251227.xlsx"
DEG_SHEET <- "all"

# 输出文件前缀 (可自定义)
OUTPUT_PREFIX <- "figures/methylation_genes_heatmap"

# 筛选条件: 可选 "CUSUS M", "CUSUS F", "CSSUS M", "CURES M", "CURES F", "CSRES M" 或 NULL (不过滤)
# 用 NULL 表示不根据 Model 过滤
FILTER_MODEL <- NULL  # e.g. "CUSUS M"

# 筛选 Region: 可选 c("AMY","HPF","PFC") 或 NULL (不过滤)
FILTER_REGIONS <- NULL  # e.g. c("AMY","HPF","PFC")

# 值类型: "log2FC" 或 "nlogp"
VALUE_TYPE <- "log2FC"

# 热图参数
CLUSTER_ROWS <- TRUE
CLUSTER_COLS <- TRUE
SHOW_ROW_NAMES <- TRUE
SHOW_COL_NAMES <- TRUE
ROW_FONTSIZE <- 10
COL_FONTSIZE <- 6

# ============================================================================
# 1. DNA 甲基化基因列表
# ============================================================================

dna_methylation_dict <- list(
  "Writers"    = c("Dnmt1", "Dnmt3a", "Dnmt3b"),
  "Co-factors" = c("Dnmt3l", "Uhrf1", "Uhrf2"),
  "Readers"    = c("Mecp2", "Mbd1", "Mbd2", "Mbd3", "Mbd4"),
  "Erasers"    = c("Tet1", "Tet2", "Tet3"),
  "BER_Repair" = c("Tdg", "Neil1", "Neil2", "Smug1", "Mutyh", "Mbd4")
)

# 展平成基因 -> category 的映射表
gene_category <- stack(dna_methylation_dict)
colnames(gene_category) <- c("Gene", "Category")

# 所有基因名 (大写，用于匹配)
all_methylation_genes <- toupper(unique(gene_category$Gene))
cat(sprintf("Total methylation-related genes to search: %d\n", length(all_methylation_genes)))

# ============================================================================
# 2. 读取 DEG 数据
# ============================================================================

cat(sprintf("Reading DEG file: %s\n", DEG_FILE))
df_deg_all <- read_excel(DEG_FILE, sheet = DEG_SHEET)
df_deg_all <- as.data.frame(df_deg_all, stringsAsFactors = FALSE)

cat(sprintf("Loaded %d rows, %d columns\n", nrow(df_deg_all), ncol(df_deg_all)))
cat(sprintf("Columns: %s\n", paste(colnames(df_deg_all), collapse = ", ")))

# ============================================================================
# 3. 数据预处理
# ============================================================================

# 将第一列作为 gene 名称 (index_col=0 的逻辑)
# read_excel 默认第一列是数据列，需要把它设成行名
if (!("Gene" %in% colnames(df_deg_all))) {
  # 如果第一列不是 "Gene"，假设第一列是 gene 名
  first_col <- colnames(df_deg_all)[1]
  cat(sprintf("Setting first column '%s' as gene index\n", first_col))
  df_deg_all$Gene <- df_deg_all[[first_col]]
}

# 统一基因名为大写用于匹配
df_deg_all$Gene_upper <- toupper(df_deg_all$Gene)

# 筛选甲基化基因
df_meth <- df_deg_all[df_deg_all$Gene_upper %in% all_methylation_genes, ]
cat(sprintf("Found %d rows for methylation genes\n", nrow(df_meth)))

if (nrow(df_meth) == 0) {
  stop("No methylation genes found in DEG data! Check gene name case/format.")
}

# 打印找到的基因
found_genes <- unique(df_meth$Gene_upper)
cat(sprintf("Matched genes: %s\n", paste(found_genes, collapse = ", ")))
missing_genes <- setdiff(all_methylation_genes, found_genes)
if (length(missing_genes) > 0) {
  cat(sprintf("Missing genes (not in DEG data): %s\n", paste(missing_genes, collapse = ", ")))
}

# 可选: 按 Model 过滤
if (!is.null(FILTER_MODEL)) {
  df_meth <- df_meth[grepl(FILTER_MODEL, df_meth$Model, fixed = TRUE), ]
  cat(sprintf("After filtering Model='%s': %d rows\n", FILTER_MODEL, nrow(df_meth)))
}

# 可选: 按 Region 过滤
if (!is.null(FILTER_REGIONS)) {
  df_meth <- df_meth[df_meth$Region %in% FILTER_REGIONS, ]
  cat(sprintf("After filtering Regions=%s: %d rows\n",
              paste(FILTER_REGIONS, collapse = ","), nrow(df_meth)))
}

# 添加 Category 信息
df_meth <- merge(df_meth, gene_category, by.x = "Gene_upper", by.y = "Gene", all.x = TRUE)

# 创建细胞类型列 (sample = Model_Region_Subclass)
# 清理 Model 名: 空格换下划线
df_meth$Model_clean <- gsub(" ", "_", df_meth$Model)

if (all(c("Region", "Subclass") %in% colnames(df_meth))) {
  df_meth$ctname <- paste(df_meth$Region, df_meth$Subclass, sep = "_")
  df_meth$ctname <- gsub("/", "-", df_meth$ctname)
  df_meth$ctname <- gsub(" ", "_", df_meth$ctname)
} else if ("Region Subclass" %in% colnames(df_meth)) {
  df_meth$ctname <- df_meth[["Region Subclass"]]
  df_meth$ctname <- gsub("/", "-", df_meth$ctname)
  df_meth$ctname <- gsub(" ", "_", df_meth$ctname)
} else {
  stop("Cannot find Region/Subclass or 'Region Subclass' columns for cell type naming")
}

# sample = Model + celltype，确保每个列唯一
df_meth$sample <- paste(df_meth$Model_clean, df_meth$ctname, sep = "|")

# term_name = gene 名称
df_meth$term_name <- df_meth$Gene

# ============================================================================
# 4. 计算 value 列
# ============================================================================

if (VALUE_TYPE == "nlogp") {
  # nlogp = -log10(FDR) * sign(log2FC)
  df_meth$FDR_num <- suppressWarnings(as.numeric(as.character(df_meth$FDR)))
  df_meth$log2FC_num <- suppressWarnings(as.numeric(as.character(df_meth$log2FC)))

  # 处理 FDR=0 的情况, 用极小值替代
  df_meth$FDR_num[!is.na(df_meth$FDR_num) & df_meth$FDR_num == 0] <- 1e-300

  df_meth$value <- -log10(df_meth$FDR_num) * sign(df_meth$log2FC_num)
  df_meth$value[is.na(df_meth$value)] <- 0

  value_label <- "-log10(FDR) * sign(log2FC)"
  legend_name <- "nlogp"
} else {
  # log2FC
  df_meth$value <- suppressWarnings(as.numeric(as.character(df_meth$log2FC)))
  df_meth$value[is.na(df_meth$value)] <- 0

  value_label <- "log2FC"
  legend_name <- "log2FC"
}

# 去重: 如果同一 gene + sample (Model|celltype) 有多个值，取平均值
df_meth_agg <- df_meth %>%
  group_by(sample, term_name, Gene_upper, Category) %>%
  summarise(
    value = mean(value, na.rm = TRUE),
    Region = dplyr::first(Region),
    Model  = dplyr::first(Model_clean),
    .groups = "drop"
  )

cat(sprintf("After aggregation: %d rows\n", nrow(df_meth_agg)))

# ============================================================================
# 5. 构建热图矩阵 (行=基因, 列=细胞类型)
# ============================================================================

heatmap_data <- df_meth_agg %>%
  select(sample, term_name, value) %>%
  pivot_wider(names_from = term_name, values_from = value, values_fn = mean)

mat <- as.data.frame(heatmap_data)
rownames(mat) <- mat$sample
mat$sample <- NULL
mat <- as.matrix(mat)
mat[is.na(mat)] <- 0

# 转置: 行=基因, 列=细胞类型
mat_t <- t(mat)

cat(sprintf("Heatmap matrix: %d genes x %d celltypes\n", nrow(mat_t), ncol(mat_t)))

# ============================================================================
# 6. 列注释 (细胞类型) - 按 Model + Region
# ============================================================================

# 从 df_meth_agg 获取每个 sample 的 Model 和 Region
sample_meta <- df_meth_agg %>%
  select(sample, Region, Model) %>%
  distinct()

col_anno_df <- data.frame(
  sample = colnames(mat_t),
  stringsAsFactors = FALSE
)
col_anno_df <- merge(col_anno_df, sample_meta, by = "sample", all.x = TRUE)
rownames(col_anno_df) <- col_anno_df$sample
col_anno_df$sample <- NULL

# 如果缺失，从 sample 名解析: "Model|Region_Subclass"
if (any(is.na(col_anno_df$Model))) {
  parsed <- strsplit(rownames(col_anno_df), "\\|")
  col_anno_df$Model  <- sapply(parsed, `[`, 1)
  col_anno_df$Region <- gsub("_.*", "", sapply(parsed, `[`, 2))
}

# ============================================================================
# 7. 行注释 (基因) - 按 Category
# ============================================================================

row_anno_df <- data.frame(
  gene = rownames(mat_t),
  stringsAsFactors = FALSE
)
gene_to_cat <- df_meth_agg %>%
  select(term_name, Category) %>%
  distinct()
colnames(gene_to_cat)[1] <- "gene"
row_anno_df <- merge(row_anno_df, gene_to_cat, by = "gene", all.x = TRUE)
rownames(row_anno_df) <- row_anno_df$gene
row_anno_df$gene <- NULL

# ============================================================================
# 8. 颜色定义 (复用项目的配色方案)
# ============================================================================

base_region <- c(
  AMY="#ff7f0e", STR="#e377c2", PFC="#8c564b", iCTX="#1f77b4",
  MB="#9467bd", TH="#bcbd22", HY="#d62728", HPF="#009E73"
)

# Category 颜色 (甲基化分类)
base_category <- c(
  "Writers"    = "#E41A1C",
  "Co-factors" = "#377EB8",
  "Readers"    = "#4DAF4A",
  "Erasers"    = "#984EA3",
  "BER_Repair" = "#FF7F00"
)

# Model 颜色 (复用项目配色)
base_model <- c(
  "CUSUS_M"  = "#00B0F0",
  "CUSUS_F"  = "#A5D86E",
  "CSSUS_M"  = "#2B6CB8",
  "CURES_M"  = "#E63CE6",
  "CURES_F"  = "#F6B1C3",
  "CSRES_M"  = "#BD551C",
  "SUS"      = "#46A6B2",
  "RES"      = "#DB6B94",
  "CUMS"     = "#A88ED3",
  "CSDS"     = "#7A4A6F"
)

# 动态调色板函数
make_palette <- function(values, base_map) {
  vals <- as.character(unique(values))
  vals <- vals[!is.na(vals)]
  base_sub <- base_map[intersect(names(base_map), vals)]
  missing <- setdiff(vals, names(base_map))
  if (length(missing) > 0) {
    extra <- grDevices::hcl(
      h = seq(15, 375, length.out = length(missing) + 1)[-1],
      c = 70, l = 65
    )
    names(extra) <- missing
    base_sub <- c(base_sub, extra[missing])
  }
  base_sub[match(vals, names(base_sub))]
}

dyn_colors <- list()
if ("Region" %in% colnames(col_anno_df) && length(na.omit(col_anno_df$Region)) > 0) {
  dyn_colors$Region <- make_palette(col_anno_df$Region, base_region)
}
if ("Model" %in% colnames(col_anno_df) && length(na.omit(col_anno_df$Model)) > 0) {
  dyn_colors$Model <- make_palette(col_anno_df$Model, base_model)
}
if ("Category" %in% colnames(row_anno_df) && length(na.omit(row_anno_df$Category)) > 0) {
  dyn_colors$Category <- make_palette(row_anno_df$Category, base_category)
}

# ============================================================================
# 9. 热图颜色映射 (蓝-白-红)
# ============================================================================

rng <- range(mat_t, finite = TRUE)
m <- max(abs(rng))
at <- seq(-m, m, length.out = 5)

col_fun <- circlize::colorRamp2(
  c(min(at), 0, max(at)),
  c("#3498db", "#FFFFFF", "#d62728")
)

legend_param <- list(
  at = at,
  labels = format(at, digits = 2),
  title = legend_name,
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 8)
)

# ============================================================================
# 10. 构建注释
# ============================================================================

# 列注释 (Model + Region)
col_annotation_list <- list()
if ("Model" %in% colnames(col_anno_df)) {
  col_annotation_list$Model <- col_anno_df$Model
}
if ("Region" %in% colnames(col_anno_df)) {
  col_annotation_list$Region <- col_anno_df$Region
}

col_ha <- NULL
col_split <- NULL
if (length(col_annotation_list) > 0) {
  col_ha <- HeatmapAnnotation(
    df = as.data.frame(col_annotation_list),
    col = dyn_colors,
    annotation_name_side = "right",
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
  )
  # 列分割按 Model
  if ("Model" %in% colnames(col_anno_df)) {
    col_split <- col_anno_df$Model
  }
}

# 行注释
row_annotation_list <- list()
if ("Category" %in% colnames(row_anno_df)) {
  row_annotation_list$Category <- row_anno_df$Category
}

row_ha <- NULL
row_split <- NULL
if (length(row_annotation_list) > 0) {
  row_ha <- rowAnnotation(
    df = as.data.frame(row_annotation_list),
    col = dyn_colors,
    annotation_name_side = "top",
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
  )
  # 行分割按 Category
  if ("Category" %in% colnames(row_anno_df)) {
    row_split <- row_anno_df$Category
  }
}

# ============================================================================
# 11. 画热图
# ============================================================================

# 自动计算图片尺寸 (宽窄版)
plot_height <- max(6, nrow(mat_t) / 3 + 4)
plot_width  <- max(6, ncol(mat_t) / 5 + 4)

output_pdf <- sprintf("%s_%s.pdf", OUTPUT_PREFIX, VALUE_TYPE)
cat(sprintf("Saving heatmap to: %s\n", output_pdf))
cat(sprintf("Plot size: %.1f x %.1f inches\n", plot_height, plot_width))

pdf(output_pdf, height = plot_height, width = plot_width)

ht <- Heatmap(
  mat_t,
  name = legend_name,
  col = col_fun,
  heatmap_legend_param = legend_param,
  show_row_names = SHOW_ROW_NAMES,
  show_column_names = SHOW_COL_NAMES,
  cluster_rows = CLUSTER_ROWS,
  cluster_columns = CLUSTER_COLS,
  cluster_column_slices = FALSE,
  row_split = row_split,
  column_split = col_split,
  left_annotation = row_ha,
  top_annotation = col_ha,
  row_names_side = "right",
  column_names_side = "bottom",
  row_names_gp = gpar(fontsize = ROW_FONTSIZE),
  column_names_gp = gpar(fontsize = COL_FONTSIZE),
  column_gap = unit(1, "mm"),
  use_raster = FALSE,
  row_title = NULL,
  column_title = sprintf("DNA Methylation Genes - %s", VALUE_TYPE)
)

ht_drawn <- draw(ht)
dev.off()

# ============================================================================
# 12. 保存数据矩阵
# ============================================================================

output_csv <- sprintf("%s_%s_data.csv", OUTPUT_PREFIX, VALUE_TYPE)
write.csv(mat_t, output_csv, quote = TRUE)
cat(sprintf("Data matrix saved to: %s\n", output_csv))

# 打印矩阵摘要
cat(sprintf("\n=== Heatmap matrix summary ===\n"))
cat(sprintf("Genes (rows): %s\n", paste(rownames(mat_t), collapse = ", ")))
cat(sprintf("Cell types (cols): %d\n", ncol(mat_t)))
cat(sprintf("Value range: [%.3f, %.3f]\n", min(mat_t), max(mat_t)))

cat("\nDone!\n")
