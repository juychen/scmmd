#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(fgsea)
})

# ========== 配置 ==========
#gmt_file <- "/data1st2/hannan_25/data/RNAexp/snRNA/geneset/from_Gprofile/mmusculus.ALL.prefixed.merged.gmt"
gmt_file="/data1st2/hannan_25/data/RNAexp/snRNA/geneset/from_Gprofile/mmusculus.BPMF.prefixed.merged.gmt"

# fgsea 参数
minSize <- 5
maxSize <- 5000
nperm   <- 10000



# ========== 自定义函数：读取 GMT（保留第2列通路名） ==========
# 返回：list(pathways=list(id=c(genes...)), meta=data.frame(id, name2, source_hint))
read_gmt_with_names <- function(gmt_path) {
  con <- file(gmt_path, open = "r")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  pathways <- list()
  name2_map <- list()
  for (ln in lines) {
    if (!nzchar(ln)) next
    parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 2) next
    pid <- parts[[1]]                   # 第1列：通路ID（如 GO:0007399、REAC:R-HSA-...、WP:WP88、CORUM:...）
    pname2 <- parts[[2]]                # 第2列：通路名（你已加了前缀，如 GO:BP_nervous system development）
    genes <- unique(parts[3:length(parts)])
    genes <- genes[nzchar(genes)]
    if (length(genes) == 0) next
    pathways[[pid]] <- genes
    name2_map[[pid]] <- pname2
  }
  meta <- data.frame(
    term_id = names(name2_map),
    term_name_2ndcol = unlist(name2_map, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  # 从第二列前缀推断来源（更可靠：你已加了 GO:BP_/GO:MF_/GO:CC_/CORUM_/REAC_/HP/WP）
  meta$source <- dplyr::case_when(
    str_starts(meta$term_name_2ndcol, "GO:BP_") ~ "GO:BP",
    str_starts(meta$term_name_2ndcol, "GO:MF_") ~ "GO:MF",
    str_starts(meta$term_name_2ndcol, "GO:CC_") ~ "GO:CC",
    str_starts(meta$term_name_2ndcol, "CORUM_") ~ "CORUM",
    str_starts(meta$term_name_2ndcol, "REAC_")  ~ "REAC",
    str_starts(meta$term_name_2ndcol, "HP")     ~ "HP",
    str_starts(meta$term_name_2ndcol, "WP")     ~ "WP",
    TRUE ~ NA_character_
  )
  list(pathways = pathways, meta = meta)
}

message("[1/6] 读取 GMT：", gmt_file)
gmt <- read_gmt_with_names(gmt_file)
#gmt=gmt[gmt$meta$source%in%c("GO:BP","GO:MF")]
pathways <- gmt$pathways
meta_map <- gmt$meta
stopifnot(length(pathways) > 0)

# 有效域大小：GMT 所含全部基因的去重数
effective_domain_size <- length(unique(unlist(pathways, use.names = FALSE)))

# ========== 工具函数 ==========
# 根据 term_id 抓 second-column 名称（去掉前缀作为 term_name），并给出 source
term_info_from_meta <- function(term_id) {
  i <- match(term_id, meta_map$term_id)
  if (is.na(i)) return(list(term_name = term_id, source = NA_character_))
  name2 <- meta_map$term_name_2ndcol[i]
  src   <- meta_map$source[i]
  # 去掉前缀（GO:BP_/GO:MF_/GO:CC_/CORUM_/REAC_/HP/WP）
  name_stripped <- name2
  name_stripped <- sub("^GO:BP_", "", name_stripped)
  name_stripped <- sub("^GO:MF_", "", name_stripped)
  name_stripped <- sub("^GO:CC_", "", name_stripped)
  name_stripped <- sub("^CORUM_",  "", name_stripped)
  name_stripped <- sub("^REAC_",   "", name_stripped)
  name_stripped <- sub("^HP_?",    "", name_stripped) # 兼容 HP 或 HP_
  name_stripped <- sub("^WP_?",    "", name_stripped) # 兼容 WP 或 WP_
  list(term_name = name_stripped, source = src)
}

format_fgsea <- function(fg, ranks_len,
                         Sex, Subclass, Region, RegionSubclass, Neurotransmitter,
                         region_Subclass_sex_direction) {
    # 统一列结构定义
  empty <- data.frame(
    query = character(0),
    significant = logical(0),
    p_value = numeric(0),
    term_size = numeric(0),
    query_size = numeric(0),
    intersection_size = numeric(0),
    precision = numeric(0),
    recall = numeric(0),
    term_id = character(0),
    source = character(0),
    term_name = character(0),
    effective_domain_size = numeric(0),
    source_order = character(0),
    parents = character(0),
    evidence_codes = character(0),
    intersection = character(0),
    Sex = character(0),
    Subclass = character(0),
    Region = character(0),
    Region.Subclass = character(0),
    Neurotransmitter = character(0),
    Direction = character(0),
    Region_Subclass_sex_direction = character(0),
    gene_type_group = character(0),
    nlog10p = numeric(0),
    ES = numeric(0),
    NES = numeric(0),
    padj = numeric(0),
    nMoreExtreme = numeric(0),
    stringsAsFactors = FALSE
  )

  if (nrow(fg) == 0) return(empty)
  
  # fgsea 返回的 pathway 名就是我们给的 list 名（此处是 term_id）
  term_id <- fg$pathway
  
  # 从 meta 获取 term_name（第2列修饰后再去前缀）与 source
  ti <- lapply(term_id, term_info_from_meta)
  term_name <- vapply(ti, `[[`, "", "term_name")
  source    <- vapply(ti, `[[`, "", "source")
  
  # 计算 intersection/precision/recall
  intersection_size <- vapply(fg$leadingEdge, length, integer(1))
  precision <- ifelse(ranks_len > 0, intersection_size / ranks_len, NA_real_)
  recall    <- ifelse(fg$size  > 0, intersection_size / fg$size, NA_real_)
  # GSEA方向由 NES 定义
  direction_by_nes <- ifelse(fg$NES > 0, "Up",
                             ifelse(fg$NES < 0, "Down", "NA"))
  intersection <- vapply(fg$leadingEdge, function(x) paste(x, collapse = ","), "")
  
  data.frame(
    query = "query_1",
    significant = fg$padj < 0.05,
    p_value = fg$pval,
    term_size = fg$size,
    query_size = ranks_len,
    intersection_size = intersection_size,
    precision = precision,
    recall = recall,
    term_id = term_id,
    source = source,
    term_name = term_name,
    effective_domain_size = effective_domain_size,
    source_order = "",
    parents = "",
    evidence_codes = "",
    intersection = intersection,
    Sex = Sex,
    Subclass = Subclass,
    Region = Region,
    `Region Subclass` = RegionSubclass,
    Neurotransmitter = Neurotransmitter,
    Direction = direction_by_nes,
    Region_Subclass_sex_direction = region_Subclass_sex_direction,
    gene_type_group = "gene",
    nlog10p = -log10(pmax(fg$padj, .Machine$double.xmin)),
    ES = fg$ES,
    NES = fg$NES,
    padj = fg$padj,
    nMoreExtreme = fg$nMoreExtreme,
    stringsAsFactors = FALSE
  )
}

run_one_file_signed_logp <- function(csv_file, pathways, out_dir) {
  message("[2/6] 读取：", csv_file)
  dt <- fread(csv_file)
  
  needed_cols <- c("Gene", "gene_name", "log2FC", "FDR",
                   "Subclass", "Region", "Sex",
                   "Neurotransmitter", "Region Subclass")
  miss <- setdiff(needed_cols, names(dt))
  if (length(miss) > 0) {
    stop("文件缺少列：", paste(miss, collapse = ", "), "\n文件：", csv_file)
  }
  
  # 基因名列（优先 gene_name）
  gene_col <- if ("gene_name" %in% names(dt)) "gene_name" else "Gene"
  dt[[gene_col]] <- str_replace_all(dt[[gene_col]], '^"+|"+$', "")
  
  # 分细胞类型（Region Subclass + Sex）
  dt <- dt %>% mutate(Region_Subclass_sex = paste0(`Region Subclass`, "_", Sex))
  celltypes <- unique(dt$Region_Subclass_sex)
  message("[3/6] 细胞类型数量：", length(celltypes))
  
  all_results <- list()
  gene_sym <- rlang::sym(gene_col)
  
  for (ct in celltypes) {
    tryCatch({
    sub <- dt %>% dplyr::filter(Region_Subclass_sex == ct)
    if (nrow(sub) == 0) next
    
    u <- function(v) { x <- unique(na.omit(v)); if (length(x)==0) NA_character_ else x[1] }
    Sex0            <- u(sub$Sex)
    Subclass0       <- u(sub$Subclass)
    Region0         <- u(sub$Region)
    RegionSubclass0 <- u(sub$`Region Subclass`)
    Neuro0          <- u(sub$Neurotransmitter)
    
    # -------- 关键：构建 sign(log2FC) * -log10(FDR) 排序（稳健版） --------
    ranks_tbl <- sub %>%
      dplyr::select(!!gene_sym, log2FC, FDR) %>%
      dplyr::filter(!is.na(log2FC) & !is.na(FDR) & !is.na(!!gene_sym)) %>%
      # 1）把 FDR 限制在 [1e-15, 1-1e-15] 区间，避免 -log10(0/1) 出 Inf / NaN
      dplyr::mutate(
        FDR = ifelse(is.na(FDR), NA_real_,
                      pmin(pmax(FDR, 1e-15), 1 - 1e-15)),
        stat = sign(log2FC) * (-log10(FDR))
      ) %>%
      # 2）去掉非有限的 stat（Inf / -Inf / NaN）
      dplyr::filter(!is.na(stat) & is.finite(stat)) %>%
      dplyr::group_by(!!gene_sym) %>%
      dplyr::slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
    
    if (nrow(ranks_tbl) == 0) {
      message("  [跳过] ", ct, "：没有有效的 log2FC+FDR（或全部被过滤为非有限值）")
      next
    }
    
    ranks_vec <- ranks_tbl$stat
    names(ranks_vec) <- ranks_tbl[[rlang::as_string(gene_sym)]]
    
    # 再兜一层底，保证传给 fgsea 的 stats 全是有限数
    finite_idx <- is.finite(ranks_vec)
    if (any(!finite_idx)) {
      bad_n <- sum(!finite_idx)
      message("  [警告] ", ct, " 中存在 ", bad_n, " 个非有限 stats，已移除。")
      ranks_vec <- ranks_vec[finite_idx]
    }
    
    if (length(ranks_vec) == 0) {
      message("  [跳过] ", ct, "：过滤非有限值后无可用基因")
      next
    }
    
    ranks_vec <- sort(ranks_vec, decreasing = TRUE)
    message("  [", ct, "] 使用 sign(log2FC)*-log10(FDR) 排序，基因数：", length(ranks_vec))
    
    if (length(ranks_vec) < minSize) {
      message("  [跳过] ", ct, "：可用基因不足（", length(ranks_vec), "）")
      next
    }
    
    fg <- suppressWarnings(
      fgsea(pathways = pathways,
            stats = ranks_vec,
            minSize = minSize,
            maxSize = maxSize,
            nperm = nperm)
    )
    fg <- as.data.frame(fg)
    
    rssd <- paste0(RegionSubclass0, "_", Sex0, "_GSEA")

    print(rssd)
    
    res <- format_fgsea(
      fg = fg,
      ranks_len = length(ranks_vec),
      Sex = Sex0,
      Subclass = Subclass0,
      Region = Region0,
      RegionSubclass = RegionSubclass0,
      Neurotransmitter = Neuro0,
      region_Subclass_sex_direction = rssd
    )
    res$Region_Subclass_sex_direction <- paste0(res$Region.Subclass,"_",res$Sex,"_",res$Direction)
    
    ct_safe <- str_replace_all(ct, "[^A-Za-z0-9_\\-\\.]", "_")
    p <- gsub(".csv","",basename(csv_file))
    dir.create(paste0(out_dir,"/",p), showWarnings = FALSE, recursive = TRUE)
    
    per_ct_out_up <- file.path(out_dir, paste0(p, "/GO_gsea_", ct_safe, "_Up.csv"))
    per_ct_out_down <- file.path(out_dir, paste0(p, "/GO_gsea_", ct_safe, "_Down.csv"))
    
    #res <- res[res$nlog10p > 2, ]
    res_up <- res[res$NES > 0, ]
    res_down <- res[res$NES < 0, ]
    
    fwrite(res_up, per_ct_out_up)
    fwrite(res_down, per_ct_out_down)
    all_results[[ct]] <- res
    
    message("  [完成] ", ct, " -> ", per_ct_out_up, " (", nrow(res_up), " 行)")
    message("  [完成] ", ct, " -> ", per_ct_out_down, " (", nrow(res_down), " 行)")
     }, error = function(e) {
    message("  [错误跳过] celltype: ", ct)
    message("    错误信息: ", conditionMessage(e))
  })
  }
  
  # 如果所有 celltype 都被跳过，避免 rbindlist 空报错
  if (length(all_results) == 0) {
    message("[4/6] 警告：", csv_file, " 中没有任何细胞类型通过过滤，跳过合并写出。")
    return(invisible(NULL))
  }
  
  p <- gsub(".csv","",basename(csv_file))
  merged_out <- file.path(out_dir, paste0(p, "/GO_gsea.csv"))
  merged <- data.table::rbindlist(all_results, use.names = TRUE, fill = TRUE)
  fwrite(merged, merged_out)
  message("[4/6] 合并结果：", merged_out, " (", nrow(merged), " 行)")
}


run_one_file <- function(csv_file, pathways, out_dir) {
  message("[2/6] 读取（主方法 MAST）：", csv_file)
  dt <- fread(csv_file)
  
  needed_cols <- c("Gene", "gene_name", "log2FC", "Subclass", "Region", "Sex",
                   "Neurotransmitter", "Region Subclass", "FDR")
  miss <- setdiff(needed_cols, names(dt))
  if (length(miss) > 0) {
    stop("MAST 文件缺少列：", paste(miss, collapse = ", "), "\n文件：", csv_file)
  }
  
  # ---------- 尝试找到对应的 wilcoxon 文件 ----------
  # 按你的路径结构改名：mast -> wilcoxon, _mast.csv -> _wilcoxon.csv
  other_csv <- csv_file
  other_csv <- sub(
    "mast_all_ngeneson_saturation_downsampled_ratio/all_ngeneson_saturation_downsampled_ratio",
    "wilcoxon_all_sample_downsampled_ratio/all_sample_downsampled_ratio",
    other_csv
  )
  other_csv <- sub("_mast\\.csv$", "_wilcoxon.csv", other_csv)
  
  has_other <- file.exists(other_csv)
  has_other <- FALSE  # 强制关闭配对文件功能，仅用单一 MAST 排序
  if (has_other) {
    message("  [2.1] 检测到配对 wilcoxon 文件：", other_csv)
    dt2 <- fread(other_csv)
    
    needed_cols2 <- c("Gene", "gene_name", "log2FC", "Subclass", "Region", "Sex",
                      "Neurotransmitter", "Region Subclass", "FDR")
    miss2 <- setdiff(needed_cols2, names(dt2))
    if (length(miss2) > 0) {
      warning("wilcoxon 文件缺少列：", paste(miss2, collapse = ", "),
              "\n文件：", other_csv, "\n该比较将退回单一 MAST 排序。")
      has_other <- FALSE
    }
  } else {
    message("  [2.1] 未找到配对 wilcoxon 文件，使用单一 MAST 排序。")
  }
  
  # 基因名列（优先 gene_name）
  gene_col <- if ("gene_name" %in% names(dt)) "gene_name" else "Gene"
  dt[[gene_col]] <- str_replace_all(dt[[gene_col]], '^"+|"+$', "")
  
  if (has_other) {
    dt2[[gene_col]] <- str_replace_all(dt2[[gene_col]], '^"+|"+$', "")
  }
  
  # 分细胞类型（Region Subclass + Sex）
  dt <- dt %>% mutate(Region_Subclass_sex = paste0(`Region Subclass`, "_", Sex))
  if (has_other) {
    dt2 <- dt2 %>% mutate(Region_Subclass_sex = paste0(`Region Subclass`, "_", Sex))
  }
  
  celltypes <- unique(dt$Region_Subclass_sex)
  message("[3/6] 细胞类型数量：", length(celltypes))
  
  all_results <- list()
  
  gene_sym <- rlang::sym(gene_col)
  
  for (ct in celltypes) {
    sub1 <- dt %>% dplyr::filter(Region_Subclass_sex == ct)
    if (nrow(sub1) == 0) next
    
    if (has_other) {
      sub2 <- dt2 %>% dplyr::filter(Region_Subclass_sex == ct)
      if (nrow(sub2) == 0) {
        # 对应细胞类型在 wilcoxon 里没有，退回单方法
        use_meta <- FALSE
      } else {
        use_meta <- TRUE
      }
    } else {
      use_meta <- FALSE
    }
    
    u <- function(v) { x <- unique(na.omit(v)); if (length(x)==0) NA_character_ else x[1] }
    Sex0            <- u(sub1$Sex)
    Subclass0       <- u(sub1$Subclass)
    Region0         <- u(sub1$Region)
    RegionSubclass0 <- u(sub1$`Region Subclass`)
    Neuro0          <- u(sub1$Neurotransmitter)
    
    # ------------------ 核心：构建排序向量 ranks_vec ------------------
    if (use_meta) {
      # 1）提取两种方法的 log2FC 和 FDR
      g1 <- sub1 %>%
        dplyr::select(!!gene_sym, log2FC1 = log2FC, p1 = FDR) %>%
        dplyr::filter(!is.na(log2FC1) & !is.na(p1) & !is.na(!!gene_sym)) %>%
        dplyr::group_by(!!gene_sym) %>%
        dplyr::slice_max(order_by = abs(log2FC1), n = 1, with_ties = FALSE) %>%
        dplyr::ungroup()
      
      g2 <- sub2 %>%
        dplyr::select(!!gene_sym, log2FC2 = log2FC, p2 = FDR) %>%
        dplyr::filter(!is.na(log2FC2) & !is.na(p2) & !is.na(!!gene_sym)) %>%
        dplyr::group_by(!!gene_sym) %>%
        dplyr::slice_max(order_by = abs(log2FC2), n = 1, with_ties = FALSE) %>%
        dplyr::ungroup()
      
      merged <- full_join(g1, g2, by = rlang::as_string(gene_sym))
      
      # ======= 这里是关键修改：更稳的 p 截断 + 过滤非有限 z_meta =======
      merged <- merged %>%
        mutate(
          # 把 p 限制在 [1e-15, 1 - 1e-15]，避免 qnorm 给 Inf
          p1 = ifelse(is.na(p1), NA_real_,
                      pmin(pmax(p1, 1e-15), 1 - 1e-15)),
          p2 = ifelse(is.na(p2), NA_real_,
                      pmin(pmax(p2, 1e-15), 1 - 1e-15)),
          
          # 计算带符号 Z 分数
          z1 = ifelse(!is.na(p1) & !is.na(log2FC1),
                      sign(log2FC1) * qnorm(p1 / 2, lower.tail = FALSE),
                      NA_real_),
          z2 = ifelse(!is.na(p2) & !is.na(log2FC2),
                      sign(log2FC2) * qnorm(p2 / 2, lower.tail = FALSE),
                      NA_real_),
          
          # meta Z
          z_meta = dplyr::case_when(
            !is.na(z1) & !is.na(z2) ~ (z1 + z2) / sqrt(2),
            !is.na(z1) &  is.na(z2) ~ z1,
            is.na(z1) & !is.na(z2) ~ z2,
            TRUE ~ NA_real_
          )
        ) %>%
        # 删除 NA 和非有限值，防止 fgsea 报错
        dplyr::filter(!is.na(z_meta) & is.finite(z_meta))
      # =======================================================
      
      ranks_vec <- merged$z_meta
      names(ranks_vec) <- merged[[rlang::as_string(gene_sym)]]
      ranks_vec <- sort(ranks_vec, decreasing = TRUE)
      message("  [", ct, "] 使用 meta Z 排序（MAST + wilcoxon），基因数：", length(ranks_vec))
      
    } else {
      # 退回：原来单一方法的 log2FC 排序
      ranks <- sub1 %>%
        dplyr::select(!!gene_sym, log2FC) %>%
        dplyr::filter(!is.na(log2FC) & !is.na(!!gene_sym)) %>%
        dplyr::group_by(!!gene_sym) %>%
        dplyr::slice_max(order_by = abs(log2FC), n = 1, with_ties = FALSE) %>%
        dplyr::ungroup()
      
      ranks_vec <- ranks$log2FC
      names(ranks_vec) <- ranks[[rlang::as_string(gene_sym)]]
      ranks_vec <- sort(ranks_vec, decreasing = TRUE)
      message("  [", ct, "] 仅使用 MAST log2FC 排序，基因数：", length(ranks_vec))
    }
    # ------------------ 排序向量构建结束 ------------------
    
    # 保险：这里再检查一次（正常情况下不会触发）
    if (any(!is.finite(ranks_vec))) {
      bad_n <- sum(!is.finite(ranks_vec))
      message("  [警告] ", ct, " 中存在 ", bad_n, " 个非有限 stats，已移除。")
      ranks_vec <- ranks_vec[is.finite(ranks_vec)]
    }
    
    if (length(ranks_vec) < minSize) {
      message("  [跳过] ", ct, "：可用基因不足（", length(ranks_vec), "）")
      next
    }
    
    fg <- suppressWarnings(
      fgsea(pathways = pathways,
            stats = ranks_vec,
            minSize = minSize,
            maxSize = maxSize,
            nperm = nperm)
    )
    fg <- as.data.frame(fg)
    
    rssd <- paste0(RegionSubclass0, "_", Sex0, "_GSEA")
    
    res <- format_fgsea(
      fg = fg,
      ranks_len = length(ranks_vec),
      Sex = Sex0,
      Subclass = Subclass0,
      Region = Region0,
      RegionSubclass = RegionSubclass0,
      Neurotransmitter = Neuro0,
      region_Subclass_sex_direction = rssd
    )
    res$Region_Subclass_sex_direction <- paste0(res$Region.Subclass,"_",res$Sex,"_",res$Direction)
    
    ct_safe <- str_replace_all(ct, "[^A-Za-z0-9_\\-\\.]", "_")
    p <- gsub(".csv","",basename(csv_file))
    dir.create(paste0(out_dir,"/",p), showWarnings = FALSE, recursive = TRUE)
    
    per_ct_out_up <- file.path(out_dir, paste0(p, "/GO_gsea_", ct_safe, "_Up.csv"))
    per_ct_out_down <- file.path(out_dir, paste0(p, "/GO_gsea_", ct_safe, "_Down.csv"))
    
    res <- res[res$p_value <0.05, ]
    res_up <- res[res$NES > 0, ]
    res_down <- res[res$NES < 0, ]
    
    fwrite(res_up, per_ct_out_up)
    fwrite(res_down, per_ct_out_down)
    all_results[[ct]] <- res
    
    message("  [完成] ", ct, " -> ", per_ct_out_up, " (", nrow(res_up), " 行)")
    message("  [完成] ", ct, " -> ", per_ct_out_down, " (", nrow(res_down), " 行)")
  }
  
  p <- gsub(".csv","",basename(csv_file))
  merged_out <- file.path(out_dir, paste0(p, "/GO_gsea.csv"))
  merged <- data.table::rbindlist(all_results, use.names = TRUE, fill = TRUE)
  fwrite(merged, merged_out)
  message("[4/6] 合并结果：", merged_out, " (", nrow(merged), " 行)")
}






# 你的 MAST 结果（按需增减/替换通配符；可放多个模式）
input_globs <- c(
  "/data2st1/junyi/output/atac1112/dar/celltype.L2/MASTNG_*_annotated.csv"
)

# 输出目录
out_dir <- "/data2st1/junyi/output/atac1112/dar/gsea/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ========== 扫描所有输入并运行 ==========
expand_glob <- function(glob) {
  files <- Sys.glob(glob)
  files[file.exists(files)]
}
all_csv <- unique(unlist(lapply(input_globs, expand_glob)))

# all_csv <- unique(
#   c(
#     unlist(lapply("/data1st2/hannan_25/data/RNAexp/snRNA/from_MK/mast_all_ngeneson_saturation_downsampled_ratio/all_ngeneson_saturation_downsampled_ratio/*csv", expand_glob))
#     #unlist(lapply("/data1st2/hannan_25/data/RNAexp/snRNA/from_MK/wilcoxon_all_sample_downsampled_ratio/all_sample_downsampled_ratio/*csv", expand_glob))
#   )
# )
# all_csv="/data1st2/hannan_25/data/RNAexp/snRNA/from_MK/mast_all_ngeneson_saturation_downsampled_ratio/all_ngeneson_saturation_downsampled_ratio/adata_SUS_4v4_downsampled_ratio_mast.csv"
# all_csv="/data1st2/hannan_25/data/RNAexp/snRNA/from_MK/mast_all_ngeneson_saturation_downsampled_ratio/all_ngeneson_saturation_downsampled_ratio/adata_CSDS_4v3_downsampled_ratio_mast.csv"
message("[5/6] 发现比较文件数：", length(all_csv))
if (length(all_csv) == 0) stop("未找到任何 CSV 输入，请检查 input_globs。")

###并行计算###
library(parallel)
ncore <- max(1, detectCores() - 1)

Sys.setenv(MKL_NUM_THREADS="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
set.seed(2025)
RNGkind("L'Ecuyer-CMRG")

mclapply(all_csv, function(f) {
  out_dir0 <- file.path(out_dir, gsub("([/]all.*)|(^.*from_MK[/])","",f))
  dir.create(out_dir0, showWarnings = FALSE, recursive = TRUE)
  run_one_file(f, pathways, out_dir0)
}, mc.cores = ncore)
message("[6/6] 全部完成（并行 ", ncore, " 核）。")

