# 在绘图前设置随机种子
set.seed(123)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggforce)

readGO<-function(f='data/GO_orderT_p0.01_nlogP_noreplication.csv'){
  go_results<- read.csv(f)
  
  go_results$region <- gsub("-", "_", go_results$region)
  go_results$supregion<-go_results$region
  go_results$supregion[go_results$region %in% c("OPC_Oligo","Astro_Epen","Immune","Vascular")]<-'NN'
  
  #先定义变量为因子
  region_order <- c("PFC", "Isocortex", "HPF", "TH","HY","MB","AMY","STR","OPC_Oligo","Astro_Epen","Immune","Vascular")
  sex_order <- c("Male", "Female")
  table(go_results$Neurotransmitter)
  Neurotransmitter_order <- c("Glut", "GABA", "Chol","Dopa")
  
  go_results$region <- factor(go_results$region, levels = region_order)
  go_results$sex <- factor(go_results$sex, levels = sex_order)
  go_results$Neurotransmitter <- factor(go_results$Neurotransmitter, levels = Neurotransmitter_order)
  go_results$supregion <- factor(go_results$supregion, levels = c("PFC", "Isocortex", "HPF", "TH","HY","MB","AMY","STR","NN"))
  
  # 使用 order() 按 region、sex、Neurotransmitter 排序数据框
  go_results <- go_results[order(go_results$region, go_results$sex, go_results$Neurotransmitter,go_results$celltype), ]
  
  go_results$term_name<-paste(go_results$source,go_results$term_name,sep = '_')
  return(go_results)
}


# Function ----------------------------------------------------------------
goheatmap<-function(go,cluster_rows=TRUE,cutoff=10,title="_GO_BP_MF_heatmap_cutoff10_rotated.pdf",height = 40, width = 10){
  heatmap_data_BP_MF <- go %>%
    select(sample, term_name, nlog10_p_val_adj) %>%
    pivot_wider(names_from = term_name, values_from = nlog10_p_val_adj)
  
  mat_BP_MF <- as.data.frame(heatmap_data_BP_MF)
  rownames(mat_BP_MF) <- mat_BP_MF$sample
  mat_BP_MF$sample <- NULL
  mat_BP_MF <- as.matrix(mat_BP_MF)
  mat_BP_MF[is.na(mat_BP_MF)] <- 0
  mat_BP_MF_clipped <- pmax(pmin(mat_BP_MF, cutoff), -cutoff)
  mat_BP_MF_clipped <- pmax(pmin(mat_BP_MF, cutoff), -cutoff)
  
  
  set.seed(111)
  
  annotation_row <- go %>%
    select(sample, region, sex) %>%
    distinct()
  
  rownames(annotation_row) <- annotation_row$sample
  annotation_row$sample <- NULL
  
  #定义热图颜色函数 (蓝-白-红)
  col_fun <- colorRamp2(c(-cutoff, 0, cutoff), c("#3498db", "#FFFFFF", "#d62728"))
  
  #定义注释颜色
  ann_colors <- list(
    region = c(
      AMY = "#1f77b4",
      Isocortex = "#ff7f0e",
      HPF = "#009E73",
      HY = "#d62728",
      MB = "#9467bd",
      PFC = "#8c564b",
      STR = "#e377c2",
      TH = "#bcbd22",            
      Astro_Epen = "#17becf",    
      OPC_Oligo = "#7f7f7f", 
      "Astro-Epen" = "#17becf",    
      "OPC-Oligo" = "#7f7f7f", 
      Immune = "#ff9896",        
      Vascular = "#c5b0d5",
      NN= "#494949ff"  # Non-neuron cells       
    ),
    sex = c(
      Male = "#0080FF",
      Female = "#E800E8",
      M = "#0080FF",
      F = "#E800E8"
    )
  )
  
  #创建行注释对象 RowAnnotation
  row_ha <- rowAnnotation(
    region = annotation_row$region,
    sex = annotation_row$sex,
    col = ann_colors,
    annotation_name_side = "top"  # 注释名放在顶端
  )
  
  # 转置矩阵
  mat_t <- t(mat_BP_MF_clipped)
  
  #转置注释
  annotation_col <- annotation_row
  rownames(annotation_col) <- rownames(annotation_row)
  
  # 创建列注释（原来是行注释数据，转成列对应）
  col_ha <- HeatmapAnnotation(
    region = annotation_col$region,
    sex = annotation_col$sex,
    col = ann_colors,
    annotation_name_side = "right"
  )
  print(dim(mat_t))
  pdf(paste0("figures/",title), height = dim(mat_t)[1]/15+8, width = dim(mat_t)[2]/15+8)
  
  p1<-Heatmap(mat_t,
              name = "nlog10_p_val_adj",
              col = col_fun,
              heatmap_legend_param = list(
                at = c(-cutoff, -cutoff/2, 0, cutoff/2, cutoff),
                labels = c(-cutoff, -cutoff/2, 0, cutoff/2, cutoff)
              ),
              show_row_names = TRUE,        # 显示旋转后行名（原来列名 GO term）
              show_column_names = TRUE,     # 显示旋转后列名（原来行名 sample）
              cluster_rows = cluster_rows,
              cluster_columns = TRUE,
              top_annotation = col_ha,      # 注释条放上方，且对应转置后列
              row_names_side = "right",     # 原来列名 GO term 在热图右侧
              column_names_side = "bottom", # 原来行名 sample 在热图下方
              row_names_gp = gpar(fontsize = 3),
              column_names_gp = gpar(fontsize = 3),  
              use_raster = FALSE
  )
  #print(p1)
  ht_temp_drawn <- draw(p1)
  row_order <- row_order(ht_temp_drawn)
  ordered_row_names <- rownames(mat_t)[row_order]
  dev.off()
  return(ordered_row_names)
}

goheatmap1<-function(go,cluster_rows=TRUE,cutoff=10,title="_GO_BP_MF_heatmap_cutoff10_rotated.pdf",height = 40, width = 10){
  heatmap_data_BP_MF <- go %>%
    select(sample, term_name, nlog10_p_val_adj) %>%
    pivot_wider(names_from = term_name, values_from = nlog10_p_val_adj)
  
  mat_BP_MF <- as.data.frame(heatmap_data_BP_MF)
  rownames(mat_BP_MF) <- mat_BP_MF$sample
  mat_BP_MF$sample <- NULL
  mat_BP_MF <- as.matrix(mat_BP_MF)
  mat_BP_MF[is.na(mat_BP_MF)] <- 0
  mat_BP_MF_clipped <- pmax(pmin(mat_BP_MF, cutoff), -cutoff)
  mat_BP_MF_clipped <- pmax(pmin(mat_BP_MF, cutoff), -cutoff)
  
  col <- go %>% distinct(sample,Neurotransmitter)
  
  set.seed(111)
  
  annotation_row <- go %>%
    select(sample, region, sex) %>%
    distinct()
  
  rownames(annotation_row) <- annotation_row$sample
  annotation_row$sample <- NULL
  
  #定义热图颜色函数 (蓝-白-红)
  col_fun <- colorRamp2(c(-cutoff, 0, cutoff), c("#3498db", "#FFFFFF", "#d62728"))
  
  #定义注释颜色
  ann_colors <- list(
    region = c(
      AMY = "#1f77b4",
      Isocortex = "#ff7f0e",
      HPF = "#009E73",
      HY = "#d62728",
      MB = "#9467bd",
      PFC = "#8c564b",
      STR = "#e377c2",
      TH = "#bcbd22",            
      Astro_Epen = "#17becf",    
      OPC_Oligo = "#7f7f7f", 
      "Astro-Epen" = "#17becf",    
      "OPC-Oligo" = "#7f7f7f", 
      Immune = "#ff9896",        
      Vascular = "#c5b0d5",
      NN= "#494949ff"  # Non-neuron cells       
       
    ),
    sex = c(
      Male = "#0080FF",
      Female = "#E800E8",
      M = "#0080FF",
      F = "#E800E8"
    )
  )
  
  #创建行注释对象 RowAnnotation
  row_ha <- rowAnnotation(
    region = annotation_row$region,
    sex = annotation_row$sex,
    col = ann_colors,
    annotation_name_side = "top"  # 注释名放在顶端
  )
  
  # 转置矩阵
  mat_t <- t(mat_BP_MF_clipped)
  
  #转置注释
  annotation_col <- annotation_row
  rownames(annotation_col) <- rownames(annotation_row)
  
  # 创建列注释（原来是行注释数据，转成列对应）
  col_ha <- HeatmapAnnotation(
    region = annotation_col$region,
    sex = annotation_col$sex,
    col = ann_colors,
    annotation_name_side = "right"
  )
  print(dim(mat_t))
  pdf(paste0("figures/",title), height = dim(mat_t)[1]/15+8, width = dim(mat_t)[2]/15+8)
  
  p1<-Heatmap(mat_t,
              name = "nlog10_p_val_adj",
              col = col_fun,
              heatmap_legend_param = list(
                at = c(-cutoff, -cutoff/2, 0, cutoff/2, cutoff),
                labels = c(-cutoff, -cutoff/2, 0, cutoff/2, cutoff)
              ),
              show_row_names = TRUE,        # 显示旋转后行名（原来列名 GO term）
              show_column_names = TRUE,     # 显示旋转后列名（原来行名 sample）
              cluster_rows = cluster_rows,
              cluster_columns = TRUE,
              column_split = col$Neurotransmitter,
              top_annotation = col_ha,      # 注释条放上方，且对应转置后列
              row_names_side = "right",     # 原来列名 GO term 在热图右侧
              column_names_side = "bottom", # 原来行名 sample 在热图下方
              row_names_gp = gpar(fontsize = 3),
              column_names_gp = gpar(fontsize = 3),  
              use_raster = FALSE
  )
  #print(p1)
  ht_temp_drawn <- draw(p1)
  row_order <- row_order(ht_temp_drawn)
  ordered_row_names <- rownames(mat_t)[row_order]
  dev.off()
  return(ordered_row_names)
}

goheatmap2<-function(go,row=NULL,col=NULL,cluster_rows=TRUE,cluster_columns = TRUE,cutoff=10,title="_GO_BP_MF_heatmap_cutoff10_rotated.pdf",height = 40, width = 10){
  heatmap_data_BP_MF <- go %>%
    select(sample, term_name, nlog10_p_val_adj) %>%
    pivot_wider(names_from = term_name, values_from = nlog10_p_val_adj)
  
  mat_BP_MF <- as.data.frame(heatmap_data_BP_MF)
  rownames(mat_BP_MF) <- mat_BP_MF$sample
  mat_BP_MF$sample <- NULL
  mat_BP_MF <- as.matrix(mat_BP_MF)
  mat_BP_MF[is.na(mat_BP_MF)] <- 0
  mat_BP_MF_clipped <- pmax(pmin(mat_BP_MF, cutoff), -cutoff)
  mat_BP_MF_clipped <- pmax(pmin(mat_BP_MF, cutoff), -cutoff)
  
  set.seed(111)
  
  annotation_row <- go %>%
    select(sample, region,Neurotransmitter, sex) %>%
    distinct()
  
  rownames(annotation_row) <- annotation_row$sample
  annotation_row$sample <- NULL
  
  if (is.null(col)) {
    annotation_row<-annotation_row[order(annotation_row$region,annotation_row$sex),]
    col_split<-NULL
  } else {
    annotation_row<-col
    col_split<-col$split
  }
  
  mat_BP_MF_clipped<-mat_BP_MF_clipped[rownames(annotation_row),]
  
  #定义热图颜色函数 (蓝-白-红)
  col_fun <- colorRamp2(c(-cutoff, 0, cutoff), c("#3498db", "#FFFFFF", "#d62728"))
  
  #定义注释颜色
  ann_colors <- list(
    neurotransmitter = c(
      'Glut'='#FFC000',
      'GABA'='#00B050',
      'Dopa'='#ff7f0e',
      'Chol'='#1f77b4',
      'Hist'='#aa40fc',
      'Sero'='#e377c2',
      'NN'='#8c564b'
    ),
    region = c(
      AMY = "#1f77b4",
      Isocortex = "#ff7f0e",
      HPF = "#009E73",
      HY = "#d62728",
      MB = "#9467bd",
      PFC = "#8c564b",
      STR = "#e377c2",
      TH = "#bcbd22",            
      Astro_Epen = "#17becf",    
      OPC_Oligo = "#7f7f7f", 
      "Astro-Epen" = "#17becf",    
      "OPC-Oligo" = "#7f7f7f", 
      Immune = "#ff9896",        
      Vascular = "#c5b0d5"       
    ),
    sex = c(
      Male = "#0080FF",
      Female = "#E800E8",
      M = "#0080FF",
      F = "#E800E8"
    ),
    # Create a named vector for modules with associated colors
    module1 = c(
      M1 = '#59F5FD',
      M2 = '#782170',
      M3 = '#D86DCD',
      M4 = '#E49EDD',
      M5 = '#F2CEEF',
      M6 = '#7E350E',
      M7 = '#BE5014',
      M8 = '#F1A983',
      M9 = '#FFC000',
      M10 = '#92D050',
      M11 = '#00B050',
      M12 = '#83E28E',
      M13 = '#215C98',
      M14 = '#4D93D9',
      M15 = '#747474',
      M16 = '#ADADAD',
      M17 = '#FFFF00'
    ),
    module = c(
      M1 = '#59F5FD',
      M2 = '#D86DCD',
      M3 = '#E49EDD',
      M4 = '#FFC000',
      M5 = '#00B050',
      M6 = '#92D050',
      M7 = '#83E28E',
      M8 = '#4D93D9',
      M9 = '#747474',
      M10 = '#ADADAD',
      M11 = '#7E350E',
      M12 = '#BE5014',
      M13 = '#FFFF00',
      M14 = '#6666FF',
      M15 = '#782170',
      other='white'
    )
  )
  
  if (is.null(row)) {
    print("The variable 'row' is NULL.")
    row_split=NULL
    row_ha=NULL
  } else {
    mat_BP_MF_clipped<-mat_BP_MF_clipped[,row$Go]
    row_split<-row$module
    #创建行注释对象 RowAnnotation
    row_ha <- rowAnnotation(
      module = row$module,
      col = ann_colors,
      annotation_name_side = "top"  # 注释名放在顶端
    )
  }
  
  # 转置矩阵
  mat_t <- t(mat_BP_MF_clipped)
  
  #转置注释
  annotation_col <- annotation_row
  rownames(annotation_col) <- rownames(annotation_row)
  
  # 创建列注释（原来是行注释数据，转成列对应）
  col_ha <- HeatmapAnnotation(
    region = annotation_row$region,
    sex = annotation_row$sex,
    neurotransmitter=annotation_row$Neurotransmitter,
    col = ann_colors,
    annotation_name_side = "right"  # 注释名放在顶端
  )
  
  print(dim(mat_t))
  pdf(paste0("figures/",title), height = dim(mat_t)[1]/15+8, width = dim(mat_t)[2]/15+8)
  
  p1<-Heatmap(mat_t,
              name = "nlog10_p_val_adj",
              col = col_fun,
              heatmap_legend_param = list(
                at = c(-cutoff, -cutoff/2, 0, cutoff/2, cutoff),
                labels = c(-cutoff, -cutoff/2, 0, cutoff/2, cutoff)
              ),
              show_row_names = TRUE,        # 显示旋转后行名（原来列名 GO term）
              show_column_names = TRUE,     # 显示旋转后列名（原来行名 sample）
              cluster_rows = cluster_rows,
              cluster_columns = cluster_columns,
              column_split = col_split,
              row_split = row_split,
              left_annotation = row_ha,
              top_annotation = col_ha,      # 注释条放上方，且对应转置后列
              row_names_side = "right",     # 原来列名 GO term 在热图右侧
              column_names_side = "bottom", # 原来行名 sample 在热图下方
              row_names_gp = gpar(fontsize = 3),
              column_names_gp = gpar(fontsize = 3),  
              use_raster = FALSE
  )
  #print(p1)
  ht_temp_drawn <- draw(p1)
  ro <- unlist(row_order(ht_temp_drawn))
  co <- unlist(column_order(ht_temp_drawn))
  #ordered_row_names <- rownames(mat_t)[row_order]
  dev.off()
  return(mat_t[ro,co])
}

create_stacked_bar_plot_no_group <- function(go_results, topnumber, title, term_order = NULL,region_order=c("PFC", "Isocortex", "HPF", "TH","HY","MB","AMY","STR","OPC_Oligo","Astro_Epen","Immune","Vascular"),
                                             height = 40, width = 10, facet=FALSE) {
  shared_terms_df <- go_results %>%
    group_by(sex, direction, term_name) %>%
    summarise(
      n_celltype = n_distinct(celltype),
      .groups = "drop"
    )

    shared_terms<-shared_terms_df
  
  # 计算脑区组成
  region_composition_df <- go_results %>%
    group_by(sex, direction, term_name, region) %>%
    summarise(n_cells_region = n_distinct(celltype), .groups = "drop") %>%
    left_join(shared_terms, by = c("sex", "direction", "term_name")) %>%
    mutate(proportion = n_cells_region / n_celltype)
  
  region_composition_sum <- region_composition_df %>%
    group_by(term_name, region) %>%
    summarise(n_cells_region = sum(n_cells_region), .groups = "drop")
  
  
  shared_terms_sum <- shared_terms_df %>%
    group_by(term_name) %>%
    summarise(n_celltype = max(n_celltype), .groups = "drop")
  
  # 计算比例
  region_composition_sum <- region_composition_sum %>%
    left_join(shared_terms_sum, by = "term_name") %>%
    mutate(proportion = n_cells_region / n_celltype)
  
  # 选topnumber个term_name
  print(dim(shared_terms_sum))
  top_terms <- shared_terms_sum %>%
    arrange(desc(n_celltype)) %>%
    head(topnumber) %>%
    pull(term_name)
  print(dim(top_terms))
  
  if (!is.null(term_order)) {
    print(length(term_order))
    term_order_filtered <- term_order[term_order %in% top_terms]
    print(length(term_order_filtered))
    term_order_final <- c(term_order_filtered, setdiff(top_terms, term_order_filtered))
  } else {
    term_order_final <- rev(top_terms)
  }
  
  
  plot_data <- region_composition_sum %>%
    filter(term_name %in% term_order_final) %>%
    mutate(
      term_name = factor(term_name, levels = term_order_final),
      region = factor(region, levels = region_order)
    ) %>%
    arrange(term_name, region)
  
  region_colors<-c(
      AMY = "#1f77b4",
      Isocortex = "#ff7f0e",
      HPF = "#009E73",
      HY = "#d62728",
      MB = "#9467bd",
      PFC = "#8c564b",
      STR = "#e377c2",
      TH = "#bcbd22",            
      Astro_Epen = "#17becf",    
      OPC_Oligo = "#7f7f7f", 
      "Astro-Epen" = "#17becf",    
      "OPC-Oligo" = "#7f7f7f", 
      Immune = "#ff9896",        
      Vascular = "#c5b0d5"       
    )
  p <- ggplot(plot_data, aes(y = term_name, x = n_cells_region, fill = region)) +
    geom_col(position = "stack", color = "white", size = 0.2) +
    scale_fill_manual(values = region_colors, name = "Brain Region") +
    scale_x_continuous(expand = c(0, 0)) +
    labs(title = title,
         y = "GO Terms",
         x = "Number of Cell Types") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.y = element_text(angle = 0, hjust = 1, size = 1),
      axis.text.x = element_text(size = 6),
      axis.title = element_text(size = 8),
      legend.position = "right",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1, "lines")  # 增加子图之间的距离
    )
  if(facet){
    p<-p+facet_wrap(~ region, nrow = 1,scales = "free_x")+
      scale_x_continuous(expand = c(0, 0), limits = c(0, max(plot_data$n_cells_region, na.rm = TRUE)))  # 设置 x 轴范围
  }
  ggsave(paste0("figures/",title), p,limitsize = FALSE, height = height/2000*length(unique(plot_data$term_name)), width = width)
  return(p)
}

create_stacked_bar_plot_updown <- function(tmp,ord,output="figures/go/BP_MF_selected_go_num.pdf"){
  gotmp <-  tmp %>%
    group_by(term_name,sex,region) %>%
    summarise(
      count=n(),
      count_UP = sum(direction == "up"),      # 统计 "UP" 的数量
      count_DOWN = sum(direction == "down"),   # 统计 "DOWN" 的数量
      .groups = "drop"
    )
  gotmp$term_name<-factor(gotmp$term_name, levels = rev(ord))
  gotmp$region<-factor(gotmp$region,ro)
  ggplot(gotmp, aes(x = term_name)) +
    geom_bar(aes(y = count_UP, fill = sex), position = "stack", stat = "identity") +
    geom_bar(aes(y = -count_DOWN, fill = sex), position = "stack", stat = "identity") +
    scale_y_continuous(expand = c(0, 0), limits = c(-30, 35),labels = abs)+
    #scale_y_continuous(labels = abs) +
    labs(y = "Number", x = "Region") +
    scale_fill_manual(values = c("#E800E8", "#0080FF"))+nature_theme+
    coord_flip() +
    geom_hline(yintercept = 0, color = "black", size = 0.5)+
    facet_wrap(~ region, nrow = 1,scales = "free_x") +theme(axis.text.y = element_blank(),panel.spacing = unit(1, "lines"))
  ggsave(
    output, 
    width = 85+25, height = 85, units = "mm", #dpi = 600, # 85 mm square + 25 mm for legend 
    bg = "transparent",
    limitsize = FALSE,
    device = cairo_pdf
  )
}
# go_results<-readGO('data/GO_orderT_p0.01_nlogP_noreplication.csv')
# # BP+MF cutoff 10 by region ---------------------------------------------------------
# go_results <- go_results[order(go_results$region, go_results$sex, go_results$Neurotransmitter,go_results$celltype), ]
# go_results_BP_MF = subset(go_results,source %in% c("GO:BP","GO:MF"))
# goheatmap(go_results_BP_MF[go_results_BP_MF$supregion=="AMY",],title = "AMY_GO_BP_MF_heatmap_cutoff10_rotated.pdf")
# for( r in c("PFC", "Isocortex", "HPF", "TH","HY","MB","AMY","STR","NN")){ 
#   goheatmap(go_results_BP_MF[go_results_BP_MF$supregion==r,],title = paste0(r,"_GO_BP_MF_heatmap_cutoff10_rotated.pdf"))
# }


# # CC cutoff 10 by region ------------------------------------------------------------
# go_results_CC = subset(go_results,source %in% c("GO:CC"))
# goheatmap(go_results_CC[go_results_CC$supregion=="AMY",],title = "AMY_GO_CC_heatmap_cutoff10_rotated.pdf")
# for( r in c("PFC", "Isocortex", "HPF", "TH","HY","MB","AMY","STR","NN")){ 
#   goheatmap(go_results_CC[go_results_CC$supregion==r,],title = paste0(r,"_GO_CC_heatmap_cutoff10_rotated.pdf"))
# }


# # Neuron BP+MF cutoff 8 and 10 ---------------------------------------------------
# go_results <- go_results[order(go_results$region, go_results$sex, go_results$Neurotransmitter,go_results$celltype), ]
# go_results_BP_MF = subset(go_results,source %in% c("GO:BP","GO:MF"))
# ord_mb8<-goheatmap(go_results_BP_MF[go_results_BP_MF$supregion!="NN",],cutoff=8,title = "Neuron_GO_BP_MF_heatmap_cutoff8_rotated.pdf",height = 100, width = 40)
# ord_mb10<-goheatmap(go_results_BP_MF[go_results_BP_MF$supregion!="NN",],cutoff=10,title = "Neuron_GO_BP_MF_heatmap_cutoff10_rotated.pdf",height = 100, width = 40)
# p_all <- create_stacked_bar_plot_no_group(go_results_BP_MF[go_results_BP_MF$supregion!="NN",],topnumber = 200000,title = "Stackedbar_BP_MF_cutoff10_ordered.pdf",term_order = rev(ord_mb10))
# p_all <- create_stacked_bar_plot_no_group(go_results_BP_MF[go_results_BP_MF$supregion!="NN",],topnumber = 200000,title = "Stackedbar_BP_MF_cutoff10_ordered_v1.pdf",term_order = rev(ord_mb10),height = 40, width = 10,facet = TRUE)

# ### common
# tmp<-go_results_BP_MF[go_results_BP_MF$supregion!="NN",]
# mb10<-goheatmap2(tmp,cutoff=10,title = "Neuron_GO_BP_MF_heatmap_cutoff10_rotated.pdf",height = 100, width = 40)
# mb10v0<-mb10[rowSums(mb10!=0)>10,colSums(mb10!=0)>10]
# mb10v0<-mb10v0[rowSums(mb10v0!=0)>10,colSums(mb10v0!=0)>10]
# tmp<-go_results_BP_MF[go_results_BP_MF$supregion!="NN",]
# tmp<- tmp[(tmp$term_name %in% rownames(mb10v0)) & (tmp$sample %in% colnames(mb10v0)),]
# mb10<-goheatmap2(tmp,cutoff=10,title = "Neuron_GO_BP_MF_heatmap_cutoff10_common_v1.pdf",height = 100, width = 40)
# tmp1<-data.frame(rownames(mb10v0))
# tmp2<-gb_selected[gb_selected$Go %in% rownames(mb10v0),]
# colnames(tmp1)<-c("Go")
# rownames(tmp1)<-tmp1$Go
# tmp1$Module<-NA
# tmp1[tmp2$Go,"Module"]<-tmp2$Module
# write.csv(tmp1,"selectedBP_MF_common_v1.csv")

# ###
# tmp1<-read.csv("selectedBP_MF_common_v2.csv")
# tmp1$Go<-factor(tmp1$Go,levels = unique(tmp1$Go))
# tmp1$module<-factor(tmp1$Module,levels = unique(tmp1$Module))
# tmp<-go_results_BP_MF[go_results_BP_MF$supregion!="NN",]
# tmp<- tmp[(tmp$term_name %in% rownames(mb10v0)) & (tmp$sample %in% colnames(mb10v0)),]
# tmp$term_name <- factor(tmp$term_name, tmp1$Go)
# tmp<- tmp[order(tmp$term_name),]
# mb10_c<-goheatmap2(tmp,row = tmp1,cluster_rows = FALSE,cutoff=10,title = "Neuron_GO_BP_MF_heatmap_cutoff10_common_v2.pdf",height = 100, width = 40)
# mb10_c<-goheatmap2(tmp,row = tmp1,cluster_rows = TRUE,cutoff=10,title = "Neuron_GO_BP_MF_heatmap_cutoff10_common_cluster_v2.pdf",height = 100, width = 40)

# tmp1<-read.csv("selectedBP_MF_common_v4.csv")
# tmp1$Go<-factor(tmp1$Go,levels = unique(tmp1$Go))
# tmp1$module<-factor(tmp1$Module,levels = unique(tmp1$Module))
# tmp<-go_results_BP_MF[go_results_BP_MF$supregion!="NN",]
# tmp<- tmp[(tmp$term_name %in% rownames(mb10v0)) & (tmp$sample %in% colnames(mb10v0)),]
# tmp$term_name <- factor(tmp$term_name, tmp1$Go)
# tmp<- tmp[order(tmp$term_name),]
# mb10_c<-goheatmap2(tmp,row = tmp1,cluster_rows = TRUE,cutoff=10,title = "Neuron_GO_BP_MF_heatmap_cutoff10_common_v3.pdf",height = 100, width = 40)
# tmp1$Go<-factor(tmp1$Go,levels = rownames(mb10_c))
# tmp1<- tmp1[order(tmp1$module,tmp1$Go),]
# write.csv(tmp1,'selectedBP_MF_common_reorder_v4.csv')

# tmp$term_name <- factor(tmp$term_name, tmp1$Go)
# tmp<- tmp[order(tmp$term_name),]
# tmp2 <- tmp %>% distinct(sample,Neurotransmitter,region, sex)
# rownames(tmp2)<-tmp2$sample
# tmp2$region<-factor(tmp2$region,levels = ro)#c("HPF","PFC","Isocortex","TH","AMY","MB","HY","STR"))
# tmp2<-tmp2[order(tmp2$region,tmp2$sex,tmp2$Neurotransmitter),]
# mb10_c<-goheatmap2(tmp,row = tmp1,col=tmp2,cluster_rows = FALSE,cluster_columns = FALSE,cutoff=10,title = "Neuron_GO_BP_MF_heatmap_cutoff10_common_v4.pdf",height = 100, width = 40)
# tmp2<-tmp2[order(tmp2$sex,tmp2$region,tmp2$Neurotransmitter),]
# mb10_c<-goheatmap2(tmp,row = tmp1,col=tmp2,cluster_rows = FALSE,cluster_columns = FALSE,cutoff=10,title = "Neuron_GO_BP_MF_heatmap_cutoff10_common_v5.pdf",height = 100, width = 40)
# tmp2<-tmp2[order(tmp2$Neurotransmitter,tmp2$sex,tmp2$region),]
# mb10_c<-goheatmap2(tmp,row = tmp1,col=tmp2,cluster_rows = FALSE,cluster_columns = FALSE,cutoff=10,title = "Neuron_GO_BP_MF_heatmap_cutoff10_common_v6.pdf",height = 100, width = 40)



# tmp3<-tmp1[tmp1$module %in% c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9"),]
# tmp3$Go<-factor(tmp3$Go,levels = unique(tmp3$Go))
# mb10_c<-goheatmap2(tmp[tmp$term_name %in% tmp3$Go,],row = tmp3,cluster_rows = TRUE,cutoff=10,title = "Neuron_GO_BP_MF_heatmap_cutoff10_common_v7.pdf",height = 100, width = 40)

# create_stacked_bar_plot_updown(tmp,rownames(mb10_c),output="figures/go/BP_MF_selected_go_num_common.pdf")


# ### specific
# tmp<-go_results_BP_MF[go_results_BP_MF$supregion!="NN",]
# tmp<- tmp[((tmp$term_name %in% rownames(mb10v0)))==FALSE,]
# mb10_s<-goheatmap2(tmp,cutoff=10,title = "Neuron_GO_BP_MF_heatmap_cutoff10_specific_v1.pdf",height = 100, width = 40)
# tmp1<-data.frame(rownames(mb10_s))
# tmp2<-gb_selected[gb_selected$Go %in% rownames(mb10_s),]
# colnames(tmp1)<-c("Go")
# rownames(tmp1)<-tmp1$Go
# tmp1$Module<-NA
# tmp1[tmp2$Go,"Module"]<-tmp2$Module
# write.csv(tmp1,"selectedBP_MF_specific_v1.csv")


# tmp1<-read.csv("selectedBP_MF_specific_v2.csv")
# tmp1$Go<-factor(tmp1$Go,levels = unique(tmp1$Go))
# tmp1$module<-factor(tmp1$Module,levels = unique(tmp1$Module))
# tmp<-go_results_BP_MF[go_results_BP_MF$supregion!="NN",]
# tmp<- tmp[((tmp$term_name %in% rownames(mb10v0)))==FALSE,]
# tmp$term_name <- factor(tmp$term_name, tmp1$Go)
# tmp<- tmp[order(tmp$term_name),]
# mb10_c<-goheatmap2(tmp,row = tmp1,cluster_rows = FALSE,cutoff=10,title = "Neuron_GO_BP_MF_heatmap_cutoff10_specific_v3.pdf",height = 100, width = 40)
# tmp2 <- tmp %>% distinct(sample,Neurotransmitter,region, sex)
# rownames(tmp2)<-tmp2$sample
# tmp2$region<-factor(tmp2$region,levels = ro)#c("HPF","PFC","Isocortex","TH","AMY","MB","HY","STR"))
# tmp2<-tmp2[order(tmp2$region,tmp2$sex,tmp2$Neurotransmitter),]
# mb10_c<-goheatmap2(tmp,row = tmp1,col=tmp2,cluster_rows = FALSE,cluster_columns = FALSE,cutoff=10,title = "Neuron_GO_BP_MF_heatmap_cutoff10_specific_v4.pdf",height = 100, width = 40)
# tmp2<-tmp2[order(tmp2$sex,tmp2$region,tmp2$Neurotransmitter),]
# mb10_c<-goheatmap2(tmp,row = tmp1,col=tmp2,cluster_rows = FALSE,cluster_columns = FALSE,cutoff=10,title = "Neuron_GO_BP_MF_heatmap_cutoff10_specific_v5.pdf",height = 100, width = 40)
# tmp2<-tmp2[order(tmp2$Neurotransmitter,tmp2$sex,tmp2$region),]
# mb10_c<-goheatmap2(tmp,row = tmp1,col=tmp2,cluster_rows = FALSE,cluster_columns = FALSE,cutoff=10,title = "Neuron_GO_BP_MF_heatmap_cutoff10_specific_v6.pdf",height = 100, width = 40)

# tmp3<-tmp1[tmp1$module %in% c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10"),]
# mb10_c<-goheatmap2(tmp[tmp$term_name %in% tmp3$Go,],row = tmp3,cluster_rows = TRUE,cutoff=10,title = "Neuron_GO_BP_MF_heatmap_cutoff10_common_v7.pdf",height = 100, width = 40)

# mb10v0<-mb10[rowSums(mb10>0)>10,colSums(mb10>0)>15]
# mb10v0<-mb10v0[rowSums(mb10v0>0)>10,colSums(mb10v0>0)>15]
# tmp<-go_results_BP_MF[go_results_BP_MF$supregion!="NN",]
# tmp<- tmp[(tmp$term_name %in% rownames(mb10v0)) & (tmp$sample %in% colnames(mb10v0)),]
# mb10<-goheatmap2(tmp,cutoff=10,title = "Neuron_GO_BP_MF_heatmap_cutoff10_rotated_v2.pdf",height = 100, width = 40)

# gotmp <-  tmp %>%
#   group_by(term_name,sex,region) %>%
#   summarise(
#     count=n(),
#     count_UP = sum(direction == "up"),      # 统计 "UP" 的数量
#     count_DOWN = sum(direction == "down"),   # 统计 "DOWN" 的数量
#     .groups = "drop"
#   )
# gotmp$term_name<-factor(gotmp$term_name, levels = rev(rownames(mb10_c)))
# gotmp$region<-factor(gotmp$region,ro)
# ggplot(gotmp, aes(x = term_name)) +
#   geom_bar(aes(y = count_UP, fill = sex), position = "stack", stat = "identity") +
#   geom_bar(aes(y = -count_DOWN, fill = sex), position = "stack", stat = "identity") +
#   scale_y_continuous(expand = c(0, 0), limits = c(-30, 35),labels = abs)+
#   #scale_y_continuous(labels = abs) +
#   labs(y = "Number", x = "Region") +
#   scale_fill_manual(values = c("#E800E8", "#0080FF"))+nature_theme+
#   coord_flip() +
#   geom_hline(yintercept = 0, color = "black", size = 0.5)+
#   facet_wrap(~ region, nrow = 1,scales = "free_x") +theme(axis.text.y = element_blank(),panel.spacing = unit(1, "lines"))
# ggsave(
#   "figures/go/BP_MF_selected_go_num.png", 
#   width = 85+25, height = 85, units = "mm", #dpi = 600, # 85 mm square + 25 mm for legend 
#   bg = "white",
#   limitsize = FALSE,
#   #device = cairo_pdf
# )



# gotmp <-  tmp %>%
#   group_by(term_name,sex,region) %>%
#   summarise(
#     count=n(),
#     count_UP = sum(direction == "up"),      # 统计 "UP" 的数量
#     count_DOWN = sum(direction == "down"),   # 统计 "DOWN" 的数量
#     .groups = "drop"
#   )
# gotmp$term_name<-factor(gotmp$term_name, levels = rev(rownames(mb10)))
# gotmp$region<-factor(gotmp$region,ro)
# ggplot(gotmp, aes(x = term_name)) +
#   geom_bar(aes(y = count_UP, fill = sex), position = "stack", stat = "identity") +
#   geom_bar(aes(y = -count_DOWN, fill = sex), position = "stack", stat = "identity") +
#   #scale_y_continuous(expand = c(0, 0), limits = c(-40, 40),labels = abs)+
#   scale_y_continuous(labels = abs) +
#   labs(y = "Number", x = "Region") +
#   scale_fill_manual(values = c("#E800E8", "#0080FF"))+nature_theme+
#   coord_flip() +
#   geom_hline(yintercept = 0, color = "black", size = 0.5)+
#   facet_wrap(~ region, nrow = 1,scales = "free_x") #+theme(axis.text.y = element_blank(),panel.spacing = unit(1, "lines"))
# ggsave(
#   "figures/go/BP_MF_selected_go_num_common.pdf", 
#   width = 85*10+25, height = 85*40, units = "mm", #dpi = 600, # 85 mm square + 25 mm for legend 
#   bg = "transparent",
#   limitsize = FALSE,
#   device = cairo_pdf
# )


# # Neuron BP only cutoff 8 and 10 ---------------------------------------------------
# go_results_BP = subset(go_results,source %in% c("GO:BP"))
# ord_bp10<-goheatmap(go_results_BP[go_results_BP$supregion!="NN",],cutoff=10,title = "Neuron_GO_BP_heatmap_cutoff10_rotated.pdf",height = 100, width = 40)
# ord_bp10<-goheatmap1(go_results_BP[go_results_BP$supregion!="NN",],cutoff=10,title = "Neuron_GO_BP_heatmap_cutoff10_rotated_v1.pdf",height = 100, width = 40)
# p_all <- create_stacked_bar_plot_no_group(go_results_BP[go_results_BP$supregion!="NN",],topnumber = 200000,title = "Stackedbar_BP_cutoff10_ordered.pdf",term_order = rev(ord_bp10))
# p_all <- create_stacked_bar_plot_no_group(go_results_BP[go_results_BP$supregion!="NN",],topnumber = 200000,title = "Stackedbar_BP_cutoff10_ordered_v1.pdf",term_order = rev(ord_bp10),height = 40, width = 10,facet = TRUE)
# ro<-c("HPF","AMY","MB","TH", "HY","STR","PFC", "Isocortex")#c("AMY", "Isocortex", "HY","MB","TH","STR","PFC", "HPF")
# p_all <- create_stacked_bar_plot_no_group(go_results_BP[go_results_BP$supregion!="NN",],topnumber = 200000,title = "Stackedbar_BP_cutoff10_ordered_v2.pdf",term_order = rev(ord_bp10),region_order = ro,height = 40, width = 10,facet = TRUE)


# # Neuron BP selected ---------------------------------------------------
# gb_selected<-read_excel("X:/huang yin/DEGs_127/go_organized/goselected.xlsx")
# gb_selected$Gosuf<-sapply(strsplit(gb_selected$Go, "_"), function(x) tail(x, 1))
# gb_selected$Gosuf<-gsub("S", "s",gsub("P", "p",gsub("−", "-", gb_selected$Gosuf)))
# gb_selected$Gosuf<-paste0("GO:BP_",gb_selected$Gosuf)
# gb_selected$Go[c(gsub("−", "-", gb_selected$Go) %in% go_results_BP$term_name) == FALSE]<-gb_selected$Gosuf[c(gsub("−", "-", gb_selected$Go) %in% go_results_BP$term_name) == FALSE]
# gb_selected$Go<-gsub("−", "-", gb_selected$Go)
# gb_selected <- gb_selected[match(unique(gb_selected$Go), gb_selected$Go),]

# go_results_BP <- subset(go_results,source %in% c("GO:BP"))
# tmp<-go_results_BP[go_results_BP$supregion!="NN",]
# tmp <- subset(tmp, term_name %in% gb_selected$Go)
# tmp$term_name<-factor(tmp$term_name,levels = gb_selected$Go)
# tmp<-tmp %>% arrange(term_name)

# ord_bp10<-goheatmap2(tmp,gb_selected,cutoff=10,title = "Neuron_GO_BP_selected_heatmap_cutoff10_cluster_v2.pdf",height = 100, width = 40)
# ord_bp10<-goheatmap2(tmp,gb_selected,cutoff=10,title = "Neuron_GO_BP_selected_heatmap_cutoff10_cluster_v3.pdf",height = 100, width = 40)
# gb_selected_ord<-gb_selected
# gb_selected_ord$Go<-factor(gb_selected_ord$Go,levels = ord_bp10)
# gb_selected_ord<-gb_selected_ord[order(gb_selected_ord$Go),]
# write.csv(gb_selected_ord[,c('Go','Module')],'selectedBP.csv')

# ord_bp10<-goheatmap(tmp,cutoff=10,title = "Neuron_GO_BP_selected_heatmap_cutoff10_cluster_v0.pdf",height = 100, width = 40)
# ord_bp10<-goheatmap(tmp,cutoff=10,cluster_rows = FALSE,title = "Neuron_GO_BP_selected_heatmap_cutoff10_nocluster_v0.pdf",height = 100, width = 40)

# ord_bp10<-goheatmap1(tmp,cutoff=10,title = "Neuron_GO_BP_selected_heatmap_cutoff10_cluster.pdf",height = 100, width = 40)
# gotmp <-  tmp %>%
#   group_by(term_name,sex,region) %>%
#   summarise(
#     count=n(),
#     count_UP = sum(direction == "up"),      # 统计 "UP" 的数量
#     count_DOWN = sum(direction == "down"),   # 统计 "DOWN" 的数量
#     .groups = "drop"
#   )
# gotmp$term_name<-factor(gotmp$term_name, levels = rev(ord_bp10))
# gotmp$region<-factor(gotmp$region,ro)
# ggplot(gotmp, aes(x = term_name)) +
#   geom_bar(aes(y = count_UP, fill = sex), position = "stack", stat = "identity") +
#   geom_bar(aes(y = -count_DOWN, fill = sex), position = "stack", stat = "identity") +
#   scale_y_continuous(expand = c(0, 0), limits = c(-40, 40),labels = abs)+
#   #scale_y_continuous(labels = abs) +
#   labs(y = "Number", x = "Region") +
#   scale_fill_manual(values = c("#E800E8", "#0080FF"))+nature_theme+
#   coord_flip() +
#   geom_hline(yintercept = 0, color = "black", size = 0.5)+
#   facet_wrap(~ region, nrow = 1,scales = "free_x") 
# ggsave(
#   "figures/go/BP_selected_go_num_by_term_gender_cluster.pdf", 
#   width = 85*10+25, height = 85*40, units = "mm", #dpi = 600, # 85 mm square + 25 mm for legend 
#   bg = "transparent",
#   limitsize = FALSE,
#   device = cairo_pdf
# )

# ord_bp10<-goheatmap1(tmp,cutoff=10,cluster_rows = FALSE,title = "Neuron_GO_BP_selected_heatmap_cutoff10_nocluster.pdf",height = 100, width = 40)
# gotmp <-  tmp %>%
#   group_by(term_name,sex,region) %>%
#   summarise(
#     count=n(),
#     count_UP = sum(direction == "up"),      # 统计 "UP" 的数量
#     count_DOWN = sum(direction == "down"),   # 统计 "DOWN" 的数量
#     .groups = "drop"
#   )
# gotmp$term_name<-factor(gotmp$term_name, levels = rev(ord_bp10))
# gotmp$region<-factor(gotmp$region,ro)
# ggplot(gotmp, aes(x = term_name)) +
#   geom_bar(aes(y = count_UP, fill = sex), position = "stack", stat = "identity") +
#   geom_bar(aes(y = -count_DOWN, fill = sex), position = "stack", stat = "identity") +
#   scale_y_continuous(expand = c(0, 0), limits = c(-40, 40),labels = abs)+
#   #scale_y_continuous(labels = abs) +
#   labs(y = "Number", x = "Region") +
#   scale_fill_manual(values = c("#E800E8", "#0080FF"))+nature_theme+
#   coord_flip() +
#   geom_hline(yintercept = 0, color = "black", size = 0.5)+
#   facet_wrap(~ region, nrow = 1,scales = "free_x") 
# ggsave(
#   "figures/go/BP_selected_go_num_by_term_gender_nocluster.pdf", 
#   width = 85*10+25, height = 85*40, units = "mm", #dpi = 600, # 85 mm square + 25 mm for legend 
#   bg = "transparent",
#   limitsize = FALSE,
#   device = cairo_pdf
# )


# # select BP new -----------------------------------------------------------
# gb_selected<-read.csv("selectedBP_v1.csv")
# gb_selected$module<-gb_selected$Module
# gb_selected$module<- factor(gb_selected$module,levels = c("M1","M2","M3","M4","M5","M6","M7","M8","M9","M10","M11","M12","M13","M14","M15","M16","M17"))
# go_results_BP <- subset(go_results,source %in% c("GO:BP"))
# tmp<-go_results_BP[go_results_BP$supregion!="NN",]
# tmp <- subset(tmp, term_name %in% gb_selected$Go)
# tmp$term_name<-factor(tmp$term_name,levels = gb_selected$Go)
# tmp<-tmp %>% arrange(term_name)

# ord_bp10<-goheatmap2(tmp,gb_selected,cluster_rows = FALSE,cutoff=10,title = "Neuron_GO_BP_selected_heatmap_cutoff10_cluster_v4.pdf",height = 100, width = 40)
# ord_bp10<-goheatmap2(tmp,gb_selected,cluster_rows = TRUE,cutoff=10,title = "Neuron_GO_BP_selected_heatmap_cutoff10_cluster_v5.pdf",height = 100, width = 40)
# gb_selected_ord<-gb_selected
# gb_selected_ord$Go<-factor(gb_selected_ord$Go,levels = ord_bp10)
# gb_selected_ord<-gb_selected_ord[order(gb_selected_ord$module,gb_selected_ord$Go),]
# write.csv(gb_selected_ord[,c('Go','Module')],'selectedBP_v5.csv')
# ord_bp10<-goheatmap2(tmp,gb_selected_ord,cluster_rows = FALSE,cutoff=10,title = "Neuron_GO_BP_selected_heatmap_cutoff10_cluster_v6.pdf",height = 100, width = 40)
# ord_bp10<-goheatmap2(tmp,gb_selected_ord,cluster_rows = FALSE,cluster_columns = FALSE,cutoff=10,title = "Neuron_GO_BP_selected_heatmap_cutoff10_cluster_v7.pdf",height = 100, width = 40)
# ord_bp10<-goheatmap2(tmp,gb_selected_ord,cluster_rows = TRUE,cluster_columns = TRUE,cutoff=10,title = "Neuron_GO_BP_selected_heatmap_cutoff10_cluster_v8.pdf",height = 100, width = 40)


# # Neuron CC only cutoff 8 and 10 -----------------------------------------------
# go_results_CC = subset(go_results,source %in% c("GO:CC"))
# ord_cc8<-goheatmap(go_results_CC[go_results_CC$supregion!="NN",],cutoff=8,title = "Neuron_GO_CC_heatmap_cutoff8_rotated.pdf",height = 40, width = 40)
# ord_cc10<-goheatmap(go_results_CC[go_results_CC$supregion!="NN",],cutoff=10,title = "Neuron_GO_CC_heatmap_cutoff10_rotated.pdf",height = 40, width = 40)
# write.csv(ord_cc10,'CC_ordered.csv')

# cc<-read.csv("CC_ordered.csv")
# cc$module<-paste0("M",sapply(cc$module, function(x) {as.numeric(sub("M", "", x)) + 1}))
# write.csv(cc,'selectedCC_v3.csv')
# tmp<-go_results_CC[go_results_CC$supregion!="NN",]
# tmp <- subset(tmp, term_name %in% cc$Go)
# tmp$term_name<-factor(tmp$term_name,levels = cc$Go)
# tmp<-tmp %>% arrange(term_name)
# ord_cc10<-goheatmap2(tmp,row = cc,cluster_rows = FALSE,cutoff=10,title = "Neuron_GO_CC_heatmap_cutoff10_rotated_v3.pdf",height = 40, width = 40)

# ord_cc10<-goheatmap2(tmp,row = cc,cutoff=10,title = "Neuron_GO_CC_heatmap_cutoff10_rotated_v2.pdf",height = 40, width = 40)
# cc$Go<-factor(cc$Go,levels = ord_cc10)
# cc<-cc[order(cc$Go),]
# write.csv(cc,'selectedCC.csv')

# ord_cc10<-goheatmap1(go_results_CC[go_results_CC$supregion!="NN",],cutoff=10,title = "Neuron_GO_CC_heatmap_cutoff10_rotated_v1.pdf",height = 40, width = 40)
# p_all <- create_stacked_bar_plot_no_group(go_results_CC[go_results_CC$supregion!="NN",],topnumber = 200000,title = "Stackedbar_CC_cutoff10_ordered.pdf",term_order = rev(ord_cc10),height = 10, width = 10)
# p_all <- create_stacked_bar_plot_no_group(go_results_CC[go_results_CC$supregion!="NN",],topnumber = 200000,title = "Stackedbar_CC_cutoff10_ordered_v1.pdf",term_order = rev(ord_cc10),height = 40, width = 10,facet = TRUE)


# # Neuron All cutoff 8 and 10 -----------------------------------------------------
# ord_all8<-goheatmap(go_results[go_results$supregion!="NN",],cutoff=8,title = "Neuron_GO_ALL_heatmap_cutoff8_rotated.pdf",height = 120, width = 40)
# ord_all10<-goheatmap1(go_results[go_results$supregion!="NN",],cutoff=10,title = "Neuron_GO_ALL_heatmap_cutoff10_rotated_v1.pdf",height = 120, width = 40)
# p_all <- create_stacked_bar_plot_no_group(go_results[go_results$supregion!="NN",],topnumber = 200000,title = "Stackedbar_ALL_cutoff10_ordered.pdf",term_order = rev(ord_all10))
# p_all <- create_stacked_bar_plot_no_group(go_results[go_results$supregion!="NN",],topnumber = 200000,title = "Stackedbar_ALL_cutoff10_ordered_v1.pdf",term_order = rev(ord_all10),height = 40, width = 10,facet = TRUE)
# p_all <- create_stacked_bar_plot_no_group(go_results[go_results$supregion!="NN",],topnumber = 200000,title = "Stackedbar_ALL_cutoff10_ordered_v2.pdf",term_order = rev(ord_all10),region_order = ro,height = 40, width = 10,facet = TRUE)



# # Kegg cutoff10 -----------------------------------------------------------
# kegg<-readGO('data/KEGG_orderT_p0.01_nlogP_noreplication.csv')
# kegg<-subset(kegg,source %in% c("KEGG"))
# kegg$supregion<-kegg$region
# kegg$supregion[kegg$region %in% c("OPC_Oligo","Astro_Epen","Immune","Vascular")]<-'NN'
# ord_all10<-goheatmap(kegg[kegg$supregion!="NN",],cutoff=10,title = "Neuron_KEGG_ALL_heatmap_cutoff10_rotated.pdf",height = 20, width = 20)
# p_all <- create_stacked_bar_plot_no_group(kegg[kegg$supregion!="NN",],topnumber = 200000,title = "Stackedbar_Neuron_KEGG_cutoff10_ordered.pdf",term_order = rev(ord_all10),height = 5, width = 5)
# p_all <- create_stacked_bar_plot_no_group(kegg[kegg$supregion!="NN",],topnumber = 200000,title = "Stackedbar_Neuron_KEGG_cutoff10_ordered_v1.pdf",term_order = rev(ord_all10),height = 20, width = 20,facet = TRUE)


# ord_all10<-goheatmap(kegg[kegg$supregion=="NN",],cutoff=10,title = "NN_ALL_heatmap_cutoff10_rotated.pdf",height = 20, width = 20)
# p_all <- create_stacked_bar_plot_no_group(kegg[kegg$supregion=="NN",],topnumber = 200000,title = "Stackedbar_NN_KEGG_cutoff10_ordered.pdf",term_order = rev(ord_all10),height = 5, width = 10)
# p_all <- create_stacked_bar_plot_no_group(kegg[kegg$supregion=="NN",],topnumber = 200000,title = "Stackedbar_NN_KEGG_cutoff10_ordered_v1.pdf",term_order = rev(ord_all10),height = 5, width = 10,facet = TRUE)
# for( r in c("PFC", "Isocortex", "HPF", "TH","HY","MB","AMY","STR")){ 
#   tmp<-goheatmap(kegg[kegg$supregion==r,],cutoff=10,title = paste0(r,"_KEGG_heatmap_cutoff10_rotated.pdf"),height = 20, width = 20)
# }



# # up and down -------------------------------------------------------------
# gotmp <-  go_results[go_results$supregion!="NN",] %>%
#   group_by(term_name,sex,region) %>%
#   summarise(
#     count=n(),
#     count_UP = sum(direction == "up"),      # 统计 "UP" 的数量
#     count_DOWN = sum(direction == "down"),   # 统计 "DOWN" 的数量
#     .groups = "drop"
#   )
# gotmp$term_name<-factor(gotmp$term_name, levels = rev(ord_all10))
# gotmp$region<-factor(gotmp$region,ro)
# ggplot(gotmp, aes(x = term_name)) +
#   geom_bar(aes(y = count_UP, fill = sex), position = "stack", stat = "identity") +
#   geom_bar(aes(y = -count_DOWN, fill = sex), position = "stack", stat = "identity") +
#   scale_y_continuous(expand = c(0, 0), limits = c(-40, 40),labels = abs)+
#   #scale_y_continuous(labels = abs) +
#   labs(y = "Number", x = "Region") +
#   scale_fill_manual(values = c("#E800E8", "#0080FF"))+nature_theme+
#   coord_flip() +
#   geom_hline(yintercept = 0, color = "black", size = 0.5)+
#   facet_wrap(~ region, nrow = 1,scales = "free_x") 
# ggsave(
#   "figures/go/ALL_go_num_by_term_gender.pdf", 
#   width = 85*10+25, height = 85*70, units = "mm", dpi = 600, # 85 mm square + 25 mm for legend 
#   bg = "transparent",
#   limitsize = FALSE,
#   device = cairo_pdf
# )

# gotmp <-  go_results_BP[go_results_BP$supregion!="NN",] %>%
#   group_by(term_name,sex,region) %>%
#   summarise(
#     count=n(),
#     count_UP = sum(direction == "up"),      # 统计 "UP" 的数量
#     count_DOWN = sum(direction == "down"),   # 统计 "DOWN" 的数量
#     .groups = "drop"
#   )
# gotmp$term_name<-factor(gotmp$term_name, levels = rev(ord_bp10))
# gotmp$region<-factor(gotmp$region,ro)
# ggplot(gotmp, aes(x = term_name)) +
#   geom_bar(aes(y = count_UP, fill = sex), position = "stack", stat = "identity") +
#   geom_bar(aes(y = -count_DOWN, fill = sex), position = "stack", stat = "identity") +
#   scale_y_continuous(expand = c(0, 0), limits = c(-40, 40),labels = abs)+
#   #scale_y_continuous(labels = abs) +
#   labs(y = "Number", x = "Region") +
#   scale_fill_manual(values = c("#E800E8", "#0080FF"))+nature_theme+
#   coord_flip() +
#   geom_hline(yintercept = 0, color = "black", size = 0.5)+
#   facet_wrap(~ region, nrow = 1,scales = "free_x") 
# ggsave(
#   "figures/go/BP_go_num_by_term_gender.png", 
#   width = 85*10+25, height = 85*40, units = "mm", #dpi = 600, # 85 mm square + 25 mm for legend 
#   bg = "white",
#   limitsize = FALSE,
#   #device = cairo_pdf
# )


# # region+neurotransmiter --------------------------------------------------
# go_results_BP = subset(go_results,source %in% c("GO:BP"))
# for( r in c("PFC", "Isocortex", "HPF", "TH","HY","MB","AMY","STR")){ 
#   for(g in c("Glut","GABA")){
#     tmp<-go_results_BP[go_results_BP$supregion==r & go_results_BP$Neurotransmitter==g,]
#     if(dim(tmp)[1]>0){
#       goheatmap(tmp,cutoff = 10,title = paste0("go1/",g,"_",r,"_GO_BP_heatmap_cutoff10_rotated.pdf"))
#     }
#   }
# }



# # Neurotransmiter ---------------------------------------------------------
# go_results_BP = subset(go_results,source %in% c("GO:BP"))
# ord_bp10<-goheatmap(go_results_BP[go_results_BP$supregion!="NN" & go_results_BP$Neurotransmitter=="Glut",],cutoff=10,title = "go1/Glut_GO_BP_heatmap_cutoff10_rotated.pdf",height = 40, width = 10)

# gotmp <-  go_results_BP[go_results_BP$supregion!="NN" & go_results_BP$Neurotransmitter=="Glut",] %>%
#   group_by(term_name,sex,region) %>%
#   summarise(
#     count=n(),
#     count_UP = sum(direction == "up"),      # 统计 "UP" 的数量
#     count_DOWN = sum(direction == "down"),   # 统计 "DOWN" 的数量
#     .groups = "drop"
#   )
# gotmp$term_name<-factor(gotmp$term_name, levels = rev(ord_bp10))
# gotmp$region<-factor(gotmp$region,ro)
# ggplot(gotmp, aes(x = term_name)) +
#   geom_bar(aes(y = count_UP, fill = sex), position = "stack", stat = "identity") +
#   geom_bar(aes(y = -count_DOWN, fill = sex), position = "stack", stat = "identity") +
#   scale_y_continuous(expand = c(0, 0), limits = c(-40, 40),labels = abs)+
#   #scale_y_continuous(labels = abs) +
#   labs(y = "Number", x = "Region") +
#   scale_fill_manual(values = c("#E800E8", "#0080FF"))+nature_theme+
#   coord_flip() +
#   geom_hline(yintercept = 0, color = "black", size = 0.5)+
#   facet_wrap(~ region, nrow = 1,scales = "free_x") 
# ggsave(
#   "figures/go1/Glut_GO_BP_num_by_term_gender.pdf", 
#   width = 85*10+25, height = 85*40, units = "mm", #dpi = 600, # 85 mm square + 25 mm for legend 
#   bg = "white",
#   limitsize = FALSE,
#   device = cairo_pdf
# )

# go_results_BP = subset(go_results,source %in% c("GO:BP"))
# ord_bp10<-goheatmap(go_results_BP[go_results_BP$supregion!="NN" & go_results_BP$Neurotransmitter=="GABA",],cutoff=10,title = "go1/GABA_GO_BP_heatmap_cutoff10_rotated.pdf",height = 40, width = 10)

# gotmp <-  go_results_BP[go_results_BP$supregion!="NN" & go_results_BP$Neurotransmitter=="GABA",] %>%
#   group_by(term_name,sex,region) %>%
#   summarise(
#     count=n(),
#     count_UP = sum(direction == "up"),      # 统计 "UP" 的数量
#     count_DOWN = sum(direction == "down"),   # 统计 "DOWN" 的数量
#     .groups = "drop"
#   )
# gotmp$term_name<-factor(gotmp$term_name, levels = rev(ord_bp10))
# gotmp$region<-factor(gotmp$region,ro)
# ggplot(gotmp, aes(x = term_name)) +
#   geom_bar(aes(y = count_UP, fill = sex), position = "stack", stat = "identity") +
#   geom_bar(aes(y = -count_DOWN, fill = sex), position = "stack", stat = "identity") +
#   scale_y_continuous(expand = c(0, 0), limits = c(-60, 60),labels = abs)+
#   #scale_y_continuous(labels = abs) +
#   labs(y = "Number", x = "Region") +
#   scale_fill_manual(values = c("#E800E8", "#0080FF"))+nature_theme+
#   coord_flip() +
#   geom_hline(yintercept = 0, color = "black", size = 0.5)+
#   facet_wrap(~ region, nrow = 1,scales = "free_x") 
# ggsave(
#   "figures/go1/GABA_GO_BP_num_by_term_gender.pdf", 
#   width = 85*10+25, height = 85*40, units = "mm", #dpi = 600, # 85 mm square + 25 mm for legend 
#   bg = "white",
#   limitsize = FALSE,
#   device = cairo_pdf
# )


# # Neurotransmiter ---------------------------------------------------------
# go_results_BP = subset(go_results,source %in% c("GO:BP"))
# ord_bp10<-goheatmap(go_results_BP[go_results_BP$supregion!="NN" & go_results_BP$Neurotransmitter=="Glut",],cutoff=10,title = "go1/Glut_GO_BP_heatmap_cutoff10_rotated.pdf",height = 40, width = 10)

# gotmp <-  go_results_BP[go_results_BP$supregion!="NN" & go_results_BP$Neurotransmitter=="Glut",] %>%
#   group_by(term_name,sex,region) %>%
#   summarise(
#     count=n(),
#     count_UP = sum(direction == "up"),      # 统计 "UP" 的数量
#     count_DOWN = sum(direction == "down"),   # 统计 "DOWN" 的数量
#     .groups = "drop"
#   )
# gotmp$term_name<-factor(gotmp$term_name, levels = rev(ord_bp10))
# gotmp$region<-factor(gotmp$region,ro)
# ggplot(gotmp, aes(x = term_name)) +
#   geom_bar(aes(y = count_UP, fill = sex), position = "stack", stat = "identity") +
#   geom_bar(aes(y = -count_DOWN, fill = sex), position = "stack", stat = "identity") +
#   scale_y_continuous(expand = c(0, 0), limits = c(-40, 40),labels = abs)+
#   #scale_y_continuous(labels = abs) +
#   labs(y = "Number", x = "Region") +
#   scale_fill_manual(values = c("#E800E8", "#0080FF"))+nature_theme+
#   coord_flip() +
#   geom_hline(yintercept = 0, color = "black", size = 0.5)+
#   facet_wrap(~ region, nrow = 1,scales = "free_x") 
# ggsave(
#   "figures/go1/Glut_GO_BP_num_by_term_gender.pdf", 
#   width = 85*10+25, height = 85*40, units = "mm", #dpi = 600, # 85 mm square + 25 mm for legend 
#   bg = "white",
#   limitsize = FALSE,
#   device = cairo_pdf
# )


# # GO sub total ------------------------------------------------------------
# library(readxl)
# go<-read_excel("X:/huang yin/to huangyin/GO BP+MF+CC-neurotransmitter.xlsx",sheet="Total")
# go$term_name<-paste(go$source,go$term_name,sep = '_')
# go$nlog10_p_val_adj<-as.numeric(go$nlog10_p_val_adj)
# go<-go %>% distinct()
# ord_bp10<-goheatmap(subset(go,source %in% c("GO:BP")),cutoff=10,title = "go2/GO_BP_heatmap_cutoff10_rotated.pdf",height = 40, width = 10)
# ord_bp10<-goheatmap(go,cutoff=10,title = "go2/GO_ALL_heatmap_cutoff10_rotated.pdf",height = 40, width = 10)


