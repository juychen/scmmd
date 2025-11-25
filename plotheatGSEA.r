set.seed(123)

library(tidyr)
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggforce)
library(rlang)
library(tools)
library(readxl)
library(grDevices)
library(grid)


read_GO_DEG <- function(f = 'data/GO_orderT_p0.01_nlogP_noreplication.csv', sex_col = sex, sheet = 1){
  # read GO or DEG results from a file
  ext <- tolower(file_ext(f))
  if (ext %in% c("xlsx", "xls")) {
    go_results <- read_excel(path = f, sheet = sheet)
    go_results <- as.data.frame(go_results, stringsAsFactors = FALSE)
  } else if (ext %in% c("csv")) {
    go_results <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop(sprintf("Unsupported file types: .%s (only .csv / .xlsx / .xls are supported)", ext))
  }
  
  # check column names
  if(any(names(go_results) == "" | is.na(names(go_results)))){
    warning("Found empty or NA column names, replacing with 'V' plus column index")
    bad_names <- which(names(go_results) == "" | is.na(names(go_results)))
    names(go_results)[bad_names] <- paste0("V", bad_names)
  }
  
  # define the order of factors
  region_order           <- c("HPF","AMY","MB","TH", "HY","STR","PFC", "iCTX",
                              "OPC_Oligo","Astro_Epen","Immune","Vascular")
  neurotransmitter_order <- c("Glut","GABA","Chol","Dopa","Hist","Sero","NN")
  sex_order           <- c("M","F")
  
  # ensure the gender format
  go_results <- go_results %>%
    mutate(across({{ sex_col }}, ~ {
      x <- str_trim(as.character(.x))
      case_when(
        str_to_lower(x) %in% c("male","m")   ~ "M",
        str_to_lower(x) %in% c("female","f") ~ "F",
        TRUE ~ NA_character_
      )
    }))
  
  # factorize the columns
  go_results <- go_results %>%
    mutate(
      Region           = factor(Region, levels = region_order),
      Neurotransmitter = factor(Neurotransmitter, levels = neurotransmitter_order),
      across({{ sex_col }}, ~ factor(.x, levels = sex_order))
    )
  
  # rank the cell types
  if (all(c("Region",{{ sex_col }},"Neurotransmitter","Subclass") %in% names(go_results))) {
    go_results <- go_results[order(go_results$Region,
                                   go_results[[sex_col]],
                                   go_results$Neurotransmitter,
                                   go_results$Subclass), ]
  }
  
  # merger the source and term_name columns if they exist
  if (all(c("source","term_name") %in% names(go_results))) {
    go_results$term_name <- paste(go_results$source, go_results$term_name, sep = "_")
  }
  
  # revise the gene_name column if it exists
  if ("Gene_name" %in% names(go_results)) {
    go_results <- go_results %>%
      mutate(
        Gene_name = Gene_name |>
          as.character() |>
          str_trim() |>
          str_replace_all('^[\'"]+|[\'"]+$', '')
      )
  }
  
  return(go_results)
}


heatmap <- function(
    go,
    row = NULL,
    col = NULL,
    sample_col = sample,       # celltype col
    term_col = term_name,    # gene col, go term col
    value_col = nlog10_p_val_adj, # value col, usually -log10(p_val_adj), avg_log2FC, etc.
    region_col = region,
    neurotransmitter_col = Neurotransmitter,
    sex_col = sex,
    #
    row_go_col = Go,           # "Go" col for go result, "gene_name" col for DEG result
    row_module_col = module,       # "module" col for go result, "ct_freq_total" col for DEG result
    col_sample_col = sample,       # col 里指向 sample 的列
    col_split_col = split,        # col 里用于列分块的列
    #
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    cutoff = NULL,                        # NULL for DEGs, 10 for GO
    title = "_GO_BP_MF_heatmap_cutoff10_rotated.pdf",
    height = 40, width = 10
){
  # prepare the data
  heatmap_data <- go %>%
    transmute(.s = {{ sample_col }},
              .t = {{ term_col }},
              .v = suppressWarnings(as.numeric(as.character({{ value_col }})))
    ) %>%
    pivot_wider(names_from = .t, values_from = .v)
  
  mat <- as.data.frame(heatmap_data)
  rownames(mat) <- mat$.s
  mat$.s <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  
  if (!is.null(cutoff)) {
    mat <- pmax(pmin(mat, cutoff), -cutoff)
  }
  
  set.seed(111)
  
  # col annotation
  annotation_row <- go %>%
    distinct(
      sample          = {{ sample_col }},
      region          = {{ region_col }},
      neurotransmitter= {{ neurotransmitter_col }},
      sex             = {{ sex_col }},
    )
  
  rownames(annotation_row) <- annotation_row$sample
  annotation_row$sample <- NULL
  
  if (is.null(col)) {
    annotation_row <- annotation_row[order(annotation_row$region, annotation_row$sex, annotation_row$neurotransmitter), , drop = FALSE]
    col_split <- NULL
  } else {
    ann_tmp <- col %>%
      mutate(sample = {{ col_sample_col }})
    rownames(ann_tmp) <- ann_tmp$sample
    col_split <- tryCatch(pull(ann_tmp, {{ col_split_col }}), error = function(e) NULL)
    ann_tmp$sample <- NULL
    annotation_row <- ann_tmp
  }
  
  # ensure the order of annotation_row matches the rownames of mat
  mat <- mat[rownames(annotation_row), , drop = FALSE]
  
  # color function for heatmap
  ann_colors <- list(
    neurotransmitter = c(
      'Glut'='#FFC000','GABA'='#00B050','Dopa'='#ff7f0e',
      'Chol'='#1f77b4','Hist'='#aa40fc','Sero'='#e377c2','NN'='#8c564b'
    ),
    region = c(
      AMY="#1f77b4", Isocortex="#ff7f0e", HPF="#009E73", HY="#d62728",
      MB="#9467bd", PFC="#8c564b", STR="#e377c2", TH="#bcbd22",
      Astro_Epen="#17becf", OPC_Oligo="#7f7f7f",
      "Astro-Epen"="#17becf", "OPC-Oligo"="#7f7f7f",
      Immune="#ff9896", Vascular="#c5b0d5"
    ),
    sex = c(M="#0080FF", F="#E800E8"),
    module = c(
      M1='#59F5FD', M2='#D86DCD', M3='#E49EDD', M4='#FFC000', M5='#00B050',
      M6='#92D050', M7='#83E28E', M8='#4D93D9', M9='#747474', M10='#ADADAD',
      M11='#7E350E', M12='#BE5014', M13='#FFFF00', M14='#6666FF', M15='#782170',
      other='white'
    )
  )
  
  ensure_palette <- function(present, base_map) {
    present <- unique(as.character(present))
    present <- present[!is.na(present)]
    
    # extract the base palette for the present categories
    base_sub <- base_map[intersect(names(base_map), present)]
    
    # if there are missing categories, generate distinct colors
    missing <- setdiff(present, names(base_map))
    if (length(missing) > 0) {
      extra <- grDevices::hcl(
        h = seq(15, 375, length.out = length(missing) + 1)[-1],
        c = 70, l = 65
      )
      names(extra) <- missing
      base_sub <- c(base_sub, extra[missing])
    }
    
    # reorder to match the order of appearance in data
    base_sub[match(present, names(base_sub))]
  }
  
  ann_colors$region <- ensure_palette(annotation_row$region, ann_colors$region)
  ann_colors$neurotransmitter <- ensure_palette(annotation_row$neurotransmitter, ann_colors$neurotransmitter)
  ann_colors$sex <- ensure_palette(annotation_row$sex, ann_colors$sex)
  
  # optional: check if the row and col have the same sample names
  row_split <- NULL
  row_ha <- NULL
  
  # row-side GO/module handling
  if (!is.null(row)) {
    nm_go  <- as_name(ensym(row_go_col))
    nm_mod <- as_name(ensym(row_module_col))
    
    # The GO column must be present
    if (!nm_go %in% names(row)) {
      stop(sprintf("`row` must contain column `%s`.", nm_go))
    }
    
    # Get GO order, coerce to character, and remove duplicates while preserving order
    go_order <- pull(row, {{ row_go_col }})
    go_order <- as.character(go_order)
    go_order <- go_order[!duplicated(go_order)]
    
    # Intersect with matrix colnames and keep the provided order
    keep_go <- go_order[go_order %in% colnames(mat)]
    if (length(keep_go) == 0L) {
      stop(sprintf("`%s` has no intersection with matrix column names. Check names and case.", nm_go))
    }
    dropped <- setdiff(go_order, keep_go)
    if (length(dropped) > 0L) {
      message(sprintf("Note: after applying the provided GO order, %d term(s) not found in the matrix were dropped.", length(dropped)))
    }
    
    # Filter and reorder matrix columns according to keep_go
    mat <- mat[, keep_go, drop = FALSE]
    
    # 5) If the module column exists, set split and row annotation; otherwise, skip
    if (nm_mod %in% names(row)) {
      mod_vec <- pull(row, {{ row_module_col }})
      # Align module vector to the filtered/ordered GO terms
      mod_vec <- mod_vec[match(keep_go, go_order)]
      
      row_split <- mod_vec
      ann_colors$module <- ensure_palette(row_split, ann_colors$module)
      row_ha <- rowAnnotation(
        module = mod_vec,
        col = ann_colors,
        annotation_name_side = "top",
        annotation_name_gp = gpar(col = NA)
      )
    } else {
      # No module column: no split and no row annotation
      row_split <- NULL
      row_ha <- NULL
    }
  }
  
  # transpose the matrix for heatmap
  mat_t <- t(mat)
  
  if (is.null(cutoff)) {
    # use the range of the matrix for color scaling
    rng <- range(mat_t, finite = TRUE)
    m <- max(abs(rng))
    # create a color function with blue for negative, white for zero, and red for positive
    at <- seq(-m, m, length.out = 5)
    col_fun <- colorRamp2(c(min(at), 0, max(at)),
                          c("#3498db", "#FFFFFF", "#d62728"))
    legend_param <- list(at = at, labels = format(at))
  } else {
    cutoff <- as.numeric(cutoff)
    at <- c(-cutoff, -cutoff/2, 0, cutoff/2, cutoff)
    col_fun <- colorRamp2(c(-cutoff, 0, cutoff),
                          c("#3498db", "#FFFFFF", "#d62728"))
    legend_param <- list(at = at, labels = at)
  }
  
  col_ha <- HeatmapAnnotation(
    region          = annotation_row$region,
    sex             = annotation_row$sex,
    neurotransmitter= annotation_row$neurotransmitter,
    col = ann_colors,
    annotation_name_side = "right"
  )
  
  # plotting
  pdf(title,
      height = nrow(mat_t)/15 + 8,
      width  = ncol(mat_t)/15 + 8)
  
  ht <- Heatmap(
    mat_t,
    name = as_name(ensym(value_col)) ,
    col  = col_fun,
    heatmap_legend_param = legend_param,
    show_row_names    = TRUE,
    show_column_names = TRUE,
    cluster_rows      = cluster_rows,
    cluster_columns   = cluster_columns,
    column_split      = col_split,
    row_split         = row_split,
    left_annotation   = row_ha,
    top_annotation    = col_ha,
    row_names_side    = "right",
    column_names_side = "bottom",
    row_names_gp      = gpar(fontsize = 3),
    column_names_gp   = gpar(fontsize = 3),
    use_raster        = FALSE
  )
  
  ht_drawn <- draw(ht)
  ro <- unlist(row_order(ht_drawn))
  co <- unlist(column_order(ht_drawn))
  dev.off()
  
  return(mat_t[ro, co, drop = FALSE])
}


heatmap1 <- function(
    go,
    row = NULL,
    col = NULL,
    sample_col = sample,            # celltype col
    term_col = term_name,           # gene col, go term col
    value_col = nlog10_p_val_adj,   # value col, usually -log10(p_val_adj), avg_log2FC, etc.
    region_col = region,
    neurotransmitter_col = NULL,
    sex_col = NULL,
    #
    row_go_col = Go,                # "Go" col for GO result, "gene_name" for DEG result
    row_module_col = module,        # "module" for GO result, "ct_freq_total" for DEG result
    col_sample_col = sample,        # in `col`, the sample column
    col_split_col = split,          # in `col`, the split column
    #
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    cutoff = NULL,                  # NULL for DEGs, e.g. 10 for GO
    title = "_GO_BP_MF_heatmap_cutoff10_rotated.pdf",
    height = 40, width = 10
){
  # ---------- 0) Prepare matrix ----------
  heatmap_data <- go %>%
    transmute(.s = {{ sample_col }},
              .t = {{ term_col }},
              .v = suppressWarnings(as.numeric(as.character({{ value_col }})))
    ) %>%
    pivot_wider(names_from = .t, values_from = .v)
  
  mat <- as.data.frame(heatmap_data)
  rownames(mat) <- mat$.s
  mat$.s <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  
  if (!is.null(cutoff)) {
    mat <- pmax(pmin(mat, cutoff), -cutoff)
  }
  
  set.seed(111)
  
  # ---------- 1) Column (sample) annotations ----------
  annotation_row <- go %>%
    distinct(
      sample          = {{ sample_col }},
      region          = {{ region_col }},
      neurotransmitter= {{ neurotransmitter_col }},
      sex             = {{ sex_col }}
    )
  
  rownames(annotation_row) <- annotation_row$sample
  annotation_row$sample <- NULL
  
  if (is.null(col)) {
    annotation_row <- annotation_row[order(annotation_row$region, annotation_row$sex, annotation_row$neurotransmitter), , drop = FALSE]
    col_split <- NULL
  } else {
    ann_tmp <- col %>% mutate(sample = {{ col_sample_col }})
    rownames(ann_tmp) <- ann_tmp$sample
    col_split <- tryCatch(pull(ann_tmp, {{ col_split_col }}), error = function(e) NULL)
    ann_tmp$sample <- NULL
    annotation_row <- ann_tmp
  }
  
  # ensure row order aligns
  mat <- mat[rownames(annotation_row), , drop = FALSE]
  
  # ---------- 2) Optional row-side GO/module handling ----------
  row_split <- NULL
  row_ha <- NULL
  
  if (!is.null(row)) {
    nm_go  <- as_name(ensym(row_go_col))
    nm_mod <- as_name(ensym(row_module_col))
    if (!nm_go %in% names(row)) stop(sprintf("`row` must contain column `%s`.", nm_go))
    
    go_order <- pull(row, {{ row_go_col }})
    go_order <- as.character(go_order)
    go_order <- go_order[!duplicated(go_order)]
    
    keep_go <- go_order[go_order %in% colnames(mat)]
    if (length(keep_go) == 0L) {
      stop(sprintf("`%s` has no intersection with matrix column names. Check names and case.", nm_go))
    }
    dropped <- setdiff(go_order, keep_go)
    if (length(dropped) > 0L) {
      message(sprintf("Note: %d term(s) in `%s` not found in matrix; dropped.", length(dropped), nm_go))
    }
    
    mat <- mat[, keep_go, drop = FALSE]
    
    if (nm_mod %in% names(row)) {
      mod_vec <- pull(row, {{ row_module_col }})
      mod_vec <- mod_vec[match(keep_go, go_order)]
      row_split <- mod_vec
      # row_ha is created later after dynamic palettes are computed
    }
  }
  
  # ---------- 3) Dynamic palettes (derive categories from data) ----------
  # Base preferred palettes
  base_region <- c(
    AMY="#1f77b4", Isocortex="#ff7f0e", HPF="#009E73", HY="#d62728",
    MB="#9467bd", PFC="#8c564b", STR="#e377c2", TH="#bcbd22",
    Astro_Epen="#17becf", OPC_Oligo="#7f7f7f",
    "Astro-Epen"="#17becf", "OPC-Oligo"="#7f7f7f",
    Immune="#ff9896", Vascular="#c5b0d5"
  )
  base_nt <- c(
    Glut='#FFC000', GABA='#00B050', Dopa='#ff7f0e',
    Chol='#1f77b4', Hist='#aa40fc', Sero='#e377c2', NN='#8c564b'
  )
  base_sex <- c(M="#0080FF", F="#E800E8")
  
  # Helper: subset base palette to present values and auto-color any unknowns
  make_palette <- function(values, base_map) {
    vals <- as.character(unique(values))
    vals <- vals[!is.na(vals)]
    base_sub <- base_map[intersect(names(base_map), vals)]
    missing <- setdiff(vals, names(base_map))
    if (length(missing) > 0) {
      # generate distinct fallback colors
      extra <- grDevices::hcl(
        h = seq(15, 375, length.out = length(missing) + 1)[-1],
        c = 70, l = 65
      )
      names(extra) <- missing
      base_sub <- c(base_sub, extra[missing])
    }
    # reorder to match the order of appearance in data
    base_sub[match(vals, names(base_sub))]
  }
  
  present_regions <- annotation_row$region
  present_nt      <- annotation_row$neurotransmitter
  present_sex     <- annotation_row$sex
  present_mods    <- if (!is.null(row_split)) row_split else character(0)
  
  dyn_colors <- list()
  if (length(na.omit(present_regions)) > 0) {
    dyn_colors$region <- make_palette(present_regions, base_region)
  }
  if (length(na.omit(present_nt)) > 0) {
    dyn_colors$neurotransmitter <- make_palette(present_nt, base_nt)
  }
  if (length(na.omit(present_sex)) > 0) {
    dyn_colors$sex <- make_palette(present_sex, base_sex)
  }
  if (length(na.omit(present_mods)) > 0) {
    # modules don't have a predefined map; create one on the fly
    dyn_colors$module <- make_palette(present_mods, base_map = c())
  }
  
  # Now that we have dynamic palettes, we can create row_ha if needed
  if (!is.null(row_split)) {
    row_ha <- rowAnnotation(
      module = row_split,
      col = dyn_colors,
      annotation_name_side = "top"
    )
  }
  
  # ---------- 4) Build color scale for the heatmap ----------
  mat_t <- t(mat)
  
  if (is.null(cutoff)) {
    rng <- range(mat_t, finite = TRUE)
    m <- max(abs(rng))
    at <- seq(-m, m, length.out = 5)
    col_fun <- circlize::colorRamp2(c(min(at), 0, max(at)),
                                    c("#3498db", "#FFFFFF", "#d62728"))
    legend_param <- list(at = at, labels = format(at))
  } else {
    cutoff <- as.numeric(cutoff)
    at <- c(-cutoff, -cutoff/2, 0, cutoff/2, cutoff)
    col_fun <- circlize::colorRamp2(c(-cutoff, 0, cutoff),
                                    c("#3498db", "#FFFFFF", "#d62728"))
    legend_param <- list(at = at, labels = at)
  }
  
  col_ha <- HeatmapAnnotation(
    region           = annotation_row$region,
    sex              = annotation_row$sex,
    neurotransmitter = annotation_row$neurotransmitter,
    col = dyn_colors,                      # << use dynamic palettes
    annotation_name_side = "right"
  )
  
  # ---------- 5) Plot ----------
  pdf(title,
      height = nrow(mat_t)/15 + 8,
      width  = ncol(mat_t)/15 + 8)
  
  ht <- Heatmap(
    mat_t,
    name = as_name(ensym(value_col)),
    col  = col_fun,
    heatmap_legend_param = legend_param,
    show_row_names    = TRUE,
    show_column_names = TRUE,
    cluster_rows      = cluster_rows,
    cluster_columns   = cluster_columns,
    column_split      = col_split,
    row_split         = row_split,
    left_annotation   = row_ha,
    top_annotation    = col_ha,
    row_names_side    = "right",
    column_names_side = "bottom",
    row_names_gp      = gpar(fontsize = 3),
    column_names_gp   = gpar(fontsize = 3),
    use_raster        = FALSE
  )
  
  ht_drawn <- draw(ht)
  ro <- unlist(row_order(ht_drawn))
  co <- unlist(column_order(ht_drawn))
  dev.off()
  
  return(mat_t[ro, co, drop = FALSE])
}


heatmap_final <- function(
    go,
    row = NULL,
    col = NULL,
    sample_col = "sample",
    term_col = "term_name",
    value_col = "value",
    value_name = "log2FC",
    region_col = NULL,
    neurotransmitter_col = NULL,
    sex_col = NULL,
    class_col = NULL,
    type_col = NULL,
    #
    row_go_col = "Gene_name",
    row_module_col = "Category",
    row_category_col = NULL,
    col_sample_col = "sample",
    col_order_priority = c("region","sex","neurotransmitter","class","type"),
    col_split_col = 'split',
    #
    cluster_rows = TRUE,
    cluster_row_slices = FALSE,
    cluster_columns = TRUE,
    cluster_column_slices = FALSE,
    distance = "euclidean",      # Distance method for clustering, default is "euclidean"
    method = "complete",         # Clustering method, default is "complete", "ward.D2"
    cutoff = NULL,
    logfc_cutoff = NULL,
    column_gap = unit(2, "mm"),  # Gap between column splits
    row_gap = unit(2, "mm"),     # Gap between row splits
    row_names_gp = 24,
    column_names_gp = 9,
    title = "_heatmap.pdf",
    height = 40, width = 10
){
  
  # ---------- 0) Prepare matrix ----------
  heatmap_data <- go %>%
    dplyr::select(
      .s = all_of(sample_col),
      .t = all_of(term_col),
      .v = all_of(value_col)
    ) %>%
    mutate(.v = suppressWarnings(as.numeric(as.character(.v)))) %>%
    pivot_wider(names_from = .t, values_from = .v)
  
  mat <- as.data.frame(heatmap_data)
  rownames(mat) <- mat$.s
  mat$.s <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  
  if (!is.null(cutoff)) {
    mat <- pmax(pmin(mat, cutoff), -cutoff)
  }
  
  set.seed(111)
  
  # ---------- 1) Column (sample) annotations ----------
  # Start with required sample column
  annotation_row <- go %>%
    dplyr::select(sample = all_of(sample_col))
  
  # Add region column if provided
  if (!is.null(region_col) && region_col %in% names(go)) {
    annotation_row <- annotation_row %>%
      mutate(region = go[[region_col]])
  }
  
  # Add optional columns only if they are provided
  if (!is.null(neurotransmitter_col) && neurotransmitter_col %in% names(go)) {
    annotation_row <- annotation_row %>%
      mutate(neurotransmitter = go[[neurotransmitter_col]])
  }
  
  if (!is.null(sex_col) && sex_col %in% names(go)) {
    annotation_row <- annotation_row %>%
      mutate(sex = go[[sex_col]])
  }
  
  if (!is.null(class_col) && class_col %in% names(go)) {
    annotation_row <- annotation_row %>%
      mutate(class = go[[class_col]])
  }
  
  if (!is.null(type_col) && type_col %in% names(go)) {
    annotation_row <- annotation_row %>%
      mutate(type = go[[type_col]])
  }
  
  # Remove duplicates
  annotation_row <- annotation_row %>% distinct()
  rownames(annotation_row) <- annotation_row$sample
  annotation_row$sample <- NULL
  
  # Enhanced col_split logic - FIXED VERSION
  col_split <- NULL
  if (!is.null(col)) {
    # Use custom column data
    ann_tmp <- col %>%
      dplyr::select(sample = all_of(col_sample_col))
    rownames(ann_tmp) <- ann_tmp$sample
    
    # Check for split column in col data or go data
    if (!is.null(col_split_col)) {
      if (col_split_col %in% names(col)) {
        split_data <- col %>%
          dplyr::select(sample = all_of(col_sample_col), split_var = all_of(col_split_col)) %>%
          distinct()
        rownames(split_data) <- split_data$sample
        col_split <- split_data[rownames(ann_tmp), "split_var"]
      } else if (col_split_col %in% names(go)) {
        # Get split values from go data
        split_data <- go %>%
          dplyr::select(sample = all_of(sample_col), split_var = all_of(col_split_col)) %>%
          distinct()
        rownames(split_data) <- split_data$sample
        col_split <- split_data[rownames(ann_tmp), "split_var"]
      }
    }
    
    ann_tmp$sample <- NULL
    annotation_row <- ann_tmp
  } else {
    # Enhanced ordering logic using col_order_priority - DO THIS FIRST
    if (!is.null(col_order_priority) && length(col_order_priority) > 0) {
      # Use specified priority order, then add remaining columns
      priority_cols <- intersect(col_order_priority, names(annotation_row))
      remaining_cols <- setdiff(names(annotation_row), col_order_priority)
      order_cols <- c(priority_cols, remaining_cols)
    } else {
      # Default fallback order
      default_order <- c("region", "sex", "neurotransmitter", "class", "type")
      order_cols <- intersect(default_order, names(annotation_row))
      remaining_cols <- setdiff(names(annotation_row), default_order)
      order_cols <- c(order_cols, remaining_cols)
    }
    
    if (length(order_cols) > 0) {
      annotation_row <- annotation_row[do.call(order, annotation_row[order_cols]), , drop = FALSE]
    }
    
    # FIXED: Generate column split directly from sorted annotation_row
    if (!is.null(col_split_col) && col_split_col %in% names(annotation_row)) {
      # Use the sorted annotation_row data directly (maintains factor order)
      col_split <- annotation_row[[col_split_col]]
    } else if (!is.null(col_split_col) && col_split_col %in% names(go)) {
      # Fallback: regenerate from go data if column not in annotation_row
      split_data <- go %>%
        dplyr::select(sample = all_of(sample_col), split_var = all_of(col_split_col)) %>%
        distinct()
      rownames(split_data) <- split_data$sample
      col_split <- split_data[rownames(annotation_row), "split_var"]
    }
  }
  
  # ensure row order aligns
  mat <- mat[rownames(annotation_row), , drop = FALSE]
  
  # ---------- Base color palettes ----------
  base_region <- c(AMY="#ff7f0e", STR="#e377c2", PFC="#8c564b", iCTX="#1f77b4", MB="#9467bd", TH="#bcbd22",
                   HY="#d62728", HPF="#009E73")
  base_nt    <- c(Glut="#FFC000", GABA="#00B050", Dopa="#ff7f0e", Chol="#1f77b4",
                  Hist="#aa40fc", Sero="#e377c2", NN="#8c564b")
  base_class <- c(Glut="#FFC000", GABA="#00B050", Dopa="#ff7f0e", Chol="#1f77b4",
                  Hist="#aa40fc", Sero="#e377c2",
                  OPC="#8194AC", Oligo="#7030A0", Astrocyte="#17becf", Microglia="#ff9896",
                  Epen="#FF7F50", Vascular="#5F8BB2", Immune="#EF6347")
  base_sex   <- c(M="#55A0FB", F="#FFA0A0")
  base_type  <- c(N="#4CAF50", NN="#FF9800")  # Type colors (N=green, NN=orange)
  
  # Enhanced column ordering logic: predefined order vs clustering
  if (!cluster_columns) {
    # Convert to factors with predefined levels (base color order)
    if ("region" %in% names(annotation_row)) {
      # Convert to factor with strict base_region order
      annotation_row$region <- factor(annotation_row$region, levels = names(base_region))
    }
    if ("neurotransmitter" %in% names(annotation_row)) {
      annotation_row$neurotransmitter <- factor(annotation_row$neurotransmitter, levels = names(base_nt))
    }
    if ("sex" %in% names(annotation_row)) {
      annotation_row$sex <- factor(annotation_row$sex, levels = names(base_sex))
    }
    if ("class" %in% names(annotation_row)) {
      annotation_row$class <- factor(annotation_row$class, levels = names(base_class))
    }
    if ("type" %in% names(annotation_row)) {
      annotation_row$type <- factor(annotation_row$type, levels = names(base_type))
    }
    
    # Reorder samples using col_order_priority with predefined factor levels
    if (!is.null(col_order_priority) && length(col_order_priority) > 0) {
      priority_cols <- intersect(col_order_priority, names(annotation_row))
      remaining_cols <- setdiff(names(annotation_row), col_order_priority)
      order_cols <- c(priority_cols, remaining_cols)
    } else {
      default_order <- c("region", "sex", "neurotransmitter", "class", "type")
      order_cols <- intersect(default_order, names(annotation_row))
      remaining_cols <- setdiff(names(annotation_row), default_order)
      order_cols <- c(order_cols, remaining_cols)
    }
    
    if (length(order_cols) > 0) {
      # Create ordering arguments with explicit factor handling
      order_args <- list()
      for (col in order_cols) {
        if (col %in% names(annotation_row)) {
          if (is.factor(annotation_row[[col]])) {
            # For factors, use as.numeric to respect level order
            order_args[[col]] <- as.numeric(annotation_row[[col]])
          } else {
            order_args[[col]] <- annotation_row[[col]]
          }
        }
      }
      
      # Apply ordering using the prepared arguments
      if (length(order_args) > 0) {
        order_idx <- do.call(order, order_args)
        annotation_row <- annotation_row[order_idx, , drop = FALSE]
        mat <- mat[rownames(annotation_row), , drop = FALSE]
        
        # FIXED: Generate col_split directly from sorted annotation_row
        if (!is.null(col_split_col) && col_split_col %in% names(annotation_row)) {
          # Use the sorted annotation_row data directly (maintains factor order)
          col_split <- annotation_row[[col_split_col]]
        } else if (!is.null(col_split) && !is.null(col_split_col) && col_split_col %in% names(go)) {
          # Fallback: regenerate from go data if column not in annotation_row
          split_data <- go %>%
            dplyr::select(sample = all_of(sample_col), split_var = all_of(col_split_col)) %>%
            distinct()
          rownames(split_data) <- split_data$sample
          col_split <- split_data[rownames(annotation_row), "split_var"]
        }
      }
    }
    
  }
  
  # ---------- 2) Enhanced row-side GO/module handling with missing gene inclusion ----------
  row_split <- NULL
  row_ha <- NULL
  row_annotations <- list()
  
  if (!is.null(row)) {
    if (!row_go_col %in% names(row)) {
      stop(sprintf("`row` must contain column `%s`.", row_go_col))
    }
    
    go_order <- pull(row, all_of(row_go_col))
    go_order <- as.character(go_order)
    go_order <- go_order[!duplicated(go_order)]
    
    # Find genes present and missing in matrix
    keep_go <- go_order[go_order %in% colnames(mat)]
    missing_go <- setdiff(go_order, keep_go)
    
    if (length(missing_go) > 0L) {
      message(sprintf("Note: %d term(s) in `%s` not found in matrix; adding with 0 values: %s",
                      length(missing_go), row_go_col, paste(head(missing_go, 5), collapse = ", ")))
      
      # Create matrix for missing genes with 0 values
      missing_mat <- matrix(0, nrow = nrow(mat), ncol = length(missing_go))
      colnames(missing_mat) <- missing_go
      rownames(missing_mat) <- rownames(mat)
      
      # Combine existing matrix with missing genes
      mat <- cbind(mat[, keep_go, drop = FALSE], missing_mat)
    } else {
      mat <- mat[, keep_go, drop = FALSE]
    }
    
    # Reorder matrix columns to match go_order
    mat <- mat[, go_order, drop = FALSE]
    
    # Handle multiple row grouping columns
    if (!is.null(row_module_col) && row_module_col %in% names(row)) {
      mod_vec <- pull(row, all_of(row_module_col))
      mod_vec <- mod_vec[match(go_order, pull(row, all_of(row_go_col)))]
      row_split <- mod_vec
      row_annotations$module <- mod_vec
    }
    
    # Handle additional category column
    if (!is.null(row_category_col) && row_category_col %in% names(row)) {
      cat_vec <- pull(row, all_of(row_category_col))
      cat_vec <- cat_vec[match(go_order, pull(row, all_of(row_go_col)))]
      row_annotations$category <- cat_vec
    }
  }
  
  # ---------- 3) Dynamic palettes ----------
  make_palette <- function(values, base_map) {
    vals <- as.character(unique(values))
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return(NULL)
    
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
  
  # Extract present values (region is now optional)
  present_regions <- if ("region" %in% names(annotation_row)) annotation_row$region else NULL
  present_nt      <- if ("neurotransmitter" %in% names(annotation_row)) annotation_row$neurotransmitter else NULL
  present_sex     <- if ("sex" %in% names(annotation_row)) annotation_row$sex else NULL
  present_class   <- if ("class" %in% names(annotation_row)) annotation_row$class else NULL
  present_type    <- if ("type" %in% names(annotation_row)) annotation_row$type else NULL
  present_mods    <- if (!is.null(row_split)) row_split else character(0)
  present_cats    <- if ("category" %in% names(row_annotations)) row_annotations$category else character(0)
  
  dyn_colors <- list()
  # Region is now optional
  if (!is.null(present_regions) && length(na.omit(present_regions)) > 0) {
    dyn_colors$region <- make_palette(present_regions, base_region)
  }
  if (!is.null(present_nt) && length(na.omit(present_nt)) > 0) {
    dyn_colors$neurotransmitter <- make_palette(present_nt, base_nt)
  }
  if (!is.null(present_sex) && length(na.omit(present_sex)) > 0) {
    dyn_colors$sex <- make_palette(present_sex, base_sex)
  }
  if (!is.null(present_class) && length(na.omit(present_class)) > 0) {
    dyn_colors$class <- make_palette(present_class, base_class)
  }
  # Add type colors
  if (!is.null(present_type) && length(na.omit(present_type)) > 0) {
    dyn_colors$type <- make_palette(present_type, base_type)
  }
  if (length(na.omit(present_mods)) > 0) {
    dyn_colors$module <- make_palette(present_mods, base_map = c())
  }
  if (length(na.omit(present_cats)) > 0) {
    dyn_colors$category <- make_palette(present_cats, base_map = c())
  }
  
  # Define consistent legend parameters for all legends
  legend_title_gp <- gpar(fontsize = 16, fontface = "bold")
  legend_labels_gp <- gpar(fontsize = 14)
  legend_grid_height <- unit(5, "mm")
  legend_grid_width <- unit(5, "mm")
  
  # ---------- 4) Build color scale with consistent styling ----------
  mat_t <- t(mat)
  
  # Create row_ha with multiple annotations if needed
  if (length(row_annotations) > 0) {
    # Build row annotation colors list
    row_colors_list <- list()
    row_legend_param_list <- list()
    
    for (ann_name in names(row_annotations)) {
      if (ann_name %in% names(dyn_colors)) {
        row_colors_list[[ann_name]] <- dyn_colors[[ann_name]]
        row_legend_param_list[[ann_name]] <- list(
          title_gp = legend_title_gp,
          labels_gp = legend_labels_gp,
          grid_height = legend_grid_height,
          grid_width = legend_grid_width
        )
      }
    }
    
    # Create rowAnnotation with all available annotations
    row_ha_args <- c(
      row_annotations,
      list(
        col = row_colors_list,
        annotation_name_side = "top",
        annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
        show_legend = TRUE,
        gp = gpar(col = NA),
        simple_anno_size = unit(0.3, "cm"),
        annotation_legend_param = row_legend_param_list
      )
    )
    
    row_ha <- do.call(rowAnnotation, row_ha_args)
  }
  
  # Keep original logFC cutoff logic - no changes to user's settings
  
  if (is.null(cutoff)) {
    rng <- range(mat_t, finite = TRUE)
    m <- max(abs(rng))
    at <- seq(-m, m, length.out = 5)
    col_fun <- circlize::colorRamp2(c(min(at), 0, max(at)),
                                    c("blue", "#FFFFFF", "red"))
    legend_param <- list(
      at = at,
      labels = format(at, digits = 3),
      title = "logFC",
      title_gp = legend_title_gp,
      labels_gp = legend_labels_gp
    )
  } else {
    cutoff <- as.numeric(cutoff)
    at <- c(-cutoff, -cutoff/2, 0, cutoff/2, cutoff)
    col_fun <- circlize::colorRamp2(c(-cutoff, 0, cutoff),
                                    c("blue", "#FFFFFF", "red"))
    legend_param <- list(
      at = at,
      labels = at,
      title = value_name,
      title_gp = legend_title_gp,
      labels_gp = legend_labels_gp
    )
  }
  
  # Build HeatmapAnnotation with consistent styling and custom order
  annotation_list <- list()
  
  # Use col_order_priority to determine the order of annotations
  if (!is.null(col_order_priority) && length(col_order_priority) > 0) {
    priority_cols <- intersect(col_order_priority, names(annotation_row))
    remaining_cols <- setdiff(names(annotation_row), col_order_priority)
    annotation_order <- c(priority_cols, remaining_cols)
  } else {
    # Default order
    default_order <- c("region", "sex", "neurotransmitter", "class", "type")
    priority_cols <- intersect(default_order, names(annotation_row))
    remaining_cols <- setdiff(names(annotation_row), default_order)
    annotation_order <- c(priority_cols, remaining_cols)
  }
  
  # Build annotation_list in the specified order
  for (ann_name in annotation_order) {
    if (ann_name %in% names(annotation_row)) {
      annotation_list[[ann_name]] <- annotation_row[[ann_name]]
    }
  }
  
  # Build annotation_legend_param dynamically based on actual annotations
  annotation_legend_param <- list()
  for (ann_name in names(annotation_list)) {
    annotation_legend_param[[ann_name]] <- list(
      title_gp = legend_title_gp,
      labels_gp = legend_labels_gp,
      grid_height = legend_grid_height,
      grid_width = legend_grid_width
    )
  }
  
  # Create annotation with consistent legend styling and custom order
  col_ha <- do.call(HeatmapAnnotation, c(
    annotation_list,
    list(
      col = dyn_colors,
      annotation_name_side = "right",
      annotation_legend_param = annotation_legend_param,
      annotation_name_gp = gpar(fontsize = 16, fontface = "bold")
    )
  ))
  
  # ---------- 5) Plot with format detection ----------
  # Detect file format from title extension
  file_ext <- tolower(tools::file_ext(title))
  
  if (file_ext == "png") {
    # PNG output
    png(title,
        height = (nrow(mat_t)/15 + 8) * 300,  # Convert to pixels (assuming 300 DPI)
        width  = (ncol(mat_t)/15 + 8) * 300,
        res = 300)  # High resolution for PNG
  } else {
    # Default to PDF output
    pdf(title,
        height = nrow(mat_t)/15 + 8,
        width  = ncol(mat_t)/15 + 8)
  }
  
  # Create cell function for custom coloring with logfc_cutoff
  cell_fun <- NULL
  logfc_legend <- NULL
  
  if (!is.null(logfc_cutoff)) {
    logfc_cutoff <- as.numeric(logfc_cutoff)
    
    cell_fun <- function(j, i, x, y, width, height, fill) {
      value <- mat_t[i, j]
      if (abs(value) <= logfc_cutoff) {
        # Grey color for values within logfc_cutoff range
        grid.rect(x = x, y = y, width = width, height = height,
                  gp = gpar(col = NA, fill = "#D3D3D3"))
      } else {
        # Use default coloring for significant values
        grid.rect(x = x, y = y, width = width, height = height,
                  gp = gpar(col = NA, fill = fill))
      }
    }
    
    # Create logfc_legend with CONSISTENT styling matching all other legends
    logfc_legend <- Legend(
      labels = sprintf("Not detected or not a DEG"),
      legend_gp = gpar(fill = "#D3D3D3"),
      title = NULL,  # No title to keep it clean
      title_gp = legend_title_gp,
      labels_gp = legend_labels_gp,
      grid_height = legend_grid_height,
      grid_width = legend_grid_width
    )
  }
  
  # Enhanced clustering logic: within-group clustering when col_split is provided
  final_cluster_columns <- cluster_columns
  
  # FIXED: Create heatmap with proper column split settings and within-group clustering
  ht <- Heatmap(
    mat_t,
    name = value_col,
    col  = col_fun,
    cell_fun = cell_fun,  # Apply custom cell coloring if logfc_cutoff is set
    heatmap_legend_param = legend_param,
    show_row_names    = TRUE,
    show_column_names = TRUE,
    cluster_rows      = cluster_rows,
    cluster_row_slices = cluster_row_slices,
    cluster_columns   = final_cluster_columns,  # Use enhanced clustering logic
    cluster_column_slices = cluster_column_slices,
    clustering_distance_rows    = distance,    # Distance method for row clustering
    clustering_method_rows      = method,      # Clustering method for rows
    clustering_distance_columns = distance,    # Distance method for column clustering
    clustering_method_columns   = method,      # Clustering method for columns
    column_split      = col_split,  # Enhanced split logic with within-group clustering
    row_split         = row_split,
    column_gap        = column_gap,  # Control column split gap width
    row_gap           = row_gap,     # Control row split gap width
    left_annotation   = row_ha,
    top_annotation    = col_ha,
    row_names_side    = "right",
    column_names_side = "bottom",
    row_names_gp      = gpar(fontsize = row_names_gp),
    column_names_gp   = gpar(fontsize = column_names_gp),
    use_raster        = FALSE,
    row_title         = NULL,
    column_title      = NULL,
    show_column_dend  = if (cluster_columns) TRUE else FALSE,
    column_title_gp   = gpar(fontsize = 0)  # FIXED: Set title font size to 0
  )
  
  # Draw heatmap with optional logfc_cutoff legend
  if (!is.null(logfc_legend)) {
    ht_drawn <- draw(ht, annotation_legend_list = list(logfc_legend))
  } else {
    ht_drawn <- draw(ht)
  }
  
  ro <- unlist(row_order(ht_drawn))
  co <- unlist(column_order(ht_drawn))
  dev.off()
  
  row_annhn=unique(row_data[,c(1,3)])
  rownames(row_annhn)=row_annhn$GO
  a=row_annhn[rownames(mat_t[ro,co]),2:1]
  colnames(a)=c("Group","GO/Cellgroup")
  aa=cbind(a,mat_t[ro,co])
  cc=cbind.data.frame("",cbind.data.frame(rownames(t(annotation_row[colnames(mat_t[ro,co]),])),t(annotation_row[colnames(mat_t[ro,co]),])))
  colnames(cc)[1:2]=colnames(aa)[1:2]
  f=rbind.data.frame(cc,colnames(aa),aa)
  colnames(f) <- NULL 
  # save matrix
  #if (grepl("nocluster", title, ignore.case = TRUE)) {
  #final_mat <- mat_t[ro, co, drop = FALSE]
  out_file <- paste0(file_path_sans_ext(title), "_plotdata.csv")
  write.csv(f, out_file, quote = TRUE,row.names = F)
  #}
  
  return(mat_t[ro, co, drop = FALSE])
}


heatmap_final3 <- function(
    go,
    row = NULL,
    col = NULL,
    sample_col = "sample",
    term_col = "term_name",
    value_col = "value",
    value_name = "logFC",
    region_col = NULL,
    neurotransmitter_col = NULL,
    sex_col = NULL,
    class_col = NULL,
    type_col = NULL,
    #
    row_go_col = "Gene_name",
    row_module_col = "Category",
    row_category_col = NULL,
    col_sample_col = "sample",
    col_order_priority = c("region","sex","neurotransmitter","class","type"),
    col_split_col = 'split',
    #
    cluster_rows = TRUE,
    cluster_row_slices = FALSE,
    cluster_columns = TRUE,
    cluster_column_slices = FALSE,
    distance = "euclidean",      # Distance method for clustering, default is "euclidean"
    method = "complete",         # Clustering method, default is "complete","ward.D2"
    cutoff = NULL,
    logfc_cutoff = NULL,
    column_gap = unit(2, "mm"),  # Gap between column splits
    row_gap = unit(2, "mm"),     # Gap between row splits
    row_names_gp = 24,
    column_names_gp = 9,
    title = "_heatmap.pdf",
    height = 40, width = 10
){
  
  # ---------- 0) Prepare matrix ----------
  heatmap_data <- go %>%
    dplyr::select(
      .s = all_of(sample_col),
      .t = all_of(term_col),
      .v = all_of(value_col)
    ) %>%
    mutate(.v = suppressWarnings(as.numeric(as.character(.v)))) %>%
    pivot_wider(names_from = .t, values_from = .v)
  
  mat <- as.data.frame(heatmap_data)
  rownames(mat) <- mat$.s
  mat$.s <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  
  if (!is.null(cutoff)) {
    mat <- pmax(pmin(mat, cutoff), -cutoff)
  }
  
  set.seed(111)
  
  # ---------- 1) Column (sample) annotations ----------
  # Start with required sample column
  annotation_row <- go %>%
    dplyr::select(sample = all_of(sample_col))
  
  # Add region column if provided
  if (!is.null(region_col) && region_col %in% names(go)) {
    annotation_row <- annotation_row %>%
      mutate(region = go[[region_col]])
  }
  
  # Add optional columns only if they are provided
  if (!is.null(neurotransmitter_col) && neurotransmitter_col %in% names(go)) {
    annotation_row <- annotation_row %>%
      mutate(neurotransmitter = go[[neurotransmitter_col]])
  }
  
  if (!is.null(sex_col) && sex_col %in% names(go)) {
    annotation_row <- annotation_row %>%
      mutate(sex = go[[sex_col]])
  }
  
  if (!is.null(class_col) && class_col %in% names(go)) {
    annotation_row <- annotation_row %>%
      mutate(class = go[[class_col]])
  }
  
  if (!is.null(type_col) && type_col %in% names(go)) {
    annotation_row <- annotation_row %>%
      mutate(type = go[[type_col]])
  }
  
  # Remove duplicates
  annotation_row <- annotation_row %>% distinct()
  rownames(annotation_row) <- annotation_row$sample
  annotation_row$sample <- NULL
  
  # Enhanced col_split logic - FIXED VERSION
  col_split <- NULL
  if (!is.null(col)) {
    # Use custom column data
    ann_tmp <- col %>%
      dplyr::select(sample = all_of(col_sample_col))
    rownames(ann_tmp) <- ann_tmp$sample
    
    # Check for split column in col data or go data
    if (!is.null(col_split_col)) {
      if (col_split_col %in% names(col)) {
        split_data <- col %>%
          dplyr::select(sample = all_of(col_sample_col), split_var = all_of(col_split_col)) %>%
          distinct()
        rownames(split_data) <- split_data$sample
        col_split <- split_data[rownames(ann_tmp), "split_var"]
      } else if (col_split_col %in% names(go)) {
        # Get split values from go data
        split_data <- go %>%
          dplyr::select(sample = all_of(sample_col), split_var = all_of(col_split_col)) %>%
          distinct()
        rownames(split_data) <- split_data$sample
        col_split <- split_data[rownames(ann_tmp), "split_var"]
      }
    }
    
    ann_tmp$sample <- NULL
    annotation_row <- ann_tmp
  } else {
    # Enhanced ordering logic using col_order_priority - DO THIS FIRST
    if (!is.null(col_order_priority) && length(col_order_priority) > 0) {
      # Use specified priority order, then add remaining columns
      priority_cols <- intersect(col_order_priority, names(annotation_row))
      remaining_cols <- setdiff(names(annotation_row), col_order_priority)
      order_cols <- c(priority_cols, remaining_cols)
    } else {
      # Default fallback order
      default_order <- c("region", "sex", "neurotransmitter", "class", "type")
      order_cols <- intersect(default_order, names(annotation_row))
      remaining_cols <- setdiff(names(annotation_row), default_order)
      order_cols <- c(order_cols, remaining_cols)
    }
    
    if (length(order_cols) > 0) {
      annotation_row <- annotation_row[do.call(order, annotation_row[order_cols]), , drop = FALSE]
    }
    
    # FIXED: Generate column split directly from sorted annotation_row
    if (!is.null(col_split_col) && col_split_col %in% names(annotation_row)) {
      # Use the sorted annotation_row data directly (maintains factor order)
      col_split <- annotation_row[[col_split_col]]
    } else if (!is.null(col_split_col) && col_split_col %in% names(go)) {
      # Fallback: regenerate from go data if column not in annotation_row
      split_data <- go %>%
        dplyr::select(sample = all_of(sample_col), split_var = all_of(col_split_col)) %>%
        distinct()
      rownames(split_data) <- split_data$sample
      col_split <- split_data[rownames(annotation_row), "split_var"]
    }
  }
  
  # ensure row order aligns
  mat <- mat[rownames(annotation_row), , drop = FALSE]
  
  # ---------- Base color palettes ----------
  base_region <- c(AMY="#ff7f0e", STR="#e377c2", PFC="#8c564b", iCTX="#1f77b4", MB="#9467bd", TH="#bcbd22",
                   HY="#d62728", HPF="#009E73")
  base_nt    <- c(Glut="#FFC000", GABA="#00B050", Dopa="#ff7f0e", Chol="#1f77b4",
                  Hist="#aa40fc", Sero="#e377c2", NN="#8c564b")
  base_class <- c(Glut="#FFC000", GABA="#00B050", Dopa="#ff7f0e", Chol="#1f77b4",
                  Hist="#aa40fc", Sero="#e377c2",
                  OPC="#8194AC", Oligo="#7030A0", Astrocyte="#17becf", Microglia="#ff9896",
                  Epen="#FF7F50", Vascular="#5F8BB2", Immune="#EF6347")
  base_sex   <- c(M="#55A0FB", F="#FFA0A0")
  base_type  <- c(N="#4CAF50", NN="#FF9800")  # Type colors (N=green, NN=orange)
  
  base_module <- c(
    G1 = '#59F5FD', "G1-1" = "#0080FF", "G1-2" = "#17becf", "G1-3" = '#59F5FD', "G1-4" = '#AFD8E6',
    "G2-1" = '#782170', "G2-2" = '#D86DCD', "G2-3" = '#E49FDD', "G2" = '#782170',
    G4 = '#FFC000', "G3-1" = '#00A050', "G3-2" = '#92D050', "G3-3" = '#83F28E',
    G5 = "red2", G6 = '#7E350E', G7 = "#0080FF", "G7-1" = "#5099FF", "G7-2" = "#0199EE",
    G8='#6868FF', "G8-1"='#A0B0F0', "G8-2"='#B5b0d5', G9="#95516E", G10 = '#2F4F4F',
    G11 = '#BE5014', G12 = "#767F50", G13= '#E49EDD', G14="#FF7F50", G15="#2F4F4F",
    G16="#8194AC", "G17-1"="#F5516E", "G17-2"="#E298A6", G19="#D9E7F5", G18="#5F8BB2",
    other= 'gray40'
  )
  
  # Enhanced column ordering logic: predefined order vs clustering
  if (!cluster_column_slices) {
    # Convert to factors with predefined levels (base color order)
    if ("region" %in% names(annotation_row)) {
      # Convert to factor with strict base_region order
      annotation_row$region <- factor(annotation_row$region, levels = names(base_region))
    }
    if ("neurotransmitter" %in% names(annotation_row)) {
      annotation_row$neurotransmitter <- factor(annotation_row$neurotransmitter, levels = names(base_nt))
    }
    if ("sex" %in% names(annotation_row)) {
      annotation_row$sex <- factor(annotation_row$sex, levels = names(base_sex))
    }
    if ("class" %in% names(annotation_row)) {
      annotation_row$class <- factor(annotation_row$class, levels = names(base_class))
    }
    if ("type" %in% names(annotation_row)) {
      annotation_row$type <- factor(annotation_row$type, levels = names(base_type))
    }
    
    # Reorder samples using col_order_priority with predefined factor levels
    if (!is.null(col_order_priority) && length(col_order_priority) > 0) {
      priority_cols <- intersect(col_order_priority, names(annotation_row))
      remaining_cols <- setdiff(names(annotation_row), col_order_priority)
      order_cols <- c(priority_cols, remaining_cols)
    } else {
      default_order <- c("region", "sex", "neurotransmitter", "class", "type")
      order_cols <- intersect(default_order, names(annotation_row))
      remaining_cols <- setdiff(names(annotation_row), default_order)
      order_cols <- c(order_cols, remaining_cols)
    }
    
    if (length(order_cols) > 0) {
      # Create ordering arguments with explicit factor handling
      order_args <- list()
      for (col in order_cols) {
        if (col %in% names(annotation_row)) {
          if (is.factor(annotation_row[[col]])) {
            # For factors, use as.numeric to respect level order
            order_args[[col]] <- as.numeric(annotation_row[[col]])
          } else {
            order_args[[col]] <- annotation_row[[col]]
          }
        }
      }
      
      # Apply ordering using the prepared arguments
      if (length(order_args) > 0) {
        order_idx <- do.call(order, order_args)
        annotation_row <- annotation_row[order_idx, , drop = FALSE]
        mat <- mat[rownames(annotation_row), , drop = FALSE]
        
        # FIXED: Generate col_split directly from sorted annotation_row
        if (!is.null(col_split_col) && col_split_col %in% names(annotation_row)) {
          # Use the sorted annotation_row data directly (maintains factor order)
          col_split <- annotation_row[[col_split_col]]
        } else if (!is.null(col_split) && !is.null(col_split_col) && col_split_col %in% names(go)) {
          # Fallback: regenerate from go data if column not in annotation_row
          split_data <- go %>%
            dplyr::select(sample = all_of(sample_col), split_var = all_of(col_split_col)) %>%
            distinct()
          rownames(split_data) <- split_data$sample
          col_split <- split_data[rownames(annotation_row), "split_var"]
        }
      }
    }
    
  }
  
  # ---------- 2) Enhanced row-side GO/module handling with missing gene inclusion ----------
  row_split <- NULL
  row_ha <- NULL
  row_annotations <- list()
  
  if (!is.null(row)) {
    if (!row_go_col %in% names(row)) {
      stop(sprintf("`row` must contain column `%s`.", row_go_col))
    }
    
    go_order <- pull(row, all_of(row_go_col))
    go_order <- as.character(go_order)
    go_order <- go_order[!duplicated(go_order)]
    
    # Find genes present and missing in matrix
    keep_go <- go_order[go_order %in% colnames(mat)]
    missing_go <- setdiff(go_order, keep_go)
    
    if (length(missing_go) > 0L) {
      message(sprintf("Note: %d term(s) in `%s` not found in matrix; adding with 0 values: %s",
                      length(missing_go), row_go_col, paste(head(missing_go, 5), collapse = ", ")))
      
      # Create matrix for missing genes with 0 values
      missing_mat <- matrix(0, nrow = nrow(mat), ncol = length(missing_go))
      colnames(missing_mat) <- missing_go
      rownames(missing_mat) <- rownames(mat)
      
      # Combine existing matrix with missing genes
      mat <- cbind(mat[, keep_go, drop = FALSE], missing_mat)
    } else {
      mat <- mat[, keep_go, drop = FALSE]
    }
    
    # Reorder matrix columns to match go_order
    mat <- mat[, go_order, drop = FALSE]
    
    # Handle multiple row grouping columns
    if (!is.null(row_module_col) && row_module_col %in% names(row)) {
      mod_vec <- pull(row, all_of(row_module_col))
      mod_vec <- mod_vec[match(go_order, pull(row, all_of(row_go_col)))]
      row_split <- mod_vec
      row_annotations$module <- mod_vec
    }
    
    # Handle additional category column
    if (!is.null(row_category_col) && row_category_col %in% names(row)) {
      cat_vec <- pull(row, all_of(row_category_col))
      cat_vec <- cat_vec[match(go_order, pull(row, all_of(row_go_col)))]
      row_annotations$category <- cat_vec
    }
  }
  
  # ---------- 3) Dynamic palettes ----------
  make_palette <- function(values, base_map) {
    vals <- as.character(unique(values))
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return(NULL)
    
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
  
  # Extract present values (region is now optional)
  present_regions <- if ("region" %in% names(annotation_row)) annotation_row$region else NULL
  present_nt      <- if ("neurotransmitter" %in% names(annotation_row)) annotation_row$neurotransmitter else NULL
  present_sex     <- if ("sex" %in% names(annotation_row)) annotation_row$sex else NULL
  present_class   <- if ("class" %in% names(annotation_row)) annotation_row$class else NULL
  present_type    <- if ("type" %in% names(annotation_row)) annotation_row$type else NULL
  present_mods    <- if (!is.null(row_split)) row_split else character(0)
  present_cats    <- if ("category" %in% names(row_annotations)) row_annotations$category else character(0)
  
  dyn_colors <- list()
  # Region is now optional
  if (!is.null(present_regions) && length(na.omit(present_regions)) > 0) {
    dyn_colors$region <- make_palette(present_regions, base_region)
  }
  if (!is.null(present_nt) && length(na.omit(present_nt)) > 0) {
    dyn_colors$neurotransmitter <- make_palette(present_nt, base_nt)
  }
  if (!is.null(present_sex) && length(na.omit(present_sex)) > 0) {
    dyn_colors$sex <- make_palette(present_sex, base_sex)
  }
  if (!is.null(present_class) && length(na.omit(present_class)) > 0) {
    dyn_colors$class <- make_palette(present_class, base_class)
  }
  # Add type colors
  if (!is.null(present_type) && length(na.omit(present_type)) > 0) {
    dyn_colors$type <- make_palette(present_type, base_type)
  }
  if (length(na.omit(present_mods)) > 0) {
    dyn_colors$module <- make_palette(present_mods, base_module)
  }
  if (length(na.omit(present_cats)) > 0) {
    dyn_colors$category <- make_palette(present_cats, base_map = c())
  }
  
  # Define consistent legend parameters for all legends
  legend_title_gp <- gpar(fontsize = 16, fontface = "bold")
  legend_labels_gp <- gpar(fontsize = 14)
  legend_grid_height <- unit(5, "mm")
  legend_grid_width <- unit(5, "mm")
  
  # ---------- 4) Build color scale with consistent styling ----------
  mat_t <- t(mat)
  
  # Create row_ha with multiple annotations if needed
  if (length(row_annotations) > 0) {
    # Build row annotation colors list
    row_colors_list <- list()
    row_legend_param_list <- list()
    
    for (ann_name in names(row_annotations)) {
      if (ann_name %in% names(dyn_colors)) {
        row_colors_list[[ann_name]] <- dyn_colors[[ann_name]]
        row_legend_param_list[[ann_name]] <- list(
          title_gp = legend_title_gp,
          labels_gp = legend_labels_gp,
          grid_height = legend_grid_height,
          grid_width = legend_grid_width
        )
      }
    }
    
    # Create rowAnnotation with all available annotations
    row_ha_args <- c(
      row_annotations,
      list(
        col = row_colors_list,
        annotation_name_side = "top",
        annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
        show_legend = TRUE,
        gp = gpar(col = NA),
        simple_anno_size = unit(1, "cm"),
        annotation_legend_param = row_legend_param_list
      )
    )
    
    row_ha <- do.call(rowAnnotation, row_ha_args)
  }
  
  # Keep original logFC cutoff logic - no changes to user's settings
  
  if (is.null(cutoff)) {
    rng <- range(mat_t, finite = TRUE)
    m <- max(abs(rng))
    at <- seq(-m, m, length.out = 5)
    col_fun <- circlize::colorRamp2(c(min(at), 0, max(at)),
                                    c("blue", "#FFFFFF", "red"))
    legend_param <- list(
      at = at,
      labels = c(-cutoff, -cutoff/2, 0, cutoff/2, cutoff),
      title = value_name,
      title_gp = legend_title_gp,
      labels_gp = legend_labels_gp
    )
  } else {
    cutoff <- as.numeric(cutoff)
    at <- c(-cutoff, -cutoff/2, 0, cutoff/2, cutoff)
    col_fun <- circlize::colorRamp2(c(-cutoff, 0, cutoff),
                                    c("blue", "#FFFFFF", "red"))
    legend_param <- list(
      at = at,
      labels = at,
      title = value_name,
      title_gp = legend_title_gp,
      labels_gp = legend_labels_gp
    )
  }
  
  # Build HeatmapAnnotation with consistent styling and custom order
  annotation_list <- list()
  
  # Use col_order_priority to determine the order of annotations
  if (!is.null(col_order_priority) && length(col_order_priority) > 0) {
    priority_cols <- intersect(col_order_priority, names(annotation_row))
    remaining_cols <- setdiff(names(annotation_row), col_order_priority)
    annotation_order <- c(priority_cols, remaining_cols)
  } else {
    # Default order
    default_order <- c("region", "sex", "neurotransmitter", "class", "type")
    priority_cols <- intersect(default_order, names(annotation_row))
    remaining_cols <- setdiff(names(annotation_row), default_order)
    annotation_order <- c(priority_cols, remaining_cols)
  }
  
  # Build annotation_list in the specified order
  for (ann_name in annotation_order) {
    if (ann_name %in% names(annotation_row)) {
      annotation_list[[ann_name]] <- annotation_row[[ann_name]]
    }
  }
  
  # Build annotation_legend_param dynamically based on actual annotations
  annotation_legend_param <- list()
  for (ann_name in names(annotation_list)) {
    annotation_legend_param[[ann_name]] <- list(
      title_gp = legend_title_gp,
      labels_gp = legend_labels_gp,
      grid_height = legend_grid_height,
      grid_width = legend_grid_width
    )
  }
  
  # Create annotation with consistent legend styling and custom order
  col_ha <- do.call(HeatmapAnnotation, c(
    annotation_list,
    list(
      col = dyn_colors,
      annotation_name_side = "right",
      annotation_legend_param = annotation_legend_param,
      annotation_name_gp = gpar(fontsize = 16, fontface = "bold"),
      simple_anno_size = unit(6, "mm")
    )
  ))
  
  col_ha@anno_list[[1]]@fun@height <- unit(10, "mm")
  col_ha@anno_list[[1]]@height <- unit(10, "mm")
  col_ha@anno_size[1] <- unit(10, "mm")
  
  
  # ---------- 5) Plot with format detection ----------
  # Detect file format from title extension
  file_ext <- tolower(tools::file_ext(title))
  
  if (file_ext == "png") {
    # PNG output
    png(title,
        height = (nrow(mat_t)/15 + 8) * 300,  # Convert to pixels (assuming 300 DPI)
        width  = (ncol(mat_t)/15 + 8) * 300,
        res = 300)  # High resolution for PNG
  } else {
    # Default to PDF output
    pdf(title,
        height = nrow(mat_t)/150 + 8,
        width  = ncol(mat_t)/15 + 8)
  }
  
  # Create cell function for custom coloring with logfc_cutoff
  cell_fun <- NULL
  logfc_legend <- NULL
  
  if (!is.null(logfc_cutoff)) {
    logfc_cutoff <- as.numeric(logfc_cutoff)
    
    cell_fun <- function(j, i, x, y, width, height, fill) {
      value <- mat_t[i, j]
      if (abs(value) <= logfc_cutoff) {
        # Grey color for values within logfc_cutoff range
        grid.rect(x = x, y = y, width = width, height = height,
                  gp = gpar(col = NA, fill = "#D3D3D3"))
      } else {
        # Use default coloring for significant values
        grid.rect(x = x, y = y, width = width, height = height,
                  gp = gpar(col = NA, fill = fill))
      }
    }
    
    # Create logfc_legend with CONSISTENT styling matching all other legends
    logfc_legend <- Legend(
      labels = sprintf("Not detected or not a DEG"),
      legend_gp = gpar(fill = "#D3D3D3"),
      title = NULL,  # No title to keep it clean
      title_gp = legend_title_gp,
      labels_gp = legend_labels_gp,
      grid_height = legend_grid_height,
      grid_width = legend_grid_width
    )
  }
  
  # Enhanced clustering logic: within-group clustering when col_split is provided
  final_cluster_columns <- cluster_columns
  
  # FIXED: Create heatmap with proper column split settings and within-group clustering
  ht <- Heatmap(
    mat_t,
    name = value_col,
    col  = col_fun,
    cell_fun = cell_fun,  # Apply custom cell coloring if logfc_cutoff is set
    heatmap_legend_param = legend_param,
    show_row_names    = TRUE,
    show_column_names = TRUE,
    cluster_rows      = cluster_rows,
    cluster_row_slices = cluster_row_slices,
    cluster_columns   = final_cluster_columns,  # Use enhanced clustering logic
    cluster_column_slices = cluster_column_slices,
    clustering_distance_rows    = distance,    # Distance method for row clustering
    clustering_method_rows      = method,      # Clustering method for rows
    clustering_distance_columns = distance,    # Distance method for column clustering
    clustering_method_columns   = method,      # Clustering method for columns
    column_split      = col_split,  # Enhanced split logic with within-group clustering
    row_split         = row_split,
    column_gap        = column_gap,  # Control column split gap width
    row_gap           = row_gap,     # Control row split gap width
    gap               = unit(1, "cm"),
    left_annotation   = row_ha,
    top_annotation    = col_ha,
    row_names_side    = "right",
    column_names_side = "bottom",
    row_names_gp      = gpar(fontsize = row_names_gp),
    column_names_gp   = gpar(fontsize = column_names_gp),
    use_raster        = FALSE,
    row_title         = NULL,
    column_title      = NULL,
    show_column_dend  = if (cluster_columns) TRUE else FALSE,
    
  )
  
  # Draw heatmap with optional logfc_cutoff legend
  if (!is.null(logfc_legend)) {
    ht_drawn <- draw(ht, annotation_legend_list = list(logfc_legend), padding = unit(c(35,5, 5, 5), "mm"))
  } else {
    ht_drawn <- draw(ht, padding = unit(c(35, 5, 5, 5), "mm"))
  }
  
  ro <- unlist(row_order(ht_drawn))
  co <- unlist(column_order(ht_drawn))
  dev.off()
  
  row_annhn=unique(row_data[,c(1,3)])
  rownames(row_annhn)=row_annhn$GO
  a=row_annhn[rownames(mat_t[ro,co]),2:1]
  colnames(a)=c("Group","GO/Cellgroup")
  aa=cbind(a,mat_t[ro,co])
  cc=cbind.data.frame("",cbind.data.frame(rownames(t(annotation_row[colnames(mat_t[ro,co]),])),t(annotation_row[colnames(mat_t[ro,co]),])))
  colnames(cc)[1:2]=colnames(aa)[1:2]
  f=rbind.data.frame(cc,colnames(aa),aa)
  colnames(f) <- NULL 
  # save matrix
  #if (grepl("nocluster", title, ignore.case = TRUE)) {
  #final_mat <- mat_t[ro, co, drop = FALSE]
  out_file <- paste0(file_path_sans_ext(title), "_plotdata.csv")
  write.csv(f, out_file, quote = TRUE,row.names = F)
  
  return(mat_t[ro, co, drop = FALSE])
}


heatmap_final2 <- function(
    go,
    row = NULL,
    col = NULL,
    sample_col = "sample",
    term_col = "term_name",
    value_col = "value",
    value_name = "logFC",
    region_col = NULL,
    neurotransmitter_col = NULL,
    sex_col = NULL,
    class_col = NULL,
    type_col = NULL,
    #
    row_go_col = "Gene_name",
    row_module_col = "Category",
    row_category_col = NULL,
    col_sample_col = "sample",
    col_order_priority = c("region","sex","neurotransmitter","class","type"),
    col_split_col = 'split',
    #
    cluster_rows = TRUE,
    cluster_row_slices = FALSE,
    cluster_columns = TRUE,
    cluster_column_slices = FALSE,
    distance = "euclidean",      # Distance method for clustering, default is "euclidean"
    method = "complete",         # Clustering method, default is "complete","ward.D2"
    cutoff = NULL,
    logfc_cutoff = NULL,
    column_gap = unit(2, "mm"),  # Gap between column splits
    row_gap = unit(2, "mm"),     # Gap between row splits
    row_names_gp = 24,
    column_names_gp = 9,
    title = "_heatmap.pdf",
    height = 40, width = 10
){
  
  # ---------- 0) Prepare matrix ----------
  heatmap_data <- go %>%
    dplyr::select(
      .s = all_of(sample_col),
      .t = all_of(term_col),
      .v = all_of(value_col)
    ) %>%
    mutate(.v = suppressWarnings(as.numeric(as.character(.v)))) %>%
    pivot_wider(names_from = .t, values_from = .v)
  
  mat <- as.data.frame(heatmap_data)
  rownames(mat) <- mat$.s
  mat$.s <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  
  if (!is.null(cutoff)) {
    mat <- pmax(pmin(mat, cutoff), -cutoff)
  }
  
  set.seed(111)
  
  # ---------- 1) Column (sample) annotations ----------
  # Start with required sample column
  annotation_row <- go %>%
    dplyr::select(sample = all_of(sample_col))
  
  # Add region column if provided
  if (!is.null(region_col) && region_col %in% names(go)) {
    annotation_row <- annotation_row %>%
      mutate(region = go[[region_col]])
  }
  
  # Add optional columns only if they are provided
  if (!is.null(neurotransmitter_col) && neurotransmitter_col %in% names(go)) {
    annotation_row <- annotation_row %>%
      mutate(neurotransmitter = go[[neurotransmitter_col]])
  }
  
  if (!is.null(sex_col) && sex_col %in% names(go)) {
    annotation_row <- annotation_row %>%
      mutate(sex = go[[sex_col]])
  }
  
  if (!is.null(class_col) && class_col %in% names(go)) {
    annotation_row <- annotation_row %>%
      mutate(class = go[[class_col]])
  }
  
  if (!is.null(type_col) && type_col %in% names(go)) {
    annotation_row <- annotation_row %>%
      mutate(type = go[[type_col]])
  }
  
  # Remove duplicates
  annotation_row <- annotation_row %>% distinct()
  rownames(annotation_row) <- annotation_row$sample
  annotation_row$sample <- NULL
  
  # Enhanced col_split logic - FIXED VERSION
  col_split <- NULL
  if (!is.null(col)) {
    # Use custom column data
    ann_tmp <- col %>%
      dplyr::select(sample = all_of(col_sample_col))
    rownames(ann_tmp) <- ann_tmp$sample
    
    # Check for split column in col data or go data
    if (!is.null(col_split_col)) {
      if (col_split_col %in% names(col)) {
        split_data <- col %>%
          dplyr::select(sample = all_of(col_sample_col), split_var = all_of(col_split_col)) %>%
          distinct()
        rownames(split_data) <- split_data$sample
        col_split <- split_data[rownames(ann_tmp), "split_var"]
      } else if (col_split_col %in% names(go)) {
        # Get split values from go data
        split_data <- go %>%
          dplyr::select(sample = all_of(sample_col), split_var = all_of(col_split_col)) %>%
          distinct()
        rownames(split_data) <- split_data$sample
        col_split <- split_data[rownames(ann_tmp), "split_var"]
      }
    }
    
    ann_tmp$sample <- NULL
    annotation_row <- ann_tmp
  } else {
    # Enhanced ordering logic using col_order_priority - DO THIS FIRST
    if (!is.null(col_order_priority) && length(col_order_priority) > 0) {
      # Use specified priority order, then add remaining columns
      priority_cols <- intersect(col_order_priority, names(annotation_row))
      remaining_cols <- setdiff(names(annotation_row), col_order_priority)
      order_cols <- c(priority_cols, remaining_cols)
    } else {
      # Default fallback order
      default_order <- c("region", "sex", "neurotransmitter", "class", "type")
      order_cols <- intersect(default_order, names(annotation_row))
      remaining_cols <- setdiff(names(annotation_row), default_order)
      order_cols <- c(order_cols, remaining_cols)
    }
    
    if (length(order_cols) > 0) {
      annotation_row <- annotation_row[do.call(order, annotation_row[order_cols]), , drop = FALSE]
    }
    
    # FIXED: Generate column split directly from sorted annotation_row
    if (!is.null(col_split_col) && col_split_col %in% names(annotation_row)) {
      # Use the sorted annotation_row data directly (maintains factor order)
      col_split <- annotation_row[[col_split_col]]
    } else if (!is.null(col_split_col) && col_split_col %in% names(go)) {
      # Fallback: regenerate from go data if column not in annotation_row
      split_data <- go %>%
        dplyr::select(sample = all_of(sample_col), split_var = all_of(col_split_col)) %>%
        distinct()
      rownames(split_data) <- split_data$sample
      col_split <- split_data[rownames(annotation_row), "split_var"]
    }
  }
  
  # ensure row order aligns
  mat <- mat[rownames(annotation_row), , drop = FALSE]
  
  # ---------- Base color palettes ----------
  base_region <- c(AMY="#ff7f0e", STR="#e377c2", PFC="#8c564b", iCTX="#1f77b4", MB="#9467bd", TH="#bcbd22",
                   HY="#d62728", HPF="#009E73")
  base_nt    <- c(Glut="#FFC000", GABA="#00B050", Dopa="#ff7f0e", Chol="#1f77b4",
                  Hist="#aa40fc", Sero="#e377c2", NN="#8c564b")
  base_class <- c(Glut="#FFC000", GABA="#00B050", Dopa="#ff7f0e", Chol="#1f77b4",
                  Hist="#aa40fc", Sero="#e377c2",
                  OPC="#8194AC", Oligo="#7030A0", Astrocyte="#17becf", Microglia="#ff9896",
                  Epen="#FF7F50", Vascular="#5F8BB2", Immune="#EF6347")
  base_sex   <- c(M="#55A0FB", F="#FFA0A0")
  base_type  <- c(N="#4CAF50", NN="#FF9800")  # Type colors (N=green, NN=orange)
  
  base_module <- c(
    G1 = '#59F5FD', "G1-1" = "#0080FF", "G1-2" = "#17becf", "G1-3" = '#59F5FD', "G1-4" = '#AFD8E6',
    "G2-1" = '#782170', "G2-2" = '#D86DCD', "G2-3" = '#E49FDD', "G2" = '#782170',
    G4 = '#FFC000', "G3-1" = '#00A050', "G3-2" = '#92D050', "G3-3" = '#83F28E',
    G5 = "red2", G6 = '#7E350E', G7 = "#0080FF", "G7-1" = "#5099FF", "G7-2" = "#0199EE",
    G8='#6868FF', "G8-1"='#A0B0F0', "G8-2"='#B5b0d5', G9="#95516E", G10 = '#2F4F4F',
    G11 = '#BE5014', G12 = "#767F50", G13= '#E49EDD', G14="#FF7F50", G15="#2F4F4F",
    G16="#8194AC", "G17-1"="#F5516E", "G17-2"="#E298A6", G19="#D9E7F5", G18="#5F8BB2",
    other= 'gray40'
  )
  
  # Enhanced column ordering logic: predefined order vs clustering
  if (!cluster_column_slices) {
    # Convert to factors with predefined levels (base color order)
    if ("region" %in% names(annotation_row)) {
      # Convert to factor with strict base_region order
      annotation_row$region <- factor(annotation_row$region, levels = names(base_region))
    }
    if ("neurotransmitter" %in% names(annotation_row)) {
      annotation_row$neurotransmitter <- factor(annotation_row$neurotransmitter, levels = names(base_nt))
    }
    if ("sex" %in% names(annotation_row)) {
      annotation_row$sex <- factor(annotation_row$sex, levels = names(base_sex))
    }
    if ("class" %in% names(annotation_row)) {
      annotation_row$class <- factor(annotation_row$class, levels = names(base_class))
    }
    if ("type" %in% names(annotation_row)) {
      annotation_row$type <- factor(annotation_row$type, levels = names(base_type))
    }
    
    # Reorder samples using col_order_priority with predefined factor levels
    if (!is.null(col_order_priority) && length(col_order_priority) > 0) {
      priority_cols <- intersect(col_order_priority, names(annotation_row))
      remaining_cols <- setdiff(names(annotation_row), col_order_priority)
      order_cols <- c(priority_cols, remaining_cols)
    } else {
      default_order <- c("region", "sex", "neurotransmitter", "class", "type")
      order_cols <- intersect(default_order, names(annotation_row))
      remaining_cols <- setdiff(names(annotation_row), default_order)
      order_cols <- c(order_cols, remaining_cols)
    }
    
    if (length(order_cols) > 0) {
      # Create ordering arguments with explicit factor handling
      order_args <- list()
      for (col in order_cols) {
        if (col %in% names(annotation_row)) {
          if (is.factor(annotation_row[[col]])) {
            # For factors, use as.numeric to respect level order
            order_args[[col]] <- as.numeric(annotation_row[[col]])
          } else {
            order_args[[col]] <- annotation_row[[col]]
          }
        }
      }
      
      # Apply ordering using the prepared arguments
      if (length(order_args) > 0) {
        order_idx <- do.call(order, order_args)
        annotation_row <- annotation_row[order_idx, , drop = FALSE]
        mat <- mat[rownames(annotation_row), , drop = FALSE]
        
        # FIXED: Generate col_split directly from sorted annotation_row
        if (!is.null(col_split_col) && col_split_col %in% names(annotation_row)) {
          # Use the sorted annotation_row data directly (maintains factor order)
          col_split <- annotation_row[[col_split_col]]
        } else if (!is.null(col_split) && !is.null(col_split_col) && col_split_col %in% names(go)) {
          # Fallback: regenerate from go data if column not in annotation_row
          split_data <- go %>%
            dplyr::select(sample = all_of(sample_col), split_var = all_of(col_split_col)) %>%
            distinct()
          rownames(split_data) <- split_data$sample
          col_split <- split_data[rownames(annotation_row), "split_var"]
        }
      }
    }
    
  }
  
  # ---------- 2) Enhanced row-side GO/module handling with missing gene inclusion ----------
  row_split <- NULL
  row_ha <- NULL
  row_annotations <- list()
  
  if (!is.null(row)) {
    if (!row_go_col %in% names(row)) {
      stop(sprintf("`row` must contain column `%s`.", row_go_col))
    }
    
    go_order <- pull(row, all_of(row_go_col))
    go_order <- as.character(go_order)
    go_order <- go_order[!duplicated(go_order)]
    
    # Find genes present and missing in matrix
    keep_go <- go_order[go_order %in% colnames(mat)]
    missing_go <- setdiff(go_order, keep_go)
    
    if (length(missing_go) > 0L) {
      message(sprintf("Note: %d term(s) in `%s` not found in matrix; adding with 0 values: %s",
                      length(missing_go), row_go_col, paste(head(missing_go, 5), collapse = ", ")))
      
      # Create matrix for missing genes with 0 values
      missing_mat <- matrix(0, nrow = nrow(mat), ncol = length(missing_go))
      colnames(missing_mat) <- missing_go
      rownames(missing_mat) <- rownames(mat)
      
      # Combine existing matrix with missing genes
      mat <- cbind(mat[, keep_go, drop = FALSE], missing_mat)
    } else {
      mat <- mat[, keep_go, drop = FALSE]
    }
    
    # Reorder matrix columns to match go_order
    mat <- mat[, go_order, drop = FALSE]
    
    # Handle multiple row grouping columns
    if (!is.null(row_module_col) && row_module_col %in% names(row)) {
      mod_vec <- pull(row, all_of(row_module_col))
      mod_vec <- mod_vec[match(go_order, pull(row, all_of(row_go_col)))]
      row_split <- mod_vec
      row_annotations$module <- mod_vec
    }
    
    # Handle additional category column
    if (!is.null(row_category_col) && row_category_col %in% names(row)) {
      cat_vec <- pull(row, all_of(row_category_col))
      cat_vec <- cat_vec[match(go_order, pull(row, all_of(row_go_col)))]
      row_annotations$category <- cat_vec
    }
  }
  
  # ---------- 3) Dynamic palettes ----------
  make_palette <- function(values, base_map) {
    vals <- as.character(unique(values))
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return(NULL)
    
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
  
  # Extract present values (region is now optional)
  present_regions <- if ("region" %in% names(annotation_row)) annotation_row$region else NULL
  present_nt      <- if ("neurotransmitter" %in% names(annotation_row)) annotation_row$neurotransmitter else NULL
  present_sex     <- if ("sex" %in% names(annotation_row)) annotation_row$sex else NULL
  present_class   <- if ("class" %in% names(annotation_row)) annotation_row$class else NULL
  present_type    <- if ("type" %in% names(annotation_row)) annotation_row$type else NULL
  present_mods    <- if (!is.null(row_split)) row_split else character(0)
  present_cats    <- if ("category" %in% names(row_annotations)) row_annotations$category else character(0)
  
  dyn_colors <- list()
  # Region is now optional
  if (!is.null(present_regions) && length(na.omit(present_regions)) > 0) {
    dyn_colors$region <- make_palette(present_regions, base_region)
  }
  if (!is.null(present_nt) && length(na.omit(present_nt)) > 0) {
    dyn_colors$neurotransmitter <- make_palette(present_nt, base_nt)
  }
  if (!is.null(present_sex) && length(na.omit(present_sex)) > 0) {
    dyn_colors$sex <- make_palette(present_sex, base_sex)
  }
  if (!is.null(present_class) && length(na.omit(present_class)) > 0) {
    dyn_colors$class <- make_palette(present_class, base_class)
  }
  # Add type colors
  if (!is.null(present_type) && length(na.omit(present_type)) > 0) {
    dyn_colors$type <- make_palette(present_type, base_type)
  }
  if (length(na.omit(present_mods)) > 0) {
    dyn_colors$module <- make_palette(present_mods, base_module)
  }
  if (length(na.omit(present_cats)) > 0) {
    dyn_colors$category <- make_palette(present_cats, base_map = c())
  }
  
  # Define consistent legend parameters for all legends
  legend_title_gp <- gpar(fontsize = 16, fontface = "bold")
  legend_labels_gp <- gpar(fontsize = 14)
  legend_grid_height <- unit(5, "mm")
  legend_grid_width <- unit(5, "mm")
  
  # ---------- 4) Build color scale with consistent styling ----------
  mat_t <- t(mat)
  
  # Create row_ha with multiple annotations if needed
  if (length(row_annotations) > 0) {
    # Build row annotation colors list
    row_colors_list <- list()
    row_legend_param_list <- list()
    
    for (ann_name in names(row_annotations)) {
      if (ann_name %in% names(dyn_colors)) {
        row_colors_list[[ann_name]] <- dyn_colors[[ann_name]]
        row_legend_param_list[[ann_name]] <- list(
          title_gp = legend_title_gp,
          labels_gp = legend_labels_gp,
          grid_height = legend_grid_height,
          grid_width = legend_grid_width
        )
      }
    }
    
    # Create rowAnnotation with all available annotations
    row_ha_args <- c(
      row_annotations,
      list(
        col = row_colors_list,
        annotation_name_side = "top",
        annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
        show_legend = TRUE,
        gp = gpar(col = NA),
        simple_anno_size = unit(1, "cm"),
        annotation_legend_param = row_legend_param_list
      )
    )
    
    row_ha <- do.call(rowAnnotation, row_ha_args)
  }
  
  # Keep original logFC cutoff logic - no changes to user's settings
  
  if (is.null(cutoff)) {
    rng <- range(mat_t, finite = TRUE)
    m <- max(abs(rng))
    at <- seq(-m, m, length.out = 5)
    col_fun <- circlize::colorRamp2(c(min(at), 0, max(at)),
                                    c("blue", "#FFFFFF", "red"))
    legend_param <- list(
      at = at,
      labels = c(-cutoff, -cutoff/2, 0, cutoff/2, cutoff),
      title = value_name,
      title_gp = legend_title_gp,
      labels_gp = legend_labels_gp
    )
  } else {
    cutoff <- as.numeric(cutoff)
    at <- c(-cutoff, -cutoff/2, 0, cutoff/2, cutoff)
    col_fun <- circlize::colorRamp2(c(-cutoff, 0, cutoff),
                                    c("blue", "#FFFFFF", "red"))
    legend_param <- list(
      at = at,
      labels = at,
      title = value_name,
      title_gp = legend_title_gp,
      labels_gp = legend_labels_gp
    )
  }
  
  # Build HeatmapAnnotation with consistent styling and custom order
  annotation_list <- list()
  
  # Use col_order_priority to determine the order of annotations
  if (!is.null(col_order_priority) && length(col_order_priority) > 0) {
    priority_cols <- intersect(col_order_priority, names(annotation_row))
    remaining_cols <- setdiff(names(annotation_row), col_order_priority)
    annotation_order <- c(priority_cols, remaining_cols)
  } else {
    # Default order
    default_order <- c("region", "sex", "neurotransmitter", "class", "type")
    priority_cols <- intersect(default_order, names(annotation_row))
    remaining_cols <- setdiff(names(annotation_row), default_order)
    annotation_order <- c(priority_cols, remaining_cols)
  }
  
  # Build annotation_list in the specified order
  for (ann_name in annotation_order) {
    if (ann_name %in% names(annotation_row)) {
      annotation_list[[ann_name]] <- annotation_row[[ann_name]]
    }
  }
  
  # Build annotation_legend_param dynamically based on actual annotations
  annotation_legend_param <- list()
  for (ann_name in names(annotation_list)) {
    annotation_legend_param[[ann_name]] <- list(
      title_gp = legend_title_gp,
      labels_gp = legend_labels_gp,
      grid_height = legend_grid_height,
      grid_width = legend_grid_width
    )
  }
  
  # Create annotation with consistent legend styling and custom order
  col_ha <- do.call(HeatmapAnnotation, c(
    annotation_list,
    list(
      col = dyn_colors,
      annotation_name_side = "right",
      annotation_legend_param = annotation_legend_param,
      annotation_name_gp = gpar(fontsize = 16, fontface = "bold"),
      simple_anno_size = unit(6, "mm")
    )
  ))
  
  col_ha@anno_list[[1]]@fun@height <- unit(10, "mm")
  col_ha@anno_list[[1]]@height <- unit(10, "mm")
  col_ha@anno_size[1] <- unit(10, "mm")
  
  
  # ---------- 5) Plot with format detection ----------
  # Detect file format from title extension
  file_ext <- tolower(tools::file_ext(title))
  
  if (file_ext == "png") {
    # PNG output
    png(title,
        height = (nrow(mat_t)/15 + 8) * 300,  # Convert to pixels (assuming 300 DPI)
        width  = (ncol(mat_t)/15 + 8) * 300,
        res = 300)  # High resolution for PNG
  } else {
    # Default to PDF output
    pdf(title,
        height = nrow(mat_t)/15 + 8,
        width  = ncol(mat_t)/15 + 8)
  }
  
  # Create cell function for custom coloring with logfc_cutoff
  cell_fun <- NULL
  logfc_legend <- NULL
  
  if (!is.null(logfc_cutoff)) {
    logfc_cutoff <- as.numeric(logfc_cutoff)
    
    cell_fun <- function(j, i, x, y, width, height, fill) {
      value <- mat_t[i, j]
      if (abs(value) <= logfc_cutoff) {
        # Grey color for values within logfc_cutoff range
        grid.rect(x = x, y = y, width = width, height = height,
                  gp = gpar(col = NA, fill = "#D3D3D3"))
      } else {
        # Use default coloring for significant values
        grid.rect(x = x, y = y, width = width, height = height,
                  gp = gpar(col = NA, fill = fill))
      }
    }
    
    # Create logfc_legend with CONSISTENT styling matching all other legends
    logfc_legend <- Legend(
      labels = sprintf("Not detected or not a DEG"),
      legend_gp = gpar(fill = "#D3D3D3"),
      title = NULL,  # No title to keep it clean
      title_gp = legend_title_gp,
      labels_gp = legend_labels_gp,
      grid_height = legend_grid_height,
      grid_width = legend_grid_width
    )
  }
  
  # Enhanced clustering logic: within-group clustering when col_split is provided
  final_cluster_columns <- cluster_columns
  
  # FIXED: Create heatmap with proper column split settings and within-group clustering
  ht <- Heatmap(
    mat_t,
    name = value_col,
    col  = col_fun,
    cell_fun = cell_fun,  # Apply custom cell coloring if logfc_cutoff is set
    heatmap_legend_param = legend_param,
    show_row_names    = TRUE,
    show_column_names = TRUE,
    cluster_rows      = cluster_rows,
    cluster_row_slices = cluster_row_slices,
    cluster_columns   = final_cluster_columns,  # Use enhanced clustering logic
    cluster_column_slices = cluster_column_slices,
    clustering_distance_rows    = distance,    # Distance method for row clustering
    clustering_method_rows      = method,      # Clustering method for rows
    clustering_distance_columns = distance,    # Distance method for column clustering
    clustering_method_columns   = method,      # Clustering method for columns
    column_split      = col_split,  # Enhanced split logic with within-group clustering
    row_split         = row_split,
    column_gap        = column_gap,  # Control column split gap width
    row_gap           = row_gap,     # Control row split gap width
    gap               = unit(1, "cm"),
    left_annotation   = row_ha,
    top_annotation    = col_ha,
    row_names_side    = "right",
    column_names_side = "bottom",
    row_names_gp      = gpar(fontsize = row_names_gp),
    column_names_gp   = gpar(fontsize = column_names_gp),
    use_raster        = FALSE,
    row_title         = NULL,
    column_title      = NULL,
    show_column_dend  = if (cluster_columns) TRUE else FALSE,
    
  )
  
  # Draw heatmap with optional logfc_cutoff legend
  if (!is.null(logfc_legend)) {
    ht_drawn <- draw(ht, annotation_legend_list = list(logfc_legend), padding = unit(c(35,5, 5, 5), "mm"))
  } else {
    ht_drawn <- draw(ht, padding = unit(c(35, 5, 5, 5), "mm"))
  }
  
  ro <- unlist(row_order(ht_drawn))
  co <- unlist(column_order(ht_drawn))
  dev.off()
  

  row_annhn=unique(row_data[,c(1,3)])
  rownames(row_annhn)=row_annhn$GO
  a=row_annhn[rownames(mat_t[ro,co]),2:1]
  colnames(a)=c("Group","GO/Cellgroup")
  aa=cbind(a,mat_t[ro,co])
  cc=cbind.data.frame("",cbind.data.frame(rownames(t(annotation_row[colnames(mat_t[ro,co]),])),t(annotation_row[colnames(mat_t[ro,co]),])))
  colnames(cc)[1:2]=colnames(aa)[1:2]
  f=rbind.data.frame(cc,colnames(aa),aa)
  colnames(f) <- NULL 
  # save matrix
  #if (grepl("nocluster", title, ignore.case = TRUE)) {
  #final_mat <- mat_t[ro, co, drop = FALSE]
  out_file <- paste0(file_path_sans_ext(title), "_plotdata.csv")
  write.csv(f, out_file, quote = TRUE,row.names = F)
  #}
  return(mat_t[ro, co, drop = FALSE])
}


filter_degs_and_category <- function(degs_astrocyte, astrocyte_category,
                                     logfc_cutoff = 0.1,
                                     agg_fun = function(x) max(x, na.rm = TRUE),
                                     n_cells = 1,
                                     n_genes = 1) {
  # 1) Keep only genes present in astrocyte_category
  degs_astrocyte <- degs_astrocyte %>%
    semi_join(astrocyte_category, by = c("term_name" = "Gene_name"))
  
  # 2) Aggregate duplicate term_name × sample
  degs_agg <- degs_astrocyte %>%
    group_by(term_name, sample) %>%
    summarise(value = agg_fun(value), .groups = "drop")
  
  # 3) Pivot to wide matrix
  mat <- degs_agg %>%
    pivot_wider(names_from = sample, values_from = value) %>%
    as.data.frame()
  rownames(mat) <- mat$term_name
  mat$term_name <- NULL
  mat <- as.matrix(mat)
  
  # 4) Apply cutoff
  mat[is.na(mat) | abs(mat) < logfc_cutoff] <- NA
  
  # 5a) Gene-level filter (require DEG in >= n_cells samples)
  keep_genes <- rowSums(!is.na(mat)) >= n_cells
  mat <- mat[keep_genes, , drop = FALSE]
  
  # 5b) Cell-level filter (require >= n_genes DEGs per sample)
  keep_cells <- colSums(!is.na(mat)) >= n_genes
  mat <- mat[, keep_cells, drop = FALSE]
  
  # 6) Extract surviving term_name & sample
  keep_terms <- rownames(mat)
  keep_samples <- colnames(mat)
  
  # 7) Subset original degs_astrocyte (keep all 7 cols)
  degs_astrocyte_filtered <- degs_astrocyte %>%
    filter(term_name %in% keep_terms, sample %in% keep_samples)
  
  # 8) Filter astrocyte_category
  astrocyte_category_filtered <- astrocyte_category %>%
    semi_join(degs_astrocyte_filtered, by = c("Gene_name" = "term_name"))
  
  return(list(
    degs_astrocyte = degs_astrocyte_filtered,      # 7-column output
    astrocyte_category = astrocyte_category_filtered
  ))
}

filter_col_and_row <- function(main_data, row_data,
                               logfc_cutoff = NULL,
                               agg_fun = function(x) max(x, na.rm = TRUE),
                               n_cells = 1,
                               n_genes = 1,
                               row_col = "Gene_name") {
  
  #   n_cells     : Minimum number of samples a gene must appear in to be kept
  #   n_genes     : Minimum number of genes required per sample to keep sample
  #   row_col     : Column name in row_data corresponding to term_name in main_data
  
  # 1) Keep only genes present in row_data
  main_data <- main_data %>%
    semi_join(row_data, by = c("term_name" = row_col))
  
  # 2) Aggregate duplicate term_name × sample
  main_agg <- main_data %>%
    group_by(term_name, sample) %>%
    summarise(value = agg_fun(value), .groups = "drop")
  
  # 3) Pivot to wide matrix
  mat <- main_agg %>%
    pivot_wider(names_from = sample, values_from = value) %>%
    as.data.frame()
  rownames(mat) <- mat$term_name
  mat$term_name <- NULL
  mat <- as.matrix(mat)
  
  # 4) Apply cutoff
  if (!is.null(logfc_cutoff)) {
    mat[is.na(mat) | abs(mat) < logfc_cutoff] <- NA
  }
  
  # 5) Iteratively filter rows and columns until stable
  repeat {
    nrow_before <- nrow(mat)
    ncol_before <- ncol(mat)
    
    # Gene-level filter: keep genes present in >= n_cells samples
    keep_genes <- rowSums(!is.na(mat)) >= n_cells
    mat <- mat[keep_genes, , drop = FALSE]
    
    # Cell-level filter: keep samples with >= n_genes DEGs
    keep_cells <- colSums(!is.na(mat)) >= n_genes
    mat <- mat[, keep_cells, drop = FALSE]
    
    # Stop if rows and columns are stable
    if (nrow(mat) == nrow_before && ncol(mat) == ncol_before) break
    
    # Stop if matrix becomes empty
    if (nrow(mat) == 0 || ncol(mat) == 0) break
  }
  
  # 6) Extract surviving term_name & sample
  keep_terms <- rownames(mat)
  keep_samples <- colnames(mat)
  
  # 7) Subset original main_data (keep all original columns)
  main_data_filtered <- main_data %>%
    filter(term_name %in% keep_terms, sample %in% keep_samples)
  
  # 8) Filter row_data using filtered term_name
  row_data_filtered <- row_data %>%
    semi_join(main_data_filtered, by = setNames("term_name", row_col))
  
  return(list(
    main_data = main_data_filtered,
    row_data  = row_data_filtered
  ))
}

#########for degGo############
go_paths=c(Sys.glob("/data2st1/junyi/output/atac1112/dar/gsea/MASTNG_dar_annotated_avg/MASTNG_dar_annotated_avg/go_gsea_heatmap/MASTNG_dar_annotated_avg_merged_up.csv"))
           

# ##GABA adjP#########################################
# for (go_path in go_paths) {
#   message("Processing: ", go_path)
  
#   go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1)
#   go_results_tmp=go_results_tmp[go_results_tmp$Region!="iCTX",]
#   go_results_tmp=go_results_tmp[go_results_tmp$source%in%c("GO:BP","GO:MF"),]
#   go_results_tmp=go_results_tmp[go_results_tmp$Neurotransmitter=="GABA",]
#   #go_results_tmp$padj=p.adjust(go_results_tmp$p_value,method = "BH")
#   # go_results_tmp$nlogp=-log10(go_results_tmp$padj)*sign(go_results_tmp$nlogp)
#   #  go_results_tmp$nlog10p=-log10(go_results_tmp$padj)
  
  
#   colnames(go_results_tmp)=gsub("Region.Subclass","Region Subclass",colnames(go_results_tmp))
#   go_results_tmp_unique <- go_results_tmp %>%
#     distinct(`Region Subclass`, Sex, GO, .keep_all = TRUE)
#   if (nrow(go_results_tmp_unique) < nrow(go_results_tmp)) {
#     go_path_unique <- sub("\\.csv$", "_unique.csv", go_path)
#     write.csv(go_results_tmp_unique, go_path_unique, row.names = FALSE)
#   }
  
#   # go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1) %>%
#   #   filter(Sex %in% sexes)
  
#   #go_results <- switch(target_data,
#   #                     "SUS" = go_results_tmp_unique
#   #)
#   go_results=go_results_tmp_unique
  
#   output_dir_path <- dirname(go_path)
#   file_name <- tools::file_path_sans_ext(basename(go_path))
#   #output_dir_path <- file.path(output_dir_path, file_name)
#   output_dir_path <- file.path(output_dir_path, paste0(file_name,"_small_onlyGOBPMF_GABA_adjP"))
  
#   if (!dir.exists(output_dir_path)) {
#     dir.create(output_dir_path, recursive = TRUE)
#   }
  
#   subsets <- list(
#     list(tag = "N", filter = function(df) df %>% filter(Neurotransmitter != "NN"))
#     #list(tag = "NN",    filter = function(df) df %>% filter(Neurotransmitter == "NN"))
#   )
  
#   group_col <- 'predicted_group' # 'Module.P','predicted_group' or 'final_group'
#   use_row_category <-TRUE
#   class_target <-  NULL # c("Microglia")
#   # go_results$Module.P <- gsub("^G", "M", go_results[[group_col]])
#   go_results <- go_results %>%
#     mutate(!!group_col := ifelse(.data[[group_col]] == "others", "other", .data[[group_col]]))
#   module_levels <- c("G1","G2-1","G2-2","G3-1","G3-2","G3-3","G4","G5","G6","G7","G7-1","G7-2","G8","G8-1","G8-2","G9","G10","G11","G12","other") # c("innate", "overlapped", "adaptive", "inflammatory")
  
#   if (use_row_category) {
#     go_group <- go_results %>%
#       dplyr::select(GO,Neurotransmitter,all_of(group_col)) %>%
#       rename(module = all_of(group_col)) %>%
#       {
#         if (!is.null(module_levels)) {
#           mutate(., module = factor(module, levels = module_levels))
#         } else {
#           mutate(., module = factor(module))
#         }
#       } %>%
#       arrange(module)
#   } else {
#     go_group <- go_results %>%
#       dplyr::select(GO,Neurotransmitter)
#   }
  
  
#   # Generate heatmaps for each subset and clustering configuration
#   # split N and NN for separate heatmaps
#   # 生成不同子集与聚类配置的热图（包含：不过滤GO/列的unfiltered版本）
#   # 提示：确保 dplyr 在最后加载，或统一使用 dplyr:: 前缀（本代码已使用前缀）
#   # library(dplyr)  # 如果你不想每次都写 dplyr::，可把这句放到所有包加载的最后
  
#   # split N and NN for separate heatmaps
#   for (subset_cfg in subsets) {
#     subset_tag <- subset_cfg$tag
#     message("Generating subset: ", subset_tag)
    
#     # ---- 构建主数据（列数据） ----
#     go_results_subset <- subset_cfg$filter(go_results) %>%
#       { if (!is.null(class_target)) dplyr::filter(., Class %in% class_target) else . } %>%
#       dplyr::mutate(
#         sample = paste(`Region Subclass`, Sex, sep = "_"),
#         type   = ifelse(Neurotransmitter == "NN", "NN", "N")
#       ) %>%
#       dplyr::select(
#         sample = sample,
#         term_name = GO,
#         value = nlogp,
#         region = Region,
#         sex = Sex,
#         class = Class,
#         type
#       )
    
#     # ---- 构建行分组（行数据） ----
#     go_group_subset <- subset_cfg$filter(go_group)
    
#     # ---- 预置三套过滤集合 ----
#     filtered_sets <- list(
#       common_main = filter_col_and_row(go_results_subset, go_group_subset,
#                                        n_cells = 3, n_genes = 3, row_col = "GO"), # >10
#       common_sup  = filter_col_and_row(go_results_subset, go_group_subset,
#                                        n_cells = 3, n_genes = 1, row_col = "GO"),
#       specific_sup = NULL
#     )
    
#     # ---- 新增：unfiltered（不过滤 GO term 或列）----
#     filtered_sets <- c(
#       list(unfiltered = list(
#         main_data = go_results_subset,
#         row_data  = go_group_subset
#       )),
#       filtered_sets
#     )
    
#     # ---- 基于 common_sup 计算 specific_sup ----
#     common_sup_genes <- unique(filtered_sets$common_sup$row_data$GO)
#     specific_sup_main <- go_results_subset %>%
#       dplyr::filter(!term_name %in% common_sup_genes)
#     specific_sup_row <- go_group_subset %>%
#       dplyr::filter(!GO %in% common_sup_genes)
#     filtered_sets$specific_sup <- list(
#       main_data = specific_sup_main,
#       row_data  = specific_sup_row
#     )
    
#     # ---- 聚类配置 ----
#     cluster_cfgs <- list(
#       list(rows = FALSE, cols = FALSE, tag = "nocluster"),
#       list(rows = FALSE, cols = TRUE,  tag = "cluster_col"),
#       list(rows = TRUE,  cols = FALSE, tag = "cluster_row"),
#       list(rows = TRUE,  cols = TRUE,  tag = "cluster_row_col")
#     )
    
#     # ---- 逐集合、逐聚类配置画图 ----
#     for (set_name in names(filtered_sets)) {
#       data_pair <- filtered_sets[[set_name]]
      
#       # 兼容列表结构或 data.frame 结构
#       main_data <- if (is.list(data_pair) && !is.null(data_pair$main_data)) data_pair$main_data else data_pair
#       row_data  <- if (is.list(data_pair) && !is.null(data_pair$row_data))  data_pair$row_data  else go_group_subset
#       main_data[main_data$type=="N","class"]=gsub("(.* )|(_.*)","",main_data[main_data$type=="N","sample"])
      
#       if (is.null(main_data) || nrow(main_data) == 0) {
#         message("Skipping empty dataset: ", set_name)
#         next
#       }
      
#       base_output <- file.path(
#         output_dir_path,
#         paste0("heatmap_", file_name, "_", subset_tag, "_", set_name)
#       )
      
#       for (cfg in cluster_cfgs) {
#         output_path <- paste0(gsub("_N_","_GABA_",base_output), "_", cfg$tag, ".pdf")
#         message("Generating heatmap with config: ", subset_tag, " + ", set_name, " + ", cfg$tag)
        
#         results <- heatmap_final3(
#           go = main_data,
#           row = row_data,
#           sample_col = "sample",
#           term_col = "term_name",
#           value_col = "value",
#           value_name = "-log10(p_val_adj)",
#           region_col = "region",
#           sex_col = if ("sex" %in% names(main_data)) "sex" else NULL,
#           class_col = if ("class" %in% names(main_data)) "class" else NULL,
#           neurotransmitter_col = NULL,
#           row_go_col = "GO",
#           row_module_col = if (use_row_category) "module" else NULL,
#           row_category_col = NULL,
#           col_order_priority = c("region", "sex", "class"),
#           col_split_col = "region",
#           cutoff = 5,
#           distance = "euclidean",
#           method = "complete", # "ward.D2" 也可
#           logfc_cutoff = NULL,
#           column_gap = grid::unit(3, "mm"),
#           row_gap = grid::unit(3, "mm"),
#           row_names_gp = 6,
#           column_names_gp = 6,
#           title = output_path,
#           cluster_rows = cfg$rows,
#           cluster_columns = cfg$cols,
#           cluster_column_slices = FALSE,
#           cluster_row_slices = FALSE
#         )
#       }
#     }
#   }
  
#   message("Completed processing for: ", go_path)
# }
# for (go_path in go_paths) {
#   message("Processing: ", go_path)
  
#   go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1)
#   go_results_tmp=go_results_tmp[go_results_tmp$Region!="iCTX",]
#   go_results_tmp=go_results_tmp[go_results_tmp$source%in%c("GO:BP","GO:MF"),]
#   go_results_tmp=go_results_tmp[go_results_tmp$Neurotransmitter=="GABA",]
#   #go_results_tmp$padj=p.adjust(go_results_tmp$p_value,method = "BH")
  
  
  
  
#   colnames(go_results_tmp)=gsub("Region.Subclass","Region Subclass",colnames(go_results_tmp))
#   go_results_tmp_unique <- go_results_tmp %>%
#     distinct(`Region Subclass`, Sex, GO, .keep_all = TRUE)
#   if (nrow(go_results_tmp_unique) < nrow(go_results_tmp)) {
#     go_path_unique <- sub("\\.csv$", "_unique.csv", go_path)
#     write.csv(go_results_tmp_unique, go_path_unique, row.names = FALSE)
#   }
  
#   # go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1) %>%
#   #   filter(Sex %in% sexes)
  
#   #go_results <- switch(target_data,
#   #                     "SUS" = go_results_tmp_unique
#   #)
#   go_results=go_results_tmp_unique
  
#   output_dir_path <- dirname(go_path)
#   file_name <- tools::file_path_sans_ext(basename(go_path))
#   #output_dir_path <- file.path(output_dir_path, file_name)
#   output_dir_path <- file.path(output_dir_path, paste0(file_name,"_onlyGOBPMF_GABA_adjP"))
  
#   if (!dir.exists(output_dir_path)) {
#     dir.create(output_dir_path, recursive = TRUE)
#   }
  
#   subsets <- list(
#     list(tag = "N", filter = function(df) df %>% filter(Neurotransmitter != "NN"))
#     #list(tag = "NN",    filter = function(df) df %>% filter(Neurotransmitter == "NN"))
#   )
  
#   group_col <- 'predicted_group' # 'Module.P','predicted_group' or 'final_group'
#   use_row_category <-TRUE
#   class_target <-  NULL # c("Microglia")
#   # go_results$Module.P <- gsub("^G", "M", go_results[[group_col]])
#   go_results <- go_results %>%
#     mutate(!!group_col := ifelse(.data[[group_col]] == "others", "other", .data[[group_col]]))
#   module_levels <- c("G1","G2-1","G2-2","G3-1","G3-2","G3-3","G4","G5","G6","G7","G7-1","G7-2","G8","G8-1","G8-2","G9","G10","G11","G12","other") # c("innate", "overlapped", "adaptive", "inflammatory")
  
#   if (use_row_category) {
#     go_group <- go_results %>%
#       dplyr::select(GO,Neurotransmitter,all_of(group_col)) %>%
#       rename(module = all_of(group_col)) %>%
#       {
#         if (!is.null(module_levels)) {
#           mutate(., module = factor(module, levels = module_levels))
#         } else {
#           mutate(., module = factor(module))
#         }
#       } %>%
#       arrange(module)
#   } else {
#     go_group <- go_results %>%
#       dplyr::select(GO,Neurotransmitter)
#   }
  
  
#   # Generate heatmaps for each subset and clustering configuration
#   # split N and NN for separate heatmaps
#   # 生成不同子集与聚类配置的热图（包含：不过滤GO/列的unfiltered版本）
#   # 提示：确保 dplyr 在最后加载，或统一使用 dplyr:: 前缀（本代码已使用前缀）
#   # library(dplyr)  # 如果你不想每次都写 dplyr::，可把这句放到所有包加载的最后
  
#   # split N and NN for separate heatmaps
#   for (subset_cfg in subsets) {
#     subset_tag <- subset_cfg$tag
#     message("Generating subset: ", subset_tag)
    
#     # ---- 构建主数据（列数据） ----
#     go_results_subset <- subset_cfg$filter(go_results) %>%
#       { if (!is.null(class_target)) dplyr::filter(., Class %in% class_target) else . } %>%
#       dplyr::mutate(
#         sample = paste(`Region Subclass`, Sex, sep = "_"),
#         type   = ifelse(Neurotransmitter == "NN", "NN", "N")
#       ) %>%
#       dplyr::select(
#         sample = sample,
#         term_name = GO,
#         value = nlogp,
#         region = Region,
#         sex = Sex,
#         class = Class,
#         type
#       )
    
#     # ---- 构建行分组（行数据） ----
#     go_group_subset <- subset_cfg$filter(go_group)
    
#     # ---- 预置三套过滤集合 ----
#     filtered_sets <- list(
#       common_main = filter_col_and_row(go_results_subset, go_group_subset,
#                                        n_cells = 3, n_genes = 3, row_col = "GO"), # >10
#       common_sup  = filter_col_and_row(go_results_subset, go_group_subset,
#                                        n_cells = 3, n_genes = 1, row_col = "GO"),
#       specific_sup = NULL
#     )
    
#     # ---- 新增：unfiltered（不过滤 GO term 或列）----
#     filtered_sets <- c(
#       list(unfiltered = list(
#         main_data = go_results_subset,
#         row_data  = go_group_subset
#       )),
#       filtered_sets
#     )
    
#     # ---- 基于 common_sup 计算 specific_sup ----
#     common_sup_genes <- unique(filtered_sets$common_sup$row_data$GO)
#     specific_sup_main <- go_results_subset %>%
#       dplyr::filter(!term_name %in% common_sup_genes)
#     specific_sup_row <- go_group_subset %>%
#       dplyr::filter(!GO %in% common_sup_genes)
#     filtered_sets$specific_sup <- list(
#       main_data = specific_sup_main,
#       row_data  = specific_sup_row
#     )
    
#     # ---- 聚类配置 ----
#     cluster_cfgs <- list(
#       list(rows = FALSE, cols = FALSE, tag = "nocluster"),
#       list(rows = FALSE, cols = TRUE,  tag = "cluster_col"),
#       list(rows = TRUE,  cols = FALSE, tag = "cluster_row"),
#       list(rows = TRUE,  cols = TRUE,  tag = "cluster_row_col")
#     )
    
#     # ---- 逐集合、逐聚类配置画图 ----
#     for (set_name in names(filtered_sets)) {
#       data_pair <- filtered_sets[[set_name]]
      
#       # 兼容列表结构或 data.frame 结构
#       main_data <- if (is.list(data_pair) && !is.null(data_pair$main_data)) data_pair$main_data else data_pair
#       row_data  <- if (is.list(data_pair) && !is.null(data_pair$row_data))  data_pair$row_data  else go_group_subset
#       main_data[main_data$type=="N","class"]=gsub("(.* )|(_.*)","",main_data[main_data$type=="N","sample"])
      
#       if (is.null(main_data) || nrow(main_data) == 0) {
#         message("Skipping empty dataset: ", set_name)
#         next
#       }
      
#       base_output <- file.path(
#         output_dir_path,
#         paste0("heatmap_", file_name, "_", subset_tag, "_", set_name)
#       )
      
#       for (cfg in cluster_cfgs) {
#         output_path <- paste0(gsub("_N_","_GABA_",base_output), "_", cfg$tag, ".pdf")
#         message("Generating heatmap with config: ", subset_tag, " + ", set_name, " + ", cfg$tag)
        
#         results <- heatmap_final2(
#           go = main_data,
#           row = row_data,
#           sample_col = "sample",
#           term_col = "term_name",
#           value_col = "value",
#           value_name = "-log10(p_val_adj)",
#           region_col = "region",
#           sex_col = if ("sex" %in% names(main_data)) "sex" else NULL,
#           class_col = if ("class" %in% names(main_data)) "class" else NULL,
#           neurotransmitter_col = NULL,
#           row_go_col = "GO",
#           row_module_col = if (use_row_category) "module" else NULL,
#           row_category_col = NULL,
#           col_order_priority = c("region", "sex", "class"),
#           col_split_col = "region",
#           cutoff = 5,
#           distance = "euclidean",
#           method = "complete", # "ward.D2" 也可
#           logfc_cutoff = NULL,
#           column_gap = grid::unit(3, "mm"),
#           row_gap = grid::unit(3, "mm"),
#           row_names_gp = 6,
#           column_names_gp = 6,
#           title = output_path,
#           cluster_rows = cfg$rows,
#           cluster_columns = cfg$cols,
#           cluster_column_slices = FALSE,
#           cluster_row_slices = FALSE
#         )
#       }
#     }
#   }
  
#   message("Completed processing for: ", go_path)
# }
# ##Glut adjP####
# for (go_path in go_paths) {
#   message("Processing: ", go_path)
  
#   go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1)
#   go_results_tmp=go_results_tmp[go_results_tmp$Region!="iCTX",]
#   go_results_tmp=go_results_tmp[go_results_tmp$source%in%c("GO:BP","GO:MF"),]
#   go_results_tmp=go_results_tmp[go_results_tmp$Neurotransmitter=="Glut",]
#   #go_results_tmp$padj=p.adjust(go_results_tmp$p_value,method = "BH")
#   #go_results_tmp$nlogp=-log10(go_results_tmp$padj)*sign(go_results_tmp$nlogp)
#   #go_results_tmp$nlog10p=-log10(go_results_tmp$padj)
  
  
  
#   colnames(go_results_tmp)=gsub("Region.Subclass","Region Subclass",colnames(go_results_tmp))
#   go_results_tmp_unique <- go_results_tmp %>%
#     distinct(`Region Subclass`, Sex, GO, .keep_all = TRUE)
#   if (nrow(go_results_tmp_unique) < nrow(go_results_tmp)) {
#     go_path_unique <- sub("\\.csv$", "_unique.csv", go_path)
#     write.csv(go_results_tmp_unique, go_path_unique, row.names = FALSE)
#   }
  
#   # go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1) %>%
#   #   filter(Sex %in% sexes)
  
#   #go_results <- switch(target_data,
#   #                     "SUS" = go_results_tmp_unique
#   #)
#   go_results=go_results_tmp_unique
  
#   output_dir_path <- dirname(go_path)
#   file_name <- tools::file_path_sans_ext(basename(go_path))
#   #output_dir_path <- file.path(output_dir_path, file_name)
#   output_dir_path <- file.path(output_dir_path, paste0(file_name,"_small_onlyGOBPMF_Glut_adjP"))
  
#   if (!dir.exists(output_dir_path)) {
#     dir.create(output_dir_path, recursive = TRUE)
#   }
  
#   subsets <- list(
#     list(tag = "N", filter = function(df) df %>% filter(Neurotransmitter != "NN"))
#     #list(tag = "NN",    filter = function(df) df %>% filter(Neurotransmitter == "NN"))
#   )
  
#   group_col <- 'predicted_group' # 'Module.P','predicted_group' or 'final_group'
#   use_row_category <-TRUE
#   class_target <-  NULL # c("Microglia")
#   # go_results$Module.P <- gsub("^G", "M", go_results[[group_col]])
#   go_results <- go_results %>%
#     mutate(!!group_col := ifelse(.data[[group_col]] == "others", "other", .data[[group_col]]))
#   module_levels <- c("G1","G2-1","G2-2","G3-1","G3-2","G3-3","G4","G5","G6","G7","G7-1","G7-2","G8","G8-1","G8-2","G9","G10","G11","G12","other") # c("innate", "overlapped", "adaptive", "inflammatory")
  
#   if (use_row_category) {
#     go_group <- go_results %>%
#       dplyr::select(GO,Neurotransmitter,all_of(group_col)) %>%
#       rename(module = all_of(group_col)) %>%
#       {
#         if (!is.null(module_levels)) {
#           mutate(., module = factor(module, levels = module_levels))
#         } else {
#           mutate(., module = factor(module))
#         }
#       } %>%
#       arrange(module)
#   } else {
#     go_group <- go_results %>%
#       dplyr::select(GO,Neurotransmitter)
#   }
  
  
#   # Generate heatmaps for each subset and clustering configuration
#   # split N and NN for separate heatmaps
#   # 生成不同子集与聚类配置的热图（包含：不过滤GO/列的unfiltered版本）
#   # 提示：确保 dplyr 在最后加载，或统一使用 dplyr:: 前缀（本代码已使用前缀）
#   # library(dplyr)  # 如果你不想每次都写 dplyr::，可把这句放到所有包加载的最后
  
#   # split N and NN for separate heatmaps
#   for (subset_cfg in subsets) {
#     subset_tag <- subset_cfg$tag
#     message("Generating subset: ", subset_tag)
    
#     # ---- 构建主数据（列数据） ----
#     go_results_subset <- subset_cfg$filter(go_results) %>%
#       { if (!is.null(class_target)) dplyr::filter(., Class %in% class_target) else . } %>%
#       dplyr::mutate(
#         sample = paste(`Region Subclass`, Sex, sep = "_"),
#         type   = ifelse(Neurotransmitter == "NN", "NN", "N")
#       ) %>%
#       dplyr::select(
#         sample = sample,
#         term_name = GO,
#         value = nlogp,
#         region = Region,
#         sex = Sex,
#         class = Class,
#         type
#       )
    
#     # ---- 构建行分组（行数据） ----
#     go_group_subset <- subset_cfg$filter(go_group)
    
#     # ---- 预置三套过滤集合 ----
#     filtered_sets <- list(
#       common_main = filter_col_and_row(go_results_subset, go_group_subset,
#                                        n_cells = 3, n_genes = 3, row_col = "GO"), # >10
#       common_sup  = filter_col_and_row(go_results_subset, go_group_subset,
#                                        n_cells = 3, n_genes = 1, row_col = "GO"),
#       specific_sup = NULL
#     )
    
#     # ---- 新增：unfiltered（不过滤 GO term 或列）----
#     filtered_sets <- c(
#       list(unfiltered = list(
#         main_data = go_results_subset,
#         row_data  = go_group_subset
#       )),
#       filtered_sets
#     )
    
#     # ---- 基于 common_sup 计算 specific_sup ----
#     common_sup_genes <- unique(filtered_sets$common_sup$row_data$GO)
#     specific_sup_main <- go_results_subset %>%
#       dplyr::filter(!term_name %in% common_sup_genes)
#     specific_sup_row <- go_group_subset %>%
#       dplyr::filter(!GO %in% common_sup_genes)
#     filtered_sets$specific_sup <- list(
#       main_data = specific_sup_main,
#       row_data  = specific_sup_row
#     )
    
#     # ---- 聚类配置 ----
#     cluster_cfgs <- list(
#       list(rows = FALSE, cols = FALSE, tag = "nocluster"),
#       list(rows = FALSE, cols = TRUE,  tag = "cluster_col"),
#       list(rows = TRUE,  cols = FALSE, tag = "cluster_row"),
#       list(rows = TRUE,  cols = TRUE,  tag = "cluster_row_col")
#     )
    
#     # ---- 逐集合、逐聚类配置画图 ----
#     for (set_name in names(filtered_sets)) {
#       data_pair <- filtered_sets[[set_name]]
      
#       # 兼容列表结构或 data.frame 结构
#       main_data <- if (is.list(data_pair) && !is.null(data_pair$main_data)) data_pair$main_data else data_pair
#       row_data  <- if (is.list(data_pair) && !is.null(data_pair$row_data))  data_pair$row_data  else go_group_subset
#       main_data[main_data$type=="N","class"]=gsub("(.* )|(_.*)","",main_data[main_data$type=="N","sample"])
      
#       if (is.null(main_data) || nrow(main_data) == 0) {
#         message("Skipping empty dataset: ", set_name)
#         next
#       }
      
#       base_output <- file.path(
#         output_dir_path,
#         paste0("heatmap_", file_name, "_", subset_tag, "_", set_name)
#       )
      
#       for (cfg in cluster_cfgs) {
#         output_path <- paste0(gsub("_N_","_Glut_",base_output), "_", cfg$tag, ".pdf")
#         message("Generating heatmap with config: ", subset_tag, " + ", set_name, " + ", cfg$tag)
        
#         results <- heatmap_final3(
#           go = main_data,
#           row = row_data,
#           sample_col = "sample",
#           term_col = "term_name",
#           value_col = "value",
#           value_name = "-log10(p_val_adj)",
#           region_col = "region",
#           sex_col = if ("sex" %in% names(main_data)) "sex" else NULL,
#           class_col = if ("class" %in% names(main_data)) "class" else NULL,
#           neurotransmitter_col = NULL,
#           row_go_col = "GO",
#           row_module_col = if (use_row_category) "module" else NULL,
#           row_category_col = NULL,
#           col_order_priority = c("region", "sex", "class"),
#           col_split_col = "region",
#           cutoff = 5,
#           distance = "euclidean",
#           method = "complete", # "ward.D2" 也可
#           logfc_cutoff = NULL,
#           column_gap = grid::unit(3, "mm"),
#           row_gap = grid::unit(3, "mm"),
#           row_names_gp = 6,
#           column_names_gp = 6,
#           title = output_path,
#           cluster_rows = cfg$rows,
#           cluster_columns = cfg$cols,
#           cluster_column_slices = FALSE,
#           cluster_row_slices = FALSE
#         )
#       }
#     }
#   }
  
#   message("Completed processing for: ", go_path)
# }
# for (go_path in go_paths) {
#   message("Processing: ", go_path)
  
#   go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1)
#   go_results_tmp=go_results_tmp[go_results_tmp$Region!="iCTX",]
#   go_results_tmp=go_results_tmp[go_results_tmp$source%in%c("GO:BP","GO:MF"),]
#   go_results_tmp=go_results_tmp[go_results_tmp$Neurotransmitter=="Glut",]
#   #go_results_tmp$padj=p.adjust(go_results_tmp$p_value,method = "BH")
#   #go_results_tmp$nlogp=-log10(go_results_tmp$padj)*sign(go_results_tmp$nlogp)
#   #go_results_tmp$nlog10p=-log10(go_results_tmp$padj)
  
  
  
#   colnames(go_results_tmp)=gsub("Region.Subclass","Region Subclass",colnames(go_results_tmp))
#   go_results_tmp_unique <- go_results_tmp %>%
#     distinct(`Region Subclass`, Sex, GO, .keep_all = TRUE)
#   if (nrow(go_results_tmp_unique) < nrow(go_results_tmp)) {
#     go_path_unique <- sub("\\.csv$", "_unique.csv", go_path)
#     write.csv(go_results_tmp_unique, go_path_unique, row.names = FALSE)
#   }
  
#   # go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1) %>%
#   #   filter(Sex %in% sexes)
  
#   #go_results <- switch(target_data,
#   #                     "SUS" = go_results_tmp_unique
#   #)
#   go_results=go_results_tmp_unique
  
#   output_dir_path <- dirname(go_path)
#   file_name <- tools::file_path_sans_ext(basename(go_path))
#   #output_dir_path <- file.path(output_dir_path, file_name)
#   output_dir_path <- file.path(output_dir_path, paste0(file_name,"_onlyGOBPMF_Glut_adjP"))
  
#   if (!dir.exists(output_dir_path)) {
#     dir.create(output_dir_path, recursive = TRUE)
#   }
  
#   subsets <- list(
#     list(tag = "N", filter = function(df) df %>% filter(Neurotransmitter != "NN"))
#     #list(tag = "NN",    filter = function(df) df %>% filter(Neurotransmitter == "NN"))
#   )
  
#   group_col <- 'predicted_group' # 'Module.P','predicted_group' or 'final_group'
#   use_row_category <-TRUE
#   class_target <-  NULL # c("Microglia")
#   # go_results$Module.P <- gsub("^G", "M", go_results[[group_col]])
#   go_results <- go_results %>%
#     mutate(!!group_col := ifelse(.data[[group_col]] == "others", "other", .data[[group_col]]))
#   module_levels <- c("G1","G2-1","G2-2","G3-1","G3-2","G3-3","G4","G5","G6","G7","G7-1","G7-2","G8","G8-1","G8-2","G9","G10","G11","G12","other") # c("innate", "overlapped", "adaptive", "inflammatory")
  
#   if (use_row_category) {
#     go_group <- go_results %>%
#       dplyr::select(GO,Neurotransmitter,all_of(group_col)) %>%
#       rename(module = all_of(group_col)) %>%
#       {
#         if (!is.null(module_levels)) {
#           mutate(., module = factor(module, levels = module_levels))
#         } else {
#           mutate(., module = factor(module))
#         }
#       } %>%
#       arrange(module)
#   } else {
#     go_group <- go_results %>%
#       dplyr::select(GO,Neurotransmitter)
#   }
  
  
#   # Generate heatmaps for each subset and clustering configuration
#   # split N and NN for separate heatmaps
#   # 生成不同子集与聚类配置的热图（包含：不过滤GO/列的unfiltered版本）
#   # 提示：确保 dplyr 在最后加载，或统一使用 dplyr:: 前缀（本代码已使用前缀）
#   # library(dplyr)  # 如果你不想每次都写 dplyr::，可把这句放到所有包加载的最后
  
#   # split N and NN for separate heatmaps
#   for (subset_cfg in subsets) {
#     subset_tag <- subset_cfg$tag
#     message("Generating subset: ", subset_tag)
    
#     # ---- 构建主数据（列数据） ----
#     go_results_subset <- subset_cfg$filter(go_results) %>%
#       { if (!is.null(class_target)) dplyr::filter(., Class %in% class_target) else . } %>%
#       dplyr::mutate(
#         sample = paste(`Region Subclass`, Sex, sep = "_"),
#         type   = ifelse(Neurotransmitter == "NN", "NN", "N")
#       ) %>%
#       dplyr::select(
#         sample = sample,
#         term_name = GO,
#         value = nlogp,
#         region = Region,
#         sex = Sex,
#         class = Class,
#         type
#       )
    
#     # ---- 构建行分组（行数据） ----
#     go_group_subset <- subset_cfg$filter(go_group)
    
#     # ---- 预置三套过滤集合 ----
#     filtered_sets <- list(
#       common_main = filter_col_and_row(go_results_subset, go_group_subset,
#                                        n_cells = 3, n_genes = 3, row_col = "GO"), # >10
#       common_sup  = filter_col_and_row(go_results_subset, go_group_subset,
#                                        n_cells = 3, n_genes = 1, row_col = "GO"),
#       specific_sup = NULL
#     )
    
#     # ---- 新增：unfiltered（不过滤 GO term 或列）----
#     filtered_sets <- c(
#       list(unfiltered = list(
#         main_data = go_results_subset,
#         row_data  = go_group_subset
#       )),
#       filtered_sets
#     )
    
#     # ---- 基于 common_sup 计算 specific_sup ----
#     common_sup_genes <- unique(filtered_sets$common_sup$row_data$GO)
#     specific_sup_main <- go_results_subset %>%
#       dplyr::filter(!term_name %in% common_sup_genes)
#     specific_sup_row <- go_group_subset %>%
#       dplyr::filter(!GO %in% common_sup_genes)
#     filtered_sets$specific_sup <- list(
#       main_data = specific_sup_main,
#       row_data  = specific_sup_row
#     )
    
#     # ---- 聚类配置 ----
#     cluster_cfgs <- list(
#       list(rows = FALSE, cols = FALSE, tag = "nocluster"),
#       list(rows = FALSE, cols = TRUE,  tag = "cluster_col"),
#       list(rows = TRUE,  cols = FALSE, tag = "cluster_row"),
#       list(rows = TRUE,  cols = TRUE,  tag = "cluster_row_col")
#     )
    
#     # ---- 逐集合、逐聚类配置画图 ----
#     for (set_name in names(filtered_sets)) {
#       data_pair <- filtered_sets[[set_name]]
      
#       # 兼容列表结构或 data.frame 结构
#       main_data <- if (is.list(data_pair) && !is.null(data_pair$main_data)) data_pair$main_data else data_pair
#       row_data  <- if (is.list(data_pair) && !is.null(data_pair$row_data))  data_pair$row_data  else go_group_subset
#       main_data[main_data$type=="N","class"]=gsub("(.* )|(_.*)","",main_data[main_data$type=="N","sample"])
      
#       if (is.null(main_data) || nrow(main_data) == 0) {
#         message("Skipping empty dataset: ", set_name)
#         next
#       }
      
#       base_output <- file.path(
#         output_dir_path,
#         paste0("heatmap_", file_name, "_", subset_tag, "_", set_name)
#       )
      
#       for (cfg in cluster_cfgs) {
#         output_path <- paste0(gsub("_N_","_Glut_",base_output), "_", cfg$tag, ".pdf")
#         message("Generating heatmap with config: ", subset_tag, " + ", set_name, " + ", cfg$tag)
        
#         results <- heatmap_final2(
#           go = main_data,
#           row = row_data,
#           sample_col = "sample",
#           term_col = "term_name",
#           value_col = "value",
#           value_name = "-log10(p_val_adj)",
#           region_col = "region",
#           sex_col = if ("sex" %in% names(main_data)) "sex" else NULL,
#           class_col = if ("class" %in% names(main_data)) "class" else NULL,
#           neurotransmitter_col = NULL,
#           row_go_col = "GO",
#           row_module_col = if (use_row_category) "module" else NULL,
#           row_category_col = NULL,
#           col_order_priority = c("region", "sex", "class"),
#           col_split_col = "region",
#           cutoff = 5,
#           distance = "euclidean",
#           method = "complete", # "ward.D2" 也可
#           logfc_cutoff = NULL,
#           column_gap = grid::unit(3, "mm"),
#           row_gap = grid::unit(3, "mm"),
#           row_names_gp = 6,
#           column_names_gp = 6,
#           title = output_path,
#           cluster_rows = cfg$rows,
#           cluster_columns = cfg$cols,
#           cluster_column_slices = FALSE,
#           cluster_row_slices = FALSE
#         )
#       }
#     }
#   }
  
#   message("Completed processing for: ", go_path)
# }


#####################

####for GSEA##############
##GABA NES#########################################
for (go_path in go_paths) {
  message("Processing: ", go_path)
  
  go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1)
  go_results_tmp=go_results_tmp[go_results_tmp$Region!="iCTX",]
  go_results_tmp=go_results_tmp[go_results_tmp$source%in%c("GO:BP","GO:MF"),]
  go_results_tmp=go_results_tmp[go_results_tmp$Neurotransmitter=="GABA",]
  #go_results_tmp$padj=p.adjust(go_results_tmp$p_value,method = "BH")
 # go_results_tmp$nlogp=-log10(go_results_tmp$padj)*sign(go_results_tmp$nlogp)
#  go_results_tmp$nlog10p=-log10(go_results_tmp$padj)
  go_results_tmp$nlogp=go_results_tmp$NES
  go_results_tmp=go_results_tmp[go_results_tmp$p_value<0.01,]
  colnames(go_results_tmp)=gsub("Region.Subclass","Region Subclass",colnames(go_results_tmp))
  go_results_tmp_unique <- go_results_tmp %>%
    distinct(`Region Subclass`, Sex, GO, .keep_all = TRUE)
  if (nrow(go_results_tmp_unique) < nrow(go_results_tmp)) {
    go_path_unique <- sub("\\.csv$", "_unique.csv", go_path)
    write.csv(go_results_tmp_unique, go_path_unique, row.names = FALSE)
  }
  
  # go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1) %>%
  #   filter(Sex %in% sexes)
  
  #go_results <- switch(target_data,
  #                     "SUS" = go_results_tmp_unique
  #)
  go_results=go_results_tmp_unique
  
  output_dir_path <- dirname(go_path)
  file_name <- tools::file_path_sans_ext(basename(go_path))
  #output_dir_path <- file.path(output_dir_path, file_name)
  output_dir_path <- file.path(output_dir_path, paste0(file_name,"_small_onlyGOBPMF_GABA_NES"))
  
  if (!dir.exists(output_dir_path)) {
    dir.create(output_dir_path, recursive = TRUE)
  }
  
  subsets <- list(
    list(tag = "N", filter = function(df) df %>% filter(Neurotransmitter != "NN"))
    #list(tag = "NN",    filter = function(df) df %>% filter(Neurotransmitter == "NN"))
  )
  
  group_col <- 'predicted_group' # 'Module.P','predicted_group' or 'final_group'
  use_row_category <-TRUE
  class_target <-  NULL # c("Microglia")
  # go_results$Module.P <- gsub("^G", "M", go_results[[group_col]])
  go_results <- go_results %>%
    mutate(!!group_col := ifelse(.data[[group_col]] == "others", "other", .data[[group_col]]))
  module_levels <- c("G1","G2-1","G2-2","G3-1","G3-2","G3-3","G4","G5","G6","G7","G7-1","G7-2","G8","G8-1","G8-2","G9","G10","G11","G12","other") # c("innate", "overlapped", "adaptive", "inflammatory")
  
  if (use_row_category) {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter,all_of(group_col)) %>%
      rename(module = all_of(group_col)) %>%
      {
        if (!is.null(module_levels)) {
          mutate(., module = factor(module, levels = module_levels))
        } else {
          mutate(., module = factor(module))
        }
      } %>%
      arrange(module)
  } else {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter)
  }
  
  
  # Generate heatmaps for each subset and clustering configuration
  # split N and NN for separate heatmaps
  # 生成不同子集与聚类配置的热图（包含：不过滤GO/列的unfiltered版本）
  # 提示：确保 dplyr 在最后加载，或统一使用 dplyr:: 前缀（本代码已使用前缀）
  # library(dplyr)  # 如果你不想每次都写 dplyr::，可把这句放到所有包加载的最后
  
  # split N and NN for separate heatmaps
  for (subset_cfg in subsets) {
    subset_tag <- subset_cfg$tag
    message("Generating subset: ", subset_tag)
    
    # ---- 构建主数据（列数据） ----
    go_results_subset <- subset_cfg$filter(go_results) %>%
      { if (!is.null(class_target)) dplyr::filter(., Class %in% class_target) else . } %>%
      dplyr::mutate(
        sample = paste(`Region Subclass`, Sex, sep = "_"),
        type   = ifelse(Neurotransmitter == "NN", "NN", "N")
      ) %>%
      dplyr::select(
        sample = sample,
        term_name = GO,
        value = nlogp,
        region = Region,
        sex = Sex,
        class = Class,
        type
      )
    
    # ---- 构建行分组（行数据） ----
    go_group_subset <- subset_cfg$filter(go_group)
    
    # ---- 预置三套过滤集合 ----
    filtered_sets <- list(
      common_main = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 3, row_col = "GO"), # >10
      common_sup  = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 1, row_col = "GO"),
      specific_sup = NULL
    )
    
    # ---- 新增：unfiltered（不过滤 GO term 或列）----
    filtered_sets <- c(
      list(unfiltered = list(
        main_data = go_results_subset,
        row_data  = go_group_subset
      )),
      filtered_sets
    )
    
    # ---- 基于 common_sup 计算 specific_sup ----
    common_sup_genes <- unique(filtered_sets$common_sup$row_data$GO)
    specific_sup_main <- go_results_subset %>%
      dplyr::filter(!term_name %in% common_sup_genes)
    specific_sup_row <- go_group_subset %>%
      dplyr::filter(!GO %in% common_sup_genes)
    filtered_sets$specific_sup <- list(
      main_data = specific_sup_main,
      row_data  = specific_sup_row
    )
    
    # ---- 聚类配置 ----
    cluster_cfgs <- list(
      list(rows = FALSE, cols = FALSE, tag = "nocluster"),
      list(rows = FALSE, cols = TRUE,  tag = "cluster_col"),
      list(rows = TRUE,  cols = FALSE, tag = "cluster_row"),
      list(rows = TRUE,  cols = TRUE,  tag = "cluster_row_col")
    )
    
    # ---- 逐集合、逐聚类配置画图 ----
    for (set_name in names(filtered_sets)) {
      data_pair <- filtered_sets[[set_name]]
      
      # 兼容列表结构或 data.frame 结构
      main_data <- if (is.list(data_pair) && !is.null(data_pair$main_data)) data_pair$main_data else data_pair
      row_data  <- if (is.list(data_pair) && !is.null(data_pair$row_data))  data_pair$row_data  else go_group_subset
      main_data[main_data$type=="N","class"]=gsub("(.* )|(_.*)","",main_data[main_data$type=="N","sample"])
      
      if (is.null(main_data) || nrow(main_data) == 0) {
        message("Skipping empty dataset: ", set_name)
        next
      }
      
      base_output <- file.path(
        output_dir_path,
        paste0("heatmap_", file_name, "_", subset_tag, "_", set_name)
      )
      
      for (cfg in cluster_cfgs) {
        output_path <- paste0(gsub("_N_","_GABA_",base_output), "_", cfg$tag, ".pdf")
        message("Generating heatmap with config: ", subset_tag, " + ", set_name, " + ", cfg$tag)
        
        results <- heatmap_final3(
          go = main_data,
          row = row_data,
          sample_col = "sample",
          term_col = "term_name",
          value_col = "value",
          value_name = "NES",
          region_col = "region",
          sex_col = if ("sex" %in% names(main_data)) "sex" else NULL,
          class_col = if ("class" %in% names(main_data)) "class" else NULL,
          neurotransmitter_col = NULL,
          row_go_col = "GO",
          row_module_col = if (use_row_category) "module" else NULL,
          row_category_col = NULL,
          col_order_priority = c("region", "sex", "class"),
          col_split_col = "region",
          cutoff = 5,
          distance = "euclidean",
          method = "complete", # "ward.D2" 也可
          logfc_cutoff = NULL,
          column_gap = grid::unit(3, "mm"),
          row_gap = grid::unit(3, "mm"),
          row_names_gp = 6,
          column_names_gp = 6,
          title = output_path,
          cluster_rows = cfg$rows,
          cluster_columns = cfg$cols,
          cluster_column_slices = FALSE,
          cluster_row_slices = FALSE
        )
      }
    }
  }
  
  message("Completed processing for: ", go_path)
}
for (go_path in go_paths) {
  message("Processing: ", go_path)
  
  go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1)
  go_results_tmp=go_results_tmp[go_results_tmp$Region!="iCTX",]
  go_results_tmp=go_results_tmp[go_results_tmp$source%in%c("GO:BP","GO:MF"),]
  go_results_tmp=go_results_tmp[go_results_tmp$Neurotransmitter=="GABA",]
  #go_results_tmp$padj=p.adjust(go_results_tmp$p_value,method = "BH")
  go_results_tmp$nlogp=go_results_tmp$NES
  
  
  go_results_tmp=go_results_tmp[go_results_tmp$p_value<0.01,]
  colnames(go_results_tmp)=gsub("Region.Subclass","Region Subclass",colnames(go_results_tmp))
  go_results_tmp_unique <- go_results_tmp %>%
    distinct(`Region Subclass`, Sex, GO, .keep_all = TRUE)
  if (nrow(go_results_tmp_unique) < nrow(go_results_tmp)) {
    go_path_unique <- sub("\\.csv$", "_unique.csv", go_path)
    write.csv(go_results_tmp_unique, go_path_unique, row.names = FALSE)
  }
  
  # go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1) %>%
  #   filter(Sex %in% sexes)
  
  #go_results <- switch(target_data,
  #                     "SUS" = go_results_tmp_unique
  #)
  go_results=go_results_tmp_unique
  
  output_dir_path <- dirname(go_path)
  file_name <- tools::file_path_sans_ext(basename(go_path))
  #output_dir_path <- file.path(output_dir_path, file_name)
  output_dir_path <- file.path(output_dir_path, paste0(file_name,"_onlyGOBPMF_GABA_NES"))
  
  if (!dir.exists(output_dir_path)) {
    dir.create(output_dir_path, recursive = TRUE)
  }
  
  subsets <- list(
    list(tag = "N", filter = function(df) df %>% filter(Neurotransmitter != "NN"))
    #list(tag = "NN",    filter = function(df) df %>% filter(Neurotransmitter == "NN"))
  )
  
  group_col <- 'predicted_group' # 'Module.P','predicted_group' or 'final_group'
  use_row_category <-TRUE
  class_target <-  NULL # c("Microglia")
  # go_results$Module.P <- gsub("^G", "M", go_results[[group_col]])
  go_results <- go_results %>%
    mutate(!!group_col := ifelse(.data[[group_col]] == "others", "other", .data[[group_col]]))
  module_levels <- c("G1","G2-1","G2-2","G3-1","G3-2","G3-3","G4","G5","G6","G7","G7-1","G7-2","G8","G8-1","G8-2","G9","G10","G11","G12","other") # c("innate", "overlapped", "adaptive", "inflammatory")
  
  if (use_row_category) {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter,all_of(group_col)) %>%
      rename(module = all_of(group_col)) %>%
      {
        if (!is.null(module_levels)) {
          mutate(., module = factor(module, levels = module_levels))
        } else {
          mutate(., module = factor(module))
        }
      } %>%
      arrange(module)
  } else {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter)
  }
  
  
  # Generate heatmaps for each subset and clustering configuration
  # split N and NN for separate heatmaps
  # 生成不同子集与聚类配置的热图（包含：不过滤GO/列的unfiltered版本）
  # 提示：确保 dplyr 在最后加载，或统一使用 dplyr:: 前缀（本代码已使用前缀）
  # library(dplyr)  # 如果你不想每次都写 dplyr::，可把这句放到所有包加载的最后
  
  # split N and NN for separate heatmaps
  for (subset_cfg in subsets) {
    subset_tag <- subset_cfg$tag
    message("Generating subset: ", subset_tag)
    
    # ---- 构建主数据（列数据） ----
    go_results_subset <- subset_cfg$filter(go_results) %>%
      { if (!is.null(class_target)) dplyr::filter(., Class %in% class_target) else . } %>%
      dplyr::mutate(
        sample = paste(`Region Subclass`, Sex, sep = "_"),
        type   = ifelse(Neurotransmitter == "NN", "NN", "N")
      ) %>%
      dplyr::select(
        sample = sample,
        term_name = GO,
        value = nlogp,
        region = Region,
        sex = Sex,
        class = Class,
        type
      )
    
    # ---- 构建行分组（行数据） ----
    go_group_subset <- subset_cfg$filter(go_group)
    
    # ---- 预置三套过滤集合 ----
    filtered_sets <- list(
      common_main = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 3, row_col = "GO"), # >10
      common_sup  = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 1, row_col = "GO"),
      specific_sup = NULL
    )
    
    # ---- 新增：unfiltered（不过滤 GO term 或列）----
    filtered_sets <- c(
      list(unfiltered = list(
        main_data = go_results_subset,
        row_data  = go_group_subset
      )),
      filtered_sets
    )
    
    # ---- 基于 common_sup 计算 specific_sup ----
    common_sup_genes <- unique(filtered_sets$common_sup$row_data$GO)
    specific_sup_main <- go_results_subset %>%
      dplyr::filter(!term_name %in% common_sup_genes)
    specific_sup_row <- go_group_subset %>%
      dplyr::filter(!GO %in% common_sup_genes)
    filtered_sets$specific_sup <- list(
      main_data = specific_sup_main,
      row_data  = specific_sup_row
    )
    
    # ---- 聚类配置 ----
    cluster_cfgs <- list(
      list(rows = FALSE, cols = FALSE, tag = "nocluster"),
      list(rows = FALSE, cols = TRUE,  tag = "cluster_col"),
      list(rows = TRUE,  cols = FALSE, tag = "cluster_row"),
      list(rows = TRUE,  cols = TRUE,  tag = "cluster_row_col")
    )
    
    # ---- 逐集合、逐聚类配置画图 ----
    for (set_name in names(filtered_sets)) {
      data_pair <- filtered_sets[[set_name]]
      
      # 兼容列表结构或 data.frame 结构
      main_data <- if (is.list(data_pair) && !is.null(data_pair$main_data)) data_pair$main_data else data_pair
      row_data  <- if (is.list(data_pair) && !is.null(data_pair$row_data))  data_pair$row_data  else go_group_subset
      main_data[main_data$type=="N","class"]=gsub("(.* )|(_.*)","",main_data[main_data$type=="N","sample"])
      
      if (is.null(main_data) || nrow(main_data) == 0) {
        message("Skipping empty dataset: ", set_name)
        next
      }
      
      base_output <- file.path(
        output_dir_path,
        paste0("heatmap_", file_name, "_", subset_tag, "_", set_name)
      )
      
      for (cfg in cluster_cfgs) {
        output_path <- paste0(gsub("_N_","_GABA_",base_output), "_", cfg$tag, ".pdf")
        message("Generating heatmap with config: ", subset_tag, " + ", set_name, " + ", cfg$tag)
        
        results <- heatmap_final2(
          go = main_data,
          row = row_data,
          sample_col = "sample",
          term_col = "term_name",
          value_col = "value",
          value_name = "NES",
          region_col = "region",
          sex_col = if ("sex" %in% names(main_data)) "sex" else NULL,
          class_col = if ("class" %in% names(main_data)) "class" else NULL,
          neurotransmitter_col = NULL,
          row_go_col = "GO",
          row_module_col = if (use_row_category) "module" else NULL,
          row_category_col = NULL,
          col_order_priority = c("region", "sex", "class"),
          col_split_col = "region",
          cutoff = 5,
          distance = "euclidean",
          method = "complete", # "ward.D2" 也可
          logfc_cutoff = NULL,
          column_gap = grid::unit(3, "mm"),
          row_gap = grid::unit(3, "mm"),
          row_names_gp = 6,
          column_names_gp = 6,
          title = output_path,
          cluster_rows = cfg$rows,
          cluster_columns = cfg$cols,
          cluster_column_slices = FALSE,
          cluster_row_slices = FALSE
        )
      }
    }
  }
  
  message("Completed processing for: ", go_path)
}
##Glut NES ####
for (go_path in go_paths) {
  message("Processing: ", go_path)
  
  go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1)
  go_results_tmp=go_results_tmp[go_results_tmp$Region!="iCTX",]
  go_results_tmp=go_results_tmp[go_results_tmp$source%in%c("GO:BP","GO:MF"),]
  go_results_tmp=go_results_tmp[go_results_tmp$Neurotransmitter=="Glut",]
  #go_results_tmp$padj=p.adjust(go_results_tmp$p_value,method = "BH")
  #go_results_tmp$nlogp=-log10(go_results_tmp$padj)*sign(go_results_tmp$nlogp)
  #go_results_tmp$nlog10p=-log10(go_results_tmp$padj)
  go_results_tmp$nlogp=go_results_tmp$NES
  
  go_results_tmp=go_results_tmp[go_results_tmp$p_value<0.01,]
  colnames(go_results_tmp)=gsub("Region.Subclass","Region Subclass",colnames(go_results_tmp))
  go_results_tmp_unique <- go_results_tmp %>%
    distinct(`Region Subclass`, Sex, GO, .keep_all = TRUE)
  if (nrow(go_results_tmp_unique) < nrow(go_results_tmp)) {
    go_path_unique <- sub("\\.csv$", "_unique.csv", go_path)
    write.csv(go_results_tmp_unique, go_path_unique, row.names = FALSE)
  }
  
  # go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1) %>%
  #   filter(Sex %in% sexes)
  
  #go_results <- switch(target_data,
  #                     "SUS" = go_results_tmp_unique
  #)
  go_results=go_results_tmp_unique
  
  output_dir_path <- dirname(go_path)
  file_name <- tools::file_path_sans_ext(basename(go_path))
  #output_dir_path <- file.path(output_dir_path, file_name)
  output_dir_path <- file.path(output_dir_path, paste0(file_name,"_small_onlyGOBPMF_Glut_NES"))
  
  if (!dir.exists(output_dir_path)) {
    dir.create(output_dir_path, recursive = TRUE)
  }
  
  subsets <- list(
    list(tag = "N", filter = function(df) df %>% filter(Neurotransmitter != "NN"))
    #list(tag = "NN",    filter = function(df) df %>% filter(Neurotransmitter == "NN"))
  )
  
  group_col <- 'predicted_group' # 'Module.P','predicted_group' or 'final_group'
  use_row_category <-TRUE
  class_target <-  NULL # c("Microglia")
  # go_results$Module.P <- gsub("^G", "M", go_results[[group_col]])
  go_results <- go_results %>%
    mutate(!!group_col := ifelse(.data[[group_col]] == "others", "other", .data[[group_col]]))
  module_levels <- c("G1","G2-1","G2-2","G3-1","G3-2","G3-3","G4","G5","G6","G7","G7-1","G7-2","G8","G8-1","G8-2","G9","G10","G11","G12","other") # c("innate", "overlapped", "adaptive", "inflammatory")
  
  if (use_row_category) {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter,all_of(group_col)) %>%
      rename(module = all_of(group_col)) %>%
      {
        if (!is.null(module_levels)) {
          mutate(., module = factor(module, levels = module_levels))
        } else {
          mutate(., module = factor(module))
        }
      } %>%
      arrange(module)
  } else {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter)
  }
  
  
  # Generate heatmaps for each subset and clustering configuration
  # split N and NN for separate heatmaps
  # 生成不同子集与聚类配置的热图（包含：不过滤GO/列的unfiltered版本）
  # 提示：确保 dplyr 在最后加载，或统一使用 dplyr:: 前缀（本代码已使用前缀）
  # library(dplyr)  # 如果你不想每次都写 dplyr::，可把这句放到所有包加载的最后
  
  # split N and NN for separate heatmaps
  for (subset_cfg in subsets) {
    subset_tag <- subset_cfg$tag
    message("Generating subset: ", subset_tag)
    
    # ---- 构建主数据（列数据） ----
    go_results_subset <- subset_cfg$filter(go_results) %>%
      { if (!is.null(class_target)) dplyr::filter(., Class %in% class_target) else . } %>%
      dplyr::mutate(
        sample = paste(`Region Subclass`, Sex, sep = "_"),
        type   = ifelse(Neurotransmitter == "NN", "NN", "N")
      ) %>%
      dplyr::select(
        sample = sample,
        term_name = GO,
        value = nlogp,
        region = Region,
        sex = Sex,
        class = Class,
        type
      )
    
    # ---- 构建行分组（行数据） ----
    go_group_subset <- subset_cfg$filter(go_group)
    
    # ---- 预置三套过滤集合 ----
    filtered_sets <- list(
      common_main = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 3, row_col = "GO"), # >10
      common_sup  = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 1, row_col = "GO"),
      specific_sup = NULL
    )
    
    # ---- 新增：unfiltered（不过滤 GO term 或列）----
    filtered_sets <- c(
      list(unfiltered = list(
        main_data = go_results_subset,
        row_data  = go_group_subset
      )),
      filtered_sets
    )
    
    # ---- 基于 common_sup 计算 specific_sup ----
    common_sup_genes <- unique(filtered_sets$common_sup$row_data$GO)
    specific_sup_main <- go_results_subset %>%
      dplyr::filter(!term_name %in% common_sup_genes)
    specific_sup_row <- go_group_subset %>%
      dplyr::filter(!GO %in% common_sup_genes)
    filtered_sets$specific_sup <- list(
      main_data = specific_sup_main,
      row_data  = specific_sup_row
    )
    
    # ---- 聚类配置 ----
    cluster_cfgs <- list(
      list(rows = FALSE, cols = FALSE, tag = "nocluster"),
      list(rows = FALSE, cols = TRUE,  tag = "cluster_col"),
      list(rows = TRUE,  cols = FALSE, tag = "cluster_row"),
      list(rows = TRUE,  cols = TRUE,  tag = "cluster_row_col")
    )
    
    # ---- 逐集合、逐聚类配置画图 ----
    for (set_name in names(filtered_sets)) {
      data_pair <- filtered_sets[[set_name]]
      
      # 兼容列表结构或 data.frame 结构
      main_data <- if (is.list(data_pair) && !is.null(data_pair$main_data)) data_pair$main_data else data_pair
      row_data  <- if (is.list(data_pair) && !is.null(data_pair$row_data))  data_pair$row_data  else go_group_subset
      main_data[main_data$type=="N","class"]=gsub("(.* )|(_.*)","",main_data[main_data$type=="N","sample"])
      
      if (is.null(main_data) || nrow(main_data) == 0) {
        message("Skipping empty dataset: ", set_name)
        next
      }
      
      base_output <- file.path(
        output_dir_path,
        paste0("heatmap_", file_name, "_", subset_tag, "_", set_name)
      )
      
      for (cfg in cluster_cfgs) {
        output_path <- paste0(gsub("_N_","_Glut_",base_output), "_", cfg$tag, ".pdf")
        message("Generating heatmap with config: ", subset_tag, " + ", set_name, " + ", cfg$tag)
        
        results <- heatmap_final3(
          go = main_data,
          row = row_data,
          sample_col = "sample",
          term_col = "term_name",
          value_col = "value",
          value_name = "NES",
          region_col = "region",
          sex_col = if ("sex" %in% names(main_data)) "sex" else NULL,
          class_col = if ("class" %in% names(main_data)) "class" else NULL,
          neurotransmitter_col = NULL,
          row_go_col = "GO",
          row_module_col = if (use_row_category) "module" else NULL,
          row_category_col = NULL,
          col_order_priority = c("region", "sex", "class"),
          col_split_col = "region",
          cutoff = 5,
          distance = "euclidean",
          method = "complete", # "ward.D2" 也可
          logfc_cutoff = NULL,
          column_gap = grid::unit(3, "mm"),
          row_gap = grid::unit(3, "mm"),
          row_names_gp = 6,
          column_names_gp = 6,
          title = output_path,
          cluster_rows = cfg$rows,
          cluster_columns = cfg$cols,
          cluster_column_slices = FALSE,
          cluster_row_slices = FALSE
        )
      }
    }
  }
  
  message("Completed processing for: ", go_path)
}
for (go_path in go_paths) {
  message("Processing: ", go_path)
  
  go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1)
  go_results_tmp=go_results_tmp[go_results_tmp$Region!="iCTX",]
  go_results_tmp=go_results_tmp[go_results_tmp$source%in%c("GO:BP","GO:MF"),]
  go_results_tmp=go_results_tmp[go_results_tmp$Neurotransmitter=="Glut",]
  #go_results_tmp$padj=p.adjust(go_results_tmp$p_value,method = "BH")
  #go_results_tmp$nlogp=-log10(go_results_tmp$padj)*sign(go_results_tmp$nlogp)
  #go_results_tmp$nlog10p=-log10(go_results_tmp$padj)
  go_results_tmp$nlogp=go_results_tmp$NES
  
  go_results_tmp=go_results_tmp[go_results_tmp$p_value<0.01,]
  colnames(go_results_tmp)=gsub("Region.Subclass","Region Subclass",colnames(go_results_tmp))
  go_results_tmp_unique <- go_results_tmp %>%
    distinct(`Region Subclass`, Sex, GO, .keep_all = TRUE)
  if (nrow(go_results_tmp_unique) < nrow(go_results_tmp)) {
    go_path_unique <- sub("\\.csv$", "_unique.csv", go_path)
    write.csv(go_results_tmp_unique, go_path_unique, row.names = FALSE)
  }
  
  # go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1) %>%
  #   filter(Sex %in% sexes)
  
  #go_results <- switch(target_data,
  #                     "SUS" = go_results_tmp_unique
  #)
  go_results=go_results_tmp_unique
  
  output_dir_path <- dirname(go_path)
  file_name <- tools::file_path_sans_ext(basename(go_path))
  #output_dir_path <- file.path(output_dir_path, file_name)
  output_dir_path <- file.path(output_dir_path, paste0(file_name,"_onlyGOBPMF_Glut_NES"))
  
  if (!dir.exists(output_dir_path)) {
    dir.create(output_dir_path, recursive = TRUE)
  }
  
  subsets <- list(
    list(tag = "N", filter = function(df) df %>% filter(Neurotransmitter != "NN"))
    #list(tag = "NN",    filter = function(df) df %>% filter(Neurotransmitter == "NN"))
  )
  
  group_col <- 'predicted_group' # 'Module.P','predicted_group' or 'final_group'
  use_row_category <-TRUE
  class_target <-  NULL # c("Microglia")
  # go_results$Module.P <- gsub("^G", "M", go_results[[group_col]])
  go_results <- go_results %>%
    mutate(!!group_col := ifelse(.data[[group_col]] == "others", "other", .data[[group_col]]))
  module_levels <- c("G1","G2-1","G2-2","G3-1","G3-2","G3-3","G4","G5","G6","G7","G7-1","G7-2","G8","G8-1","G8-2","G9","G10","G11","G12","other") # c("innate", "overlapped", "adaptive", "inflammatory")
  
  if (use_row_category) {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter,all_of(group_col)) %>%
      rename(module = all_of(group_col)) %>%
      {
        if (!is.null(module_levels)) {
          mutate(., module = factor(module, levels = module_levels))
        } else {
          mutate(., module = factor(module))
        }
      } %>%
      arrange(module)
  } else {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter)
  }
  
  
  # Generate heatmaps for each subset and clustering configuration
  # split N and NN for separate heatmaps
  # 生成不同子集与聚类配置的热图（包含：不过滤GO/列的unfiltered版本）
  # 提示：确保 dplyr 在最后加载，或统一使用 dplyr:: 前缀（本代码已使用前缀）
  # library(dplyr)  # 如果你不想每次都写 dplyr::，可把这句放到所有包加载的最后
  
  # split N and NN for separate heatmaps
  for (subset_cfg in subsets) {
    subset_tag <- subset_cfg$tag
    message("Generating subset: ", subset_tag)
    
    # ---- 构建主数据（列数据） ----
    go_results_subset <- subset_cfg$filter(go_results) %>%
      { if (!is.null(class_target)) dplyr::filter(., Class %in% class_target) else . } %>%
      dplyr::mutate(
        sample = paste(`Region Subclass`, Sex, sep = "_"),
        type   = ifelse(Neurotransmitter == "NN", "NN", "N")
      ) %>%
      dplyr::select(
        sample = sample,
        term_name = GO,
        value = nlogp,
        region = Region,
        sex = Sex,
        class = Class,
        type
      )
    
    # ---- 构建行分组（行数据） ----
    go_group_subset <- subset_cfg$filter(go_group)
    
    # ---- 预置三套过滤集合 ----
    filtered_sets <- list(
      common_main = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 3, row_col = "GO"), # >10
      common_sup  = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 1, row_col = "GO"),
      specific_sup = NULL
    )
    
    # ---- 新增：unfiltered（不过滤 GO term 或列）----
    filtered_sets <- c(
      list(unfiltered = list(
        main_data = go_results_subset,
        row_data  = go_group_subset
      )),
      filtered_sets
    )
    
    # ---- 基于 common_sup 计算 specific_sup ----
    common_sup_genes <- unique(filtered_sets$common_sup$row_data$GO)
    specific_sup_main <- go_results_subset %>%
      dplyr::filter(!term_name %in% common_sup_genes)
    specific_sup_row <- go_group_subset %>%
      dplyr::filter(!GO %in% common_sup_genes)
    filtered_sets$specific_sup <- list(
      main_data = specific_sup_main,
      row_data  = specific_sup_row
    )
    
    # ---- 聚类配置 ----
    cluster_cfgs <- list(
      list(rows = FALSE, cols = FALSE, tag = "nocluster"),
      list(rows = FALSE, cols = TRUE,  tag = "cluster_col"),
      list(rows = TRUE,  cols = FALSE, tag = "cluster_row"),
      list(rows = TRUE,  cols = TRUE,  tag = "cluster_row_col")
    )
    
    # ---- 逐集合、逐聚类配置画图 ----
    for (set_name in names(filtered_sets)) {
      data_pair <- filtered_sets[[set_name]]
      
      # 兼容列表结构或 data.frame 结构
      main_data <- if (is.list(data_pair) && !is.null(data_pair$main_data)) data_pair$main_data else data_pair
      row_data  <- if (is.list(data_pair) && !is.null(data_pair$row_data))  data_pair$row_data  else go_group_subset
      main_data[main_data$type=="N","class"]=gsub("(.* )|(_.*)","",main_data[main_data$type=="N","sample"])
      
      if (is.null(main_data) || nrow(main_data) == 0) {
        message("Skipping empty dataset: ", set_name)
        next
      }
      
      base_output <- file.path(
        output_dir_path,
        paste0("heatmap_", file_name, "_", subset_tag, "_", set_name)
      )
      
      for (cfg in cluster_cfgs) {
        output_path <- paste0(gsub("_N_","_Glut_",base_output), "_", cfg$tag, ".pdf")
        message("Generating heatmap with config: ", subset_tag, " + ", set_name, " + ", cfg$tag)
        
        results <- heatmap_final2(
          go = main_data,
          row = row_data,
          sample_col = "sample",
          term_col = "term_name",
          value_col = "value",
          value_name = "NES",
          region_col = "region",
          sex_col = if ("sex" %in% names(main_data)) "sex" else NULL,
          class_col = if ("class" %in% names(main_data)) "class" else NULL,
          neurotransmitter_col = NULL,
          row_go_col = "GO",
          row_module_col = if (use_row_category) "module" else NULL,
          row_category_col = NULL,
          col_order_priority = c("region", "sex", "class"),
          col_split_col = "region",
          cutoff = 5,
          distance = "euclidean",
          method = "complete", # "ward.D2" 也可
          logfc_cutoff = NULL,
          column_gap = grid::unit(3, "mm"),
          row_gap = grid::unit(3, "mm"),
          row_names_gp = 6,
          column_names_gp = 6,
          title = output_path,
          cluster_rows = cfg$rows,
          cluster_columns = cfg$cols,
          cluster_column_slices = FALSE,
          cluster_row_slices = FALSE
        )
      }
    }
  }
  
  message("Completed processing for: ", go_path)
}






################other NES#################
######################small_onlyGOBPMF##############
for (go_path in go_paths) {
  message("Processing: ", go_path)
  
  go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1)
  go_results_tmp=go_results_tmp[go_results_tmp$Region!="iCTX",]
  go_results_tmp=go_results_tmp[go_results_tmp$source%in%c("GO:BP","GO:MF"),]
  #go_results_tmp$padj=p.adjust(go_results_tmp$p_value,method = "BH")
  go_results_tmp$nlogp=go_results_tmp$NES
  
  go_results_tmp=go_results_tmp[go_results_tmp$p_value<0.01,]
  colnames(go_results_tmp)=gsub("Region.Subclass","Region Subclass",colnames(go_results_tmp))
  go_results_tmp_unique <- go_results_tmp %>%
    distinct(`Region Subclass`, Sex, GO, .keep_all = TRUE)
  if (nrow(go_results_tmp_unique) < nrow(go_results_tmp)) {
    go_path_unique <- sub("\\.csv$", "_unique.csv", go_path)
    write.csv(go_results_tmp_unique, go_path_unique, row.names = FALSE)
  }
  
  # go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1) %>%
  #   filter(Sex %in% sexes)
  
  #go_results <- switch(target_data,
  #                     "SUS" = go_results_tmp_unique
  #)
  go_results=go_results_tmp_unique
  
  output_dir_path <- dirname(go_path)
  file_name <- tools::file_path_sans_ext(basename(go_path))
  #output_dir_path <- file.path(output_dir_path, file_name)
  output_dir_path <- file.path(output_dir_path, paste0(file_name,"_small_onlyGOBPMF_NES"))
  
  if (!dir.exists(output_dir_path)) {
    dir.create(output_dir_path, recursive = TRUE)
  }
  
  subsets <- list(
    list(tag = "N", filter = function(df) df %>% filter(Neurotransmitter != "NN")),
    list(tag = "NN",    filter = function(df) df %>% filter(Neurotransmitter == "NN"))
  )
  
  group_col <- 'predicted_group' # 'Module.P','predicted_group' or 'final_group'
  use_row_category <-TRUE
  class_target <-  NULL # c("Microglia")
  # go_results$Module.P <- gsub("^G", "M", go_results[[group_col]])
  go_results <- go_results %>%
    mutate(!!group_col := ifelse(.data[[group_col]] == "others", "other", .data[[group_col]]))
  module_levels <- c("G1","G2-1","G2-2","G3-1","G3-2","G3-3","G4","G5","G6","G7","G7-1","G7-2","G8","G8-1","G8-2","G9","G10","G11","G12","other") # c("innate", "overlapped", "adaptive", "inflammatory")
  
  if (use_row_category) {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter,all_of(group_col)) %>%
      rename(module = all_of(group_col)) %>%
      {
        if (!is.null(module_levels)) {
          mutate(., module = factor(module, levels = module_levels))
        } else {
          mutate(., module = factor(module))
        }
      } %>%
      arrange(module)
  } else {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter)
  }
  
  
  # Generate heatmaps for each subset and clustering configuration
  # split N and NN for separate heatmaps
  # 生成不同子集与聚类配置的热图（包含：不过滤GO/列的unfiltered版本）
  # 提示：确保 dplyr 在最后加载，或统一使用 dplyr:: 前缀（本代码已使用前缀）
  # library(dplyr)  # 如果你不想每次都写 dplyr::，可把这句放到所有包加载的最后
  
  # split N and NN for separate heatmaps
  for (subset_cfg in subsets) {
    subset_tag <- subset_cfg$tag
    message("Generating subset: ", subset_tag)
    
    # ---- 构建主数据（列数据） ----
    go_results_subset <- subset_cfg$filter(go_results) %>%
      { if (!is.null(class_target)) dplyr::filter(., Class %in% class_target) else . } %>%
      dplyr::mutate(
        sample = paste(`Region Subclass`, Sex, sep = "_"),
        type   = ifelse(Neurotransmitter == "NN", "NN", "N")
      ) %>%
      dplyr::select(
        sample = sample,
        term_name = GO,
        value = nlogp,
        region = Region,
        sex = Sex,
        class = Class,
        type
      )
    
    # ---- 构建行分组（行数据） ----
    go_group_subset <- subset_cfg$filter(go_group)
    
    # ---- 预置三套过滤集合 ----
    filtered_sets <- list(
      common_main = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 3, row_col = "GO"), # >10
      common_sup  = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 1, row_col = "GO"),
      specific_sup = NULL
    )
    
    # ---- 新增：unfiltered（不过滤 GO term 或列）----
    filtered_sets <- c(
      list(unfiltered = list(
        main_data = go_results_subset,
        row_data  = go_group_subset
      )),
      filtered_sets
    )
    
    # ---- 基于 common_sup 计算 specific_sup ----
    common_sup_genes <- unique(filtered_sets$common_sup$row_data$GO)
    specific_sup_main <- go_results_subset %>%
      dplyr::filter(!term_name %in% common_sup_genes)
    specific_sup_row <- go_group_subset %>%
      dplyr::filter(!GO %in% common_sup_genes)
    filtered_sets$specific_sup <- list(
      main_data = specific_sup_main,
      row_data  = specific_sup_row
    )
    
    # ---- 聚类配置 ----
    cluster_cfgs <- list(
      list(rows = FALSE, cols = FALSE, tag = "nocluster"),
      list(rows = FALSE, cols = TRUE,  tag = "cluster_col"),
      list(rows = TRUE,  cols = FALSE, tag = "cluster_row"),
      list(rows = TRUE,  cols = TRUE,  tag = "cluster_row_col")
    )
    
    # ---- 逐集合、逐聚类配置画图 ----
    for (set_name in names(filtered_sets)) {
      data_pair <- filtered_sets[[set_name]]
      
      # 兼容列表结构或 data.frame 结构
      main_data <- if (is.list(data_pair) && !is.null(data_pair$main_data)) data_pair$main_data else data_pair
      row_data  <- if (is.list(data_pair) && !is.null(data_pair$row_data))  data_pair$row_data  else go_group_subset
      main_data[main_data$type=="N","class"]=gsub("(.* )|(_.*)","",main_data[main_data$type=="N","sample"])
      
      if (is.null(main_data) || nrow(main_data) == 0) {
        message("Skipping empty dataset: ", set_name)
        next
      }
      
      base_output <- file.path(
        output_dir_path,
        paste0("heatmap_", file_name, "_", subset_tag, "_", set_name)
      )
      
      for (cfg in cluster_cfgs) {
        output_path <- paste0(base_output, "_", cfg$tag, ".pdf")
        message("Generating heatmap with config: ", subset_tag, " + ", set_name, " + ", cfg$tag)
        
        results <- heatmap_final3(
          go = main_data,
          row = row_data,
          sample_col = "sample",
          term_col = "term_name",
          value_col = "value",
          value_name = "NES",
          region_col = "region",
          sex_col = if ("sex" %in% names(main_data)) "sex" else NULL,
          class_col = if ("class" %in% names(main_data)) "class" else NULL,
          neurotransmitter_col = NULL,
          row_go_col = "GO",
          row_module_col = if (use_row_category) "module" else NULL,
          row_category_col = NULL,
          col_order_priority = c("region", "sex", "class"),
          col_split_col = "region",
          cutoff = 5,
          distance = "euclidean",
          method = "complete", # "ward.D2" 也可
          logfc_cutoff = NULL,
          column_gap = grid::unit(3, "mm"),
          row_gap = grid::unit(3, "mm"),
          row_names_gp = 6,
          column_names_gp = 6,
          title = output_path,
          cluster_rows = cfg$rows,
          cluster_columns = cfg$cols,
          cluster_column_slices = FALSE,
          cluster_row_slices = FALSE
        )
      }
    }
  }
  
  message("Completed processing for: ", go_path)
}
#########################onlyGo BPMF###########
for (go_path in go_paths) {
  message("Processing: ", go_path)
  
  go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1)
  go_results_tmp=go_results_tmp[go_results_tmp$Region!="iCTX",]
  go_results_tmp=go_results_tmp[go_results_tmp$source%in%c("GO:BP","GO:MF"),]
  #go_results_tmp$padj=p.adjust(go_results_tmp$p_value,method = "BH")
  go_results_tmp$nlogp=go_results_tmp$NES
  
  go_results_tmp=go_results_tmp[go_results_tmp$p_value<0.01,]
  colnames(go_results_tmp)=gsub("Region.Subclass","Region Subclass",colnames(go_results_tmp))
  go_results_tmp_unique <- go_results_tmp %>%
    distinct(`Region Subclass`, Sex, GO, .keep_all = TRUE)
  if (nrow(go_results_tmp_unique) < nrow(go_results_tmp)) {
    go_path_unique <- sub("\\.csv$", "_unique.csv", go_path)
    write.csv(go_results_tmp_unique, go_path_unique, row.names = FALSE)
  }
  
  # go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1) %>%
  #   filter(Sex %in% sexes)
  
  #go_results <- switch(target_data,
  #                     "SUS" = go_results_tmp_unique
  #)
  go_results=go_results_tmp_unique
  
  output_dir_path <- dirname(go_path)
  file_name <- tools::file_path_sans_ext(basename(go_path))
  #output_dir_path <- file.path(output_dir_path, file_name)
  output_dir_path <- file.path(output_dir_path, paste0(file_name,"_onlyGOBMF_NES"))
  
  if (!dir.exists(output_dir_path)) {
    dir.create(output_dir_path, recursive = TRUE)
  }
  
  subsets <- list(
    list(tag = "N", filter = function(df) df %>% filter(Neurotransmitter != "NN")),
    list(tag = "NN",    filter = function(df) df %>% filter(Neurotransmitter == "NN"))
  )
  
  group_col <- 'predicted_group' # 'Module.P','predicted_group' or 'final_group'
  use_row_category <-TRUE
  class_target <-  NULL # c("Microglia")
  # go_results$Module.P <- gsub("^G", "M", go_results[[group_col]])
  go_results <- go_results %>%
    mutate(!!group_col := ifelse(.data[[group_col]] == "others", "other", .data[[group_col]]))
  module_levels <- c("G1","G2-1","G2-2","G3-1","G3-2","G3-3","G4","G5","G6","G7","G7-1","G7-2","G8","G8-1","G8-2","G9","G10","G11","G12","other") # c("innate", "overlapped", "adaptive", "inflammatory")
  
  if (use_row_category) {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter,all_of(group_col)) %>%
      rename(module = all_of(group_col)) %>%
      {
        if (!is.null(module_levels)) {
          mutate(., module = factor(module, levels = module_levels))
        } else {
          mutate(., module = factor(module))
        }
      } %>%
      arrange(module)
  } else {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter)
  }
  
  
  # Generate heatmaps for each subset and clustering configuration
  # split N and NN for separate heatmaps
  # 生成不同子集与聚类配置的热图（包含：不过滤GO/列的unfiltered版本）
  # 提示：确保 dplyr 在最后加载，或统一使用 dplyr:: 前缀（本代码已使用前缀）
  # library(dplyr)  # 如果你不想每次都写 dplyr::，可把这句放到所有包加载的最后
  
  # split N and NN for separate heatmaps
  for (subset_cfg in subsets) {
    subset_tag <- subset_cfg$tag
    message("Generating subset: ", subset_tag)
    
    # ---- 构建主数据（列数据） ----
    go_results_subset <- subset_cfg$filter(go_results) %>%
      { if (!is.null(class_target)) dplyr::filter(., Class %in% class_target) else . } %>%
      dplyr::mutate(
        sample = paste(`Region Subclass`, Sex, sep = "_"),
        type   = ifelse(Neurotransmitter == "NN", "NN", "N")
      ) %>%
      dplyr::select(
        sample = sample,
        term_name = GO,
        value = nlogp,
        region = Region,
        sex = Sex,
        class = Class,
        type
      )
    
    # ---- 构建行分组（行数据） ----
    go_group_subset <- subset_cfg$filter(go_group)
    
    # ---- 预置三套过滤集合 ----
    filtered_sets <- list(
      common_main = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 3, row_col = "GO"), # >10
      common_sup  = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 1, row_col = "GO"),
      specific_sup = NULL
    )
    
    # ---- 新增：unfiltered（不过滤 GO term 或列）----
    filtered_sets <- c(
      list(unfiltered = list(
        main_data = go_results_subset,
        row_data  = go_group_subset
      )),
      filtered_sets
    )
    
    # ---- 基于 common_sup 计算 specific_sup ----
    common_sup_genes <- unique(filtered_sets$common_sup$row_data$GO)
    specific_sup_main <- go_results_subset %>%
      dplyr::filter(!term_name %in% common_sup_genes)
    specific_sup_row <- go_group_subset %>%
      dplyr::filter(!GO %in% common_sup_genes)
    filtered_sets$specific_sup <- list(
      main_data = specific_sup_main,
      row_data  = specific_sup_row
    )
    
    # ---- 聚类配置 ----
    cluster_cfgs <- list(
      list(rows = FALSE, cols = FALSE, tag = "nocluster"),
      list(rows = FALSE, cols = TRUE,  tag = "cluster_col"),
      list(rows = TRUE,  cols = FALSE, tag = "cluster_row"),
      list(rows = TRUE,  cols = TRUE,  tag = "cluster_row_col")
    )
    
    # ---- 逐集合、逐聚类配置画图 ----
    for (set_name in names(filtered_sets)) {
      data_pair <- filtered_sets[[set_name]]
      
      # 兼容列表结构或 data.frame 结构
      main_data <- if (is.list(data_pair) && !is.null(data_pair$main_data)) data_pair$main_data else data_pair
      row_data  <- if (is.list(data_pair) && !is.null(data_pair$row_data))  data_pair$row_data  else go_group_subset
      main_data[main_data$type=="N","class"]=gsub("(.* )|(_.*)","",main_data[main_data$type=="N","sample"])
      
      if (is.null(main_data) || nrow(main_data) == 0) {
        message("Skipping empty dataset: ", set_name)
        next
      }
      
      base_output <- file.path(
        output_dir_path,
        paste0("heatmap_", file_name, "_", subset_tag, "_", set_name)
      )
      
      for (cfg in cluster_cfgs) {
        output_path <- paste0(base_output, "_", cfg$tag, ".pdf")
        message("Generating heatmap with config: ", subset_tag, " + ", set_name, " + ", cfg$tag)
        
        results <- heatmap_final2(
          go = main_data,
          row = row_data,
          sample_col = "sample",
          term_col = "term_name",
          value_col = "value",
          value_name = "NES",
          region_col = "region",
          sex_col = if ("sex" %in% names(main_data)) "sex" else NULL,
          class_col = if ("class" %in% names(main_data)) "class" else NULL,
          neurotransmitter_col = NULL,
          row_go_col = "GO",
          row_module_col = if (use_row_category) "module" else NULL,
          row_category_col = NULL,
          col_order_priority = c("region", "sex", "class"),
          col_split_col = "region",
          cutoff = 5,
          distance = "euclidean",
          method = "complete", # "ward.D2" 也可
          logfc_cutoff = NULL,
          column_gap = grid::unit(3, "mm"),
          row_gap = grid::unit(3, "mm"),
          row_names_gp = 6,
          column_names_gp = 6,
          title = output_path,
          cluster_rows = cfg$rows,
          cluster_columns = cfg$cols,
          cluster_column_slices = FALSE,
          cluster_row_slices = FALSE
        )
      }
    }
  }
  
  message("Completed processing for: ", go_path)
}
#########################alldatabase###########
for (go_path in go_paths) {
  message("Processing: ", go_path)
  
  go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1)
  go_results_tmp=go_results_tmp[go_results_tmp$Region!="iCTX",]
  go_results_tmp=go_results_tmp[go_results_tmp$p_value<0.01,]
  go_results_tmp$nlogp=go_results_tmp$NES
  
  #go_results_tmp=go_results_tmp[go_results_tmp$source%in%c("GO:BP","GO:MF"),]
  colnames(go_results_tmp)=gsub("Region.Subclass","Region Subclass",colnames(go_results_tmp))
  go_results_tmp_unique <- go_results_tmp %>%
    distinct(`Region Subclass`, Sex, GO, .keep_all = TRUE)
  if (nrow(go_results_tmp_unique) < nrow(go_results_tmp)) {
    go_path_unique <- sub("\\.csv$", "_unique.csv", go_path)
    write.csv(go_results_tmp_unique, go_path_unique, row.names = FALSE)
  }
  
  # go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1) %>%
  #   filter(Sex %in% sexes)
  
  #go_results <- switch(target_data,
  #                     "SUS" = go_results_tmp_unique
  #)
  go_results=go_results_tmp_unique
  
  output_dir_path <- dirname(go_path)
  file_name <- tools::file_path_sans_ext(basename(go_path))
  #output_dir_path <- file.path(output_dir_path, file_name)
  output_dir_path <- file.path(output_dir_path, paste0(file_name,"_alldatabase_NES"))
  
  if (!dir.exists(output_dir_path)) {
    dir.create(output_dir_path, recursive = TRUE)
  }
  
  subsets <- list(
    list(tag = "N", filter = function(df) df %>% filter(Neurotransmitter != "NN")),
    list(tag = "NN",    filter = function(df) df %>% filter(Neurotransmitter == "NN"))
  )
  
  group_col <- 'predicted_group' # 'Module.P','predicted_group' or 'final_group'
  use_row_category <-TRUE
  class_target <-  NULL # c("Microglia")
  # go_results$Module.P <- gsub("^G", "M", go_results[[group_col]])
  go_results <- go_results %>%
    mutate(!!group_col := ifelse(.data[[group_col]] == "others", "other", .data[[group_col]]))
  module_levels <- c("G1","G2-1","G2-2","G3-1","G3-2","G3-3","G4","G5","G6","G7","G7-1","G7-2","G8","G8-1","G8-2","G9","G10","G11","G12","other") # c("innate", "overlapped", "adaptive", "inflammatory")
  
  if (use_row_category) {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter,all_of(group_col)) %>%
      rename(module = all_of(group_col)) %>%
      {
        if (!is.null(module_levels)) {
          mutate(., module = factor(module, levels = module_levels))
        } else {
          mutate(., module = factor(module))
        }
      } %>%
      arrange(module)
  } else {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter)
  }
  
  
  # Generate heatmaps for each subset and clustering configuration
  # split N and NN for separate heatmaps
  # 生成不同子集与聚类配置的热图（包含：不过滤GO/列的unfiltered版本）
  # 提示：确保 dplyr 在最后加载，或统一使用 dplyr:: 前缀（本代码已使用前缀）
  # library(dplyr)  # 如果你不想每次都写 dplyr::，可把这句放到所有包加载的最后
  
  # split N and NN for separate heatmaps
  for (subset_cfg in subsets) {
    subset_tag <- subset_cfg$tag
    message("Generating subset: ", subset_tag)
    
    # ---- 构建主数据（列数据） ----
    go_results_subset <- subset_cfg$filter(go_results) %>%
      { if (!is.null(class_target)) dplyr::filter(., Class %in% class_target) else . } %>%
      dplyr::mutate(
        sample = paste(`Region Subclass`, Sex, sep = "_"),
        type   = ifelse(Neurotransmitter == "NN", "NN", "N")
      ) %>%
      dplyr::select(
        sample = sample,
        term_name = GO,
        value = nlogp,
        region = Region,
        sex = Sex,
        class = Class,
        type
      )
    
    # ---- 构建行分组（行数据） ----
    go_group_subset <- subset_cfg$filter(go_group)
    
    # ---- 预置三套过滤集合 ----
    filtered_sets <- list(
      common_main = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 3, row_col = "GO"), # >10
      common_sup  = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 1, row_col = "GO"),
      specific_sup = NULL
    )
    
    # ---- 新增：unfiltered（不过滤 GO term 或列）----
    filtered_sets <- c(
      list(unfiltered = list(
        main_data = go_results_subset,
        row_data  = go_group_subset
      )),
      filtered_sets
    )
    
    # ---- 基于 common_sup 计算 specific_sup ----
    common_sup_genes <- unique(filtered_sets$common_sup$row_data$GO)
    specific_sup_main <- go_results_subset %>%
      dplyr::filter(!term_name %in% common_sup_genes)
    specific_sup_row <- go_group_subset %>%
      dplyr::filter(!GO %in% common_sup_genes)
    filtered_sets$specific_sup <- list(
      main_data = specific_sup_main,
      row_data  = specific_sup_row
    )
    
    # ---- 聚类配置 ----
    cluster_cfgs <- list(
      list(rows = FALSE, cols = FALSE, tag = "nocluster"),
      list(rows = FALSE, cols = TRUE,  tag = "cluster_col"),
      list(rows = TRUE,  cols = FALSE, tag = "cluster_row"),
      list(rows = TRUE,  cols = TRUE,  tag = "cluster_row_col")
    )
    
    # ---- 逐集合、逐聚类配置画图 ----
    for (set_name in names(filtered_sets)) {
      data_pair <- filtered_sets[[set_name]]
      
      # 兼容列表结构或 data.frame 结构
      main_data <- if (is.list(data_pair) && !is.null(data_pair$main_data)) data_pair$main_data else data_pair
      row_data  <- if (is.list(data_pair) && !is.null(data_pair$row_data))  data_pair$row_data  else go_group_subset
      main_data[main_data$type=="N","class"]=gsub("(.* )|(_.*)","",main_data[main_data$type=="N","sample"])
      
      if (is.null(main_data) || nrow(main_data) == 0) {
        message("Skipping empty dataset: ", set_name)
        next
      }
      
      base_output <- file.path(
        output_dir_path,
        paste0("heatmap_", file_name, "_", subset_tag, "_", set_name)
      )
      
      for (cfg in cluster_cfgs) {
        output_path <- paste0(base_output, "_", cfg$tag, ".pdf")
        message("Generating heatmap with config: ", subset_tag, " + ", set_name, " + ", cfg$tag)
        
        results <- heatmap_final2(
          go = main_data,
          row = row_data,
          sample_col = "sample",
          term_col = "term_name",
          value_col = "value",
          value_name = "NES",
          region_col = "region",
          sex_col = if ("sex" %in% names(main_data)) "sex" else NULL,
          class_col = if ("class" %in% names(main_data)) "class" else NULL,
          neurotransmitter_col = NULL,
          row_go_col = "GO",
          row_module_col = if (use_row_category) "module" else NULL,
          row_category_col = NULL,
          col_order_priority = c("region", "sex", "class"),
          col_split_col = "region",
          cutoff = 5,
          distance = "euclidean",
          method = "complete", # "ward.D2" 也可
          logfc_cutoff = NULL,
          column_gap = grid::unit(3, "mm"),
          row_gap = grid::unit(3, "mm"),
          row_names_gp = 6,
          column_names_gp = 6,
          title = output_path,
          cluster_rows = cfg$rows,
          cluster_columns = cfg$cols,
          cluster_column_slices = FALSE,
          cluster_row_slices = FALSE
        )
      }
    }
  }
  
  message("Completed processing for: ", go_path)
}

##############################################alldatabase_small
for (go_path in go_paths) {
  message("Processing: ", go_path)
  
  go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1)
  go_results_tmp[go_results_tmp$Region!="iCTX",]
  colnames(go_results_tmp)=gsub("Region.Subclass","Region Subclass",colnames(go_results_tmp))
  go_results_tmp$nlogp=go_results_tmp$NES
  go_results_tmp_unique <- go_results_tmp %>%
    distinct(`Region Subclass`, Sex, GO, .keep_all = TRUE)
  if (nrow(go_results_tmp_unique) < nrow(go_results_tmp)) {
    go_path_unique <- sub("\\.csv$", "_unique.csv", go_path)
    write.csv(go_results_tmp_unique, go_path_unique, row.names = FALSE)
  }
  
  # go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1) %>%
  #   filter(Sex %in% sexes)
  
  #go_results <- switch(target_data,
  #                     "SUS" = go_results_tmp_unique
  #)
  go_results=go_results_tmp_unique
  
  output_dir_path <- dirname(go_path)
  file_name <- tools::file_path_sans_ext(basename(go_path))
  output_dir_path <- file.path(output_dir_path, file_name)
  output_dir_path1 <- file.path(output_dir_path, paste0(file_name,"_alldatabasesmall"))
  
  if (!dir.exists(output_dir_path)) {
    dir.create(output_dir_path, recursive = TRUE)
  }
  
  subsets <- list(
    list(tag = "N", filter = function(df) df %>% filter(Neurotransmitter != "NN")),
    list(tag = "NN",    filter = function(df) df %>% filter(Neurotransmitter == "NN"))
  )
  
  group_col <- 'predicted_group' # 'Module.P','predicted_group' or 'final_group'
  use_row_category <-TRUE
  class_target <-  NULL # c("Microglia")
  # go_results$Module.P <- gsub("^G", "M", go_results[[group_col]])
  go_results <- go_results %>%
    mutate(!!group_col := ifelse(.data[[group_col]] == "others", "other", .data[[group_col]]))
  module_levels <- c("G1","G2-1","G2-2","G3-1","G3-2","G3-3","G4","G5","G6","G7","G7-1","G7-2","G8","G8-1","G8-2","G9","G10","G11","G12","other") # c("innate", "overlapped", "adaptive", "inflammatory")
  
  if (use_row_category) {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter,all_of(group_col)) %>%
      rename(module = all_of(group_col)) %>%
      {
        if (!is.null(module_levels)) {
          mutate(., module = factor(module, levels = module_levels))
        } else {
          mutate(., module = factor(module))
        }
      } %>%
      arrange(module)
  } else {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter)
  }
  
  
  # Generate heatmaps for each subset and clustering configuration
  # split N and NN for separate heatmaps
  # 生成不同子集与聚类配置的热图（包含：不过滤GO/列的unfiltered版本）
  # 提示：确保 dplyr 在最后加载，或统一使用 dplyr:: 前缀（本代码已使用前缀）
  # library(dplyr)  # 如果你不想每次都写 dplyr::，可把这句放到所有包加载的最后
  
  # split N and NN for separate heatmaps
  for (subset_cfg in subsets) {
    subset_tag <- subset_cfg$tag
    message("Generating subset: ", subset_tag)
    
    # ---- 构建主数据（列数据） ----
    go_results_subset <- subset_cfg$filter(go_results) %>%
      { if (!is.null(class_target)) dplyr::filter(., Class %in% class_target) else . } %>%
      dplyr::mutate(
        sample = paste(`Region Subclass`, Sex, sep = "_"),
        type   = ifelse(Neurotransmitter == "NN", "NN", "N")
      ) %>%
      dplyr::select(
        sample = sample,
        term_name = GO,
        value = nlogp,
        region = Region,
        sex = Sex,
        class = Class,
        type
      )
    
    # ---- 构建行分组（行数据） ----
    go_group_subset <- subset_cfg$filter(go_group)
    
    # ---- 预置三套过滤集合 ----
    filtered_sets <- list(
      common_main = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 3, row_col = "GO"), # >10
      common_sup  = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 1, row_col = "GO"),
      specific_sup = NULL
    )
    
    # ---- 新增：unfiltered（不过滤 GO term 或列）----
    filtered_sets <- c(
      list(unfiltered = list(
        main_data = go_results_subset,
        row_data  = go_group_subset
      )),
      filtered_sets
    )
    
    # ---- 基于 common_sup 计算 specific_sup ----
    common_sup_genes <- unique(filtered_sets$common_sup$row_data$GO)
    specific_sup_main <- go_results_subset %>%
      dplyr::filter(!term_name %in% common_sup_genes)
    specific_sup_row <- go_group_subset %>%
      dplyr::filter(!GO %in% common_sup_genes)
    filtered_sets$specific_sup <- list(
      main_data = specific_sup_main,
      row_data  = specific_sup_row
    )
    
    # ---- 聚类配置 ----
    cluster_cfgs <- list(
      list(rows = FALSE, cols = FALSE, tag = "nocluster"),
      list(rows = FALSE, cols = TRUE,  tag = "cluster_col"),
      list(rows = TRUE,  cols = FALSE, tag = "cluster_row"),
      list(rows = TRUE,  cols = TRUE,  tag = "cluster_row_col")
    )
    
    # ---- 逐集合、逐聚类配置画图 ----
    for (set_name in names(filtered_sets)) {
      data_pair <- filtered_sets[[set_name]]
      
      # 兼容列表结构或 data.frame 结构
      main_data <- if (is.list(data_pair) && !is.null(data_pair$main_data)) data_pair$main_data else data_pair
      row_data  <- if (is.list(data_pair) && !is.null(data_pair$row_data))  data_pair$row_data  else go_group_subset
      main_data[main_data$type=="N","class"]=gsub("(.* )|(_.*)","",main_data[main_data$type=="N","sample"])
      
      if (is.null(main_data) || nrow(main_data) == 0) {
        message("Skipping empty dataset: ", set_name)
        next
      }
      
      base_output <- file.path(
        output_dir_path,
        paste0("heatmap_", file_name, "_", subset_tag, "_", set_name)
      )
      
      for (cfg in cluster_cfgs) {
        output_path <- paste0(base_output, "_", cfg$tag, ".pdf")
        message("Generating heatmap with config: ", subset_tag, " + ", set_name, " + ", cfg$tag)
        
        results <- heatmap_final2(
          go = main_data,
          row = row_data,
          sample_col = "sample",
          term_col = "term_name",
          value_col = "value",
          value_name = "NES",
          region_col = "region",
          sex_col = if ("sex" %in% names(main_data)) "sex" else NULL,
          class_col = if ("class" %in% names(main_data)) "class" else NULL,
          neurotransmitter_col = NULL,
          row_go_col = "GO",
          row_module_col = if (use_row_category) "module" else NULL,
          row_category_col = NULL,
          col_order_priority = c("region", "sex", "class"),
          col_split_col = "region",
          cutoff = 5,
          distance = "euclidean",
          method = "complete", # "ward.D2" 也可
          logfc_cutoff = NULL,
          column_gap = grid::unit(3, "mm"),
          row_gap = grid::unit(3, "mm"),
          row_names_gp = 6,
          column_names_gp = 6,
          title = output_path,
          cluster_rows = cfg$rows,
          cluster_columns = cfg$cols,
          cluster_column_slices = FALSE,
          cluster_row_slices = FALSE
        )
      }
    }
  }
  
  message("Completed processing for: ", go_path)
}





######################small_onlyGO BPMF ##############
for (go_path in go_paths) {
  message("Processing: ", go_path)
  
  go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1)
  go_results_tmp=go_results_tmp[go_results_tmp$Region!="iCTX",]
  go_results_tmp=go_results_tmp[go_results_tmp$source%in%c("GO:BP","GO:MF"),]
  #go_results_tmp$padj=p.adjust(go_results_tmp$p_value,method = "BH")
  go_results_tmp$nlogp=-log10(go_results_tmp$padj)*sign(go_results_tmp$nlogp)
  go_results_tmp$nlog10p=-log10(go_results_tmp$padj)
  
  #go_results_tmp=go_results_tmp[go_results_tmp$p_value<0.01,]
  colnames(go_results_tmp)=gsub("Region.Subclass","Region Subclass",colnames(go_results_tmp))
  go_results_tmp_unique <- go_results_tmp %>%
    distinct(`Region Subclass`, Sex, GO, .keep_all = TRUE)
  if (nrow(go_results_tmp_unique) < nrow(go_results_tmp)) {
    go_path_unique <- sub("\\.csv$", "_unique.csv", go_path)
    write.csv(go_results_tmp_unique, go_path_unique, row.names = FALSE)
  }
  
  # go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1) %>%
  #   filter(Sex %in% sexes)
  
  #go_results <- switch(target_data,
  #                     "SUS" = go_results_tmp_unique
  #)
  go_results=go_results_tmp_unique
  
  output_dir_path <- dirname(go_path)
  file_name <- tools::file_path_sans_ext(basename(go_path))
  #output_dir_path <- file.path(output_dir_path, file_name)
  output_dir_path <- file.path(output_dir_path, paste0(file_name,"_small_onlyGOBPMF"))
  
  if (!dir.exists(output_dir_path)) {
    dir.create(output_dir_path, recursive = TRUE)
  }
  
  subsets <- list(
    list(tag = "N", filter = function(df) df %>% filter(Neurotransmitter != "NN")),
    list(tag = "NN",    filter = function(df) df %>% filter(Neurotransmitter == "NN"))
  )
  
  group_col <- 'predicted_group' # 'Module.P','predicted_group' or 'final_group'
  use_row_category <-TRUE
  class_target <-  NULL # c("Microglia")
  # go_results$Module.P <- gsub("^G", "M", go_results[[group_col]])
  go_results <- go_results %>%
    mutate(!!group_col := ifelse(.data[[group_col]] == "others", "other", .data[[group_col]]))
  module_levels <- c("G1","G2-1","G2-2","G3-1","G3-2","G3-3","G4","G5","G6","G7","G7-1","G7-2","G8","G8-1","G8-2","G9","G10","G11","G12","other") # c("innate", "overlapped", "adaptive", "inflammatory")
  
  if (use_row_category) {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter,all_of(group_col)) %>%
      rename(module = all_of(group_col)) %>%
      {
        if (!is.null(module_levels)) {
          mutate(., module = factor(module, levels = module_levels))
        } else {
          mutate(., module = factor(module))
        }
      } %>%
      arrange(module)
  } else {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter)
  }
  
  
  # Generate heatmaps for each subset and clustering configuration
  # split N and NN for separate heatmaps
  # 生成不同子集与聚类配置的热图（包含：不过滤GO/列的unfiltered版本）
  # 提示：确保 dplyr 在最后加载，或统一使用 dplyr:: 前缀（本代码已使用前缀）
  # library(dplyr)  # 如果你不想每次都写 dplyr::，可把这句放到所有包加载的最后
  
  # split N and NN for separate heatmaps
  for (subset_cfg in subsets) {
    subset_tag <- subset_cfg$tag
    message("Generating subset: ", subset_tag)
    
    # ---- 构建主数据（列数据） ----
    go_results_subset <- subset_cfg$filter(go_results) %>%
      { if (!is.null(class_target)) dplyr::filter(., Class %in% class_target) else . } %>%
      dplyr::mutate(
        sample = paste(`Region Subclass`, Sex, sep = "_"),
        type   = ifelse(Neurotransmitter == "NN", "NN", "N")
      ) %>%
      dplyr::select(
        sample = sample,
        term_name = GO,
        value = nlogp,
        region = Region,
        sex = Sex,
        class = Class,
        type
      )
    
    # ---- 构建行分组（行数据） ----
    go_group_subset <- subset_cfg$filter(go_group)
    
    # ---- 预置三套过滤集合 ----
    filtered_sets <- list(
      common_main = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 3, row_col = "GO"), # >10
      common_sup  = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 1, row_col = "GO"),
      specific_sup = NULL
    )
    
    # ---- 新增：unfiltered（不过滤 GO term 或列）----
    filtered_sets <- c(
      list(unfiltered = list(
        main_data = go_results_subset,
        row_data  = go_group_subset
      )),
      filtered_sets
    )
    
    # ---- 基于 common_sup 计算 specific_sup ----
    common_sup_genes <- unique(filtered_sets$common_sup$row_data$GO)
    specific_sup_main <- go_results_subset %>%
      dplyr::filter(!term_name %in% common_sup_genes)
    specific_sup_row <- go_group_subset %>%
      dplyr::filter(!GO %in% common_sup_genes)
    filtered_sets$specific_sup <- list(
      main_data = specific_sup_main,
      row_data  = specific_sup_row
    )
    
    # ---- 聚类配置 ----
    cluster_cfgs <- list(
      list(rows = FALSE, cols = FALSE, tag = "nocluster"),
      list(rows = FALSE, cols = TRUE,  tag = "cluster_col"),
      list(rows = TRUE,  cols = FALSE, tag = "cluster_row"),
      list(rows = TRUE,  cols = TRUE,  tag = "cluster_row_col")
    )
    
    # ---- 逐集合、逐聚类配置画图 ----
    for (set_name in names(filtered_sets)) {
      data_pair <- filtered_sets[[set_name]]
      
      # 兼容列表结构或 data.frame 结构
      main_data <- if (is.list(data_pair) && !is.null(data_pair$main_data)) data_pair$main_data else data_pair
      row_data  <- if (is.list(data_pair) && !is.null(data_pair$row_data))  data_pair$row_data  else go_group_subset
      main_data[main_data$type=="N","class"]=gsub("(.* )|(_.*)","",main_data[main_data$type=="N","sample"])
      
      if (is.null(main_data) || nrow(main_data) == 0) {
        message("Skipping empty dataset: ", set_name)
        next
      }
      
      base_output <- file.path(
        output_dir_path,
        paste0("heatmap_", file_name, "_", subset_tag, "_", set_name)
      )
      
      for (cfg in cluster_cfgs) {
        output_path <- paste0(base_output, "_", cfg$tag, ".pdf")
        message("Generating heatmap with config: ", subset_tag, " + ", set_name, " + ", cfg$tag)
        
        results <- heatmap_final3(
          go = main_data,
          row = row_data,
          sample_col = "sample",
          term_col = "term_name",
          value_col = "value",
          value_name = "-log10(p_val_adj)",
          region_col = "region",
          sex_col = if ("sex" %in% names(main_data)) "sex" else NULL,
          class_col = if ("class" %in% names(main_data)) "class" else NULL,
          neurotransmitter_col = NULL,
          row_go_col = "GO",
          row_module_col = if (use_row_category) "module" else NULL,
          row_category_col = NULL,
          col_order_priority = c("region", "sex", "class"),
          col_split_col = "region",
          cutoff = 10,
          distance = "euclidean",
          method = "complete", # "ward.D2" 也可
          logfc_cutoff = NULL,
          column_gap = grid::unit(3, "mm"),
          row_gap = grid::unit(3, "mm"),
          row_names_gp = 6,
          column_names_gp = 6,
          title = output_path,
          cluster_rows = cfg$rows,
          cluster_columns = cfg$cols,
          cluster_column_slices = FALSE,
          cluster_row_slices = FALSE
        )
      }
    }
  }
  
  message("Completed processing for: ", go_path)
}
#########################onlyGo BPMF###########
for (go_path in go_paths) {
  message("Processing: ", go_path)
  
  go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1)
  go_results_tmp=go_results_tmp[go_results_tmp$Region!="iCTX",]
  go_results_tmp=go_results_tmp[go_results_tmp$source%in%c("GO:BP","GO:MF"),]
  #go_results_tmp$padj=p.adjust(go_results_tmp$p_value,method = "BH")
  go_results_tmp$nlogp=-log10(go_results_tmp$padj)*sign(go_results_tmp$nlogp)
  go_results_tmp$nlog10p=-log10(go_results_tmp$padj)
  
  #go_results_tmp=go_results_tmp[go_results_tmp$p_value<0.01,]
  colnames(go_results_tmp)=gsub("Region.Subclass","Region Subclass",colnames(go_results_tmp))
  go_results_tmp_unique <- go_results_tmp %>%
    distinct(`Region Subclass`, Sex, GO, .keep_all = TRUE)
  if (nrow(go_results_tmp_unique) < nrow(go_results_tmp)) {
    go_path_unique <- sub("\\.csv$", "_unique.csv", go_path)
    write.csv(go_results_tmp_unique, go_path_unique, row.names = FALSE)
  }
  
  # go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1) %>%
  #   filter(Sex %in% sexes)
  
  #go_results <- switch(target_data,
  #                     "SUS" = go_results_tmp_unique
  #)
  go_results=go_results_tmp_unique
  
  output_dir_path <- dirname(go_path)
  file_name <- tools::file_path_sans_ext(basename(go_path))
  #output_dir_path <- file.path(output_dir_path, file_name)
  output_dir_path <- file.path(output_dir_path, paste0(file_name,"_onlyGOMF"))
  
  if (!dir.exists(output_dir_path)) {
    dir.create(output_dir_path, recursive = TRUE)
  }
  
  subsets <- list(
    list(tag = "N", filter = function(df) df %>% filter(Neurotransmitter != "NN")),
    list(tag = "NN",    filter = function(df) df %>% filter(Neurotransmitter == "NN"))
  )
  
  group_col <- 'predicted_group' # 'Module.P','predicted_group' or 'final_group'
  use_row_category <-TRUE
  class_target <-  NULL # c("Microglia")
  # go_results$Module.P <- gsub("^G", "M", go_results[[group_col]])
  go_results <- go_results %>%
    mutate(!!group_col := ifelse(.data[[group_col]] == "others", "other", .data[[group_col]]))
  module_levels <- c("G1","G2-1","G2-2","G3-1","G3-2","G3-3","G4","G5","G6","G7","G7-1","G7-2","G8","G8-1","G8-2","G9","G10","G11","G12","other") # c("innate", "overlapped", "adaptive", "inflammatory")
  
  if (use_row_category) {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter,all_of(group_col)) %>%
      rename(module = all_of(group_col)) %>%
      {
        if (!is.null(module_levels)) {
          mutate(., module = factor(module, levels = module_levels))
        } else {
          mutate(., module = factor(module))
        }
      } %>%
      arrange(module)
  } else {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter)
  }
  
  
  # Generate heatmaps for each subset and clustering configuration
  # split N and NN for separate heatmaps
  # 生成不同子集与聚类配置的热图（包含：不过滤GO/列的unfiltered版本）
  # 提示：确保 dplyr 在最后加载，或统一使用 dplyr:: 前缀（本代码已使用前缀）
  # library(dplyr)  # 如果你不想每次都写 dplyr::，可把这句放到所有包加载的最后
  
  # split N and NN for separate heatmaps
  for (subset_cfg in subsets) {
    subset_tag <- subset_cfg$tag
    message("Generating subset: ", subset_tag)
    
    # ---- 构建主数据（列数据） ----
    go_results_subset <- subset_cfg$filter(go_results) %>%
      { if (!is.null(class_target)) dplyr::filter(., Class %in% class_target) else . } %>%
      dplyr::mutate(
        sample = paste(`Region Subclass`, Sex, sep = "_"),
        type   = ifelse(Neurotransmitter == "NN", "NN", "N")
      ) %>%
      dplyr::select(
        sample = sample,
        term_name = GO,
        value = nlogp,
        region = Region,
        sex = Sex,
        class = Class,
        type
      )
    
    # ---- 构建行分组（行数据） ----
    go_group_subset <- subset_cfg$filter(go_group)
    
    # ---- 预置三套过滤集合 ----
    filtered_sets <- list(
      common_main = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 3, row_col = "GO"), # >10
      common_sup  = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 1, row_col = "GO"),
      specific_sup = NULL
    )
    
    # ---- 新增：unfiltered（不过滤 GO term 或列）----
    filtered_sets <- c(
      list(unfiltered = list(
        main_data = go_results_subset,
        row_data  = go_group_subset
      )),
      filtered_sets
    )
    
    # ---- 基于 common_sup 计算 specific_sup ----
    common_sup_genes <- unique(filtered_sets$common_sup$row_data$GO)
    specific_sup_main <- go_results_subset %>%
      dplyr::filter(!term_name %in% common_sup_genes)
    specific_sup_row <- go_group_subset %>%
      dplyr::filter(!GO %in% common_sup_genes)
    filtered_sets$specific_sup <- list(
      main_data = specific_sup_main,
      row_data  = specific_sup_row
    )
    
    # ---- 聚类配置 ----
    cluster_cfgs <- list(
      list(rows = FALSE, cols = FALSE, tag = "nocluster"),
      list(rows = FALSE, cols = TRUE,  tag = "cluster_col"),
      list(rows = TRUE,  cols = FALSE, tag = "cluster_row"),
      list(rows = TRUE,  cols = TRUE,  tag = "cluster_row_col")
    )
    
    # ---- 逐集合、逐聚类配置画图 ----
    for (set_name in names(filtered_sets)) {
      data_pair <- filtered_sets[[set_name]]
      
      # 兼容列表结构或 data.frame 结构
      main_data <- if (is.list(data_pair) && !is.null(data_pair$main_data)) data_pair$main_data else data_pair
      row_data  <- if (is.list(data_pair) && !is.null(data_pair$row_data))  data_pair$row_data  else go_group_subset
      main_data[main_data$type=="N","class"]=gsub("(.* )|(_.*)","",main_data[main_data$type=="N","sample"])
      
      if (is.null(main_data) || nrow(main_data) == 0) {
        message("Skipping empty dataset: ", set_name)
        next
      }
      
      base_output <- file.path(
        output_dir_path,
        paste0("heatmap_", file_name, "_", subset_tag, "_", set_name)
      )
      
      for (cfg in cluster_cfgs) {
        output_path <- paste0(base_output, "_", cfg$tag, ".pdf")
        message("Generating heatmap with config: ", subset_tag, " + ", set_name, " + ", cfg$tag)
        
        results <- heatmap_final2(
          go = main_data,
          row = row_data,
          sample_col = "sample",
          term_col = "term_name",
          value_col = "value",
          value_name = "-log10(p_val_adj)",
          region_col = "region",
          sex_col = if ("sex" %in% names(main_data)) "sex" else NULL,
          class_col = if ("class" %in% names(main_data)) "class" else NULL,
          neurotransmitter_col = NULL,
          row_go_col = "GO",
          row_module_col = if (use_row_category) "module" else NULL,
          row_category_col = NULL,
          col_order_priority = c("region", "sex", "class"),
          col_split_col = "region",
          cutoff = 10,
          distance = "euclidean",
          method = "complete", # "ward.D2" 也可
          logfc_cutoff = NULL,
          column_gap = grid::unit(3, "mm"),
          row_gap = grid::unit(3, "mm"),
          row_names_gp = 6,
          column_names_gp = 6,
          title = output_path,
          cluster_rows = cfg$rows,
          cluster_columns = cfg$cols,
          cluster_column_slices = FALSE,
          cluster_row_slices = FALSE
        )
      }
    }
  }
  
  message("Completed processing for: ", go_path)
}



######################small_onlyGO##############
for (go_path in go_paths) {
  message("Processing: ", go_path)
  
  go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1)
  go_results_tmp=go_results_tmp[go_results_tmp$Region!="iCTX",]
  go_results_tmp=go_results_tmp[go_results_tmp$source%in%c("GO:BP","GO:MF"),]
  #go_results_tmp$padj=p.adjust(go_results_tmp$p_value,method = "BH")
  go_results_tmp$nlogp=-log10(go_results_tmp$padj)*sign(go_results_tmp$nlogp)
  go_results_tmp$nlog10p=-log10(go_results_tmp$padj)
  
  go_results_tmp=go_results_tmp[go_results_tmp$padj<0.05,]
  colnames(go_results_tmp)=gsub("Region.Subclass","Region Subclass",colnames(go_results_tmp))
  go_results_tmp_unique <- go_results_tmp %>%
    distinct(`Region Subclass`, Sex, GO, .keep_all = TRUE)
  if (nrow(go_results_tmp_unique) < nrow(go_results_tmp)) {
    go_path_unique <- sub("\\.csv$", "_unique.csv", go_path)
    write.csv(go_results_tmp_unique, go_path_unique, row.names = FALSE)
  }
  
  # go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1) %>%
  #   filter(Sex %in% sexes)
  
  #go_results <- switch(target_data,
  #                     "SUS" = go_results_tmp_unique
  #)
  go_results=go_results_tmp_unique
  
  output_dir_path <- dirname(go_path)
  file_name <- tools::file_path_sans_ext(basename(go_path))
  #output_dir_path <- file.path(output_dir_path, file_name)
  output_dir_path <- file.path(output_dir_path, paste0(file_name,"_small_onlyGOadj"))
  
  if (!dir.exists(output_dir_path)) {
    dir.create(output_dir_path, recursive = TRUE)
  }
  
  subsets <- list(
    list(tag = "N", filter = function(df) df %>% filter(Neurotransmitter != "NN")),
    list(tag = "NN",    filter = function(df) df %>% filter(Neurotransmitter == "NN"))
  )
  
  group_col <- 'predicted_group' # 'Module.P','predicted_group' or 'final_group'
  use_row_category <-TRUE
  class_target <-  NULL # c("Microglia")
  # go_results$Module.P <- gsub("^G", "M", go_results[[group_col]])
  go_results <- go_results %>%
    mutate(!!group_col := ifelse(.data[[group_col]] == "others", "other", .data[[group_col]]))
  module_levels <- c("G1","G2-1","G2-2","G3-1","G3-2","G3-3","G4","G5","G6","G7","G7-1","G7-2","G8","G8-1","G8-2","G9","G10","G11","G12","other") # c("innate", "overlapped", "adaptive", "inflammatory")
  
  if (use_row_category) {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter,all_of(group_col)) %>%
      rename(module = all_of(group_col)) %>%
      {
        if (!is.null(module_levels)) {
          mutate(., module = factor(module, levels = module_levels))
        } else {
          mutate(., module = factor(module))
        }
      } %>%
      arrange(module)
  } else {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter)
  }
  
  
  # Generate heatmaps for each subset and clustering configuration
  # split N and NN for separate heatmaps
  # 生成不同子集与聚类配置的热图（包含：不过滤GO/列的unfiltered版本）
  # 提示：确保 dplyr 在最后加载，或统一使用 dplyr:: 前缀（本代码已使用前缀）
  # library(dplyr)  # 如果你不想每次都写 dplyr::，可把这句放到所有包加载的最后
  
  # split N and NN for separate heatmaps
  for (subset_cfg in subsets) {
    subset_tag <- subset_cfg$tag
    message("Generating subset: ", subset_tag)
    
    # ---- 构建主数据（列数据） ----
    go_results_subset <- subset_cfg$filter(go_results) %>%
      { if (!is.null(class_target)) dplyr::filter(., Class %in% class_target) else . } %>%
      dplyr::mutate(
        sample = paste(`Region Subclass`, Sex, sep = "_"),
        type   = ifelse(Neurotransmitter == "NN", "NN", "N")
      ) %>%
      dplyr::select(
        sample = sample,
        term_name = GO,
        value = nlogp,
        region = Region,
        sex = Sex,
        class = Class,
        type
      )
    
    # ---- 构建行分组（行数据） ----
    go_group_subset <- subset_cfg$filter(go_group)
    
    # ---- 预置三套过滤集合 ----
    filtered_sets <- list(
      common_main = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 3, row_col = "GO"), # >10
      common_sup  = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 1, row_col = "GO"),
      specific_sup = NULL
    )
    
    # ---- 新增：unfiltered（不过滤 GO term 或列）----
    filtered_sets <- c(
      list(unfiltered = list(
        main_data = go_results_subset,
        row_data  = go_group_subset
      )),
      filtered_sets
    )
    
    # ---- 基于 common_sup 计算 specific_sup ----
    common_sup_genes <- unique(filtered_sets$common_sup$row_data$GO)
    specific_sup_main <- go_results_subset %>%
      dplyr::filter(!term_name %in% common_sup_genes)
    specific_sup_row <- go_group_subset %>%
      dplyr::filter(!GO %in% common_sup_genes)
    filtered_sets$specific_sup <- list(
      main_data = specific_sup_main,
      row_data  = specific_sup_row
    )
    
    # ---- 聚类配置 ----
    cluster_cfgs <- list(
      list(rows = FALSE, cols = FALSE, tag = "nocluster"),
      list(rows = FALSE, cols = TRUE,  tag = "cluster_col"),
      list(rows = TRUE,  cols = FALSE, tag = "cluster_row"),
      list(rows = TRUE,  cols = TRUE,  tag = "cluster_row_col")
    )
    
    # ---- 逐集合、逐聚类配置画图 ----
    for (set_name in names(filtered_sets)) {
      data_pair <- filtered_sets[[set_name]]
      
      # 兼容列表结构或 data.frame 结构
      main_data <- if (is.list(data_pair) && !is.null(data_pair$main_data)) data_pair$main_data else data_pair
      row_data  <- if (is.list(data_pair) && !is.null(data_pair$row_data))  data_pair$row_data  else go_group_subset
      main_data[main_data$type=="N","class"]=gsub("(.* )|(_.*)","",main_data[main_data$type=="N","sample"])
      
      if (is.null(main_data) || nrow(main_data) == 0) {
        message("Skipping empty dataset: ", set_name)
        next
      }
      
      base_output <- file.path(
        output_dir_path,
        paste0("heatmap_", file_name, "_", subset_tag, "_", set_name)
      )
      
      for (cfg in cluster_cfgs) {
        output_path <- paste0(base_output, "_", cfg$tag, ".pdf")
        message("Generating heatmap with config: ", subset_tag, " + ", set_name, " + ", cfg$tag)
        
        results <- heatmap_final3(
          go = main_data,
          row = row_data,
          sample_col = "sample",
          term_col = "term_name",
          value_col = "value",
          value_name = "-log10(p_val_adj)",
          region_col = "region",
          sex_col = if ("sex" %in% names(main_data)) "sex" else NULL,
          class_col = if ("class" %in% names(main_data)) "class" else NULL,
          neurotransmitter_col = NULL,
          row_go_col = "GO",
          row_module_col = if (use_row_category) "module" else NULL,
          row_category_col = NULL,
          col_order_priority = c("region", "sex", "class"),
          col_split_col = "region",
          cutoff = 10,
          distance = "euclidean",
          method = "complete", # "ward.D2" 也可
          logfc_cutoff = NULL,
          column_gap = grid::unit(3, "mm"),
          row_gap = grid::unit(3, "mm"),
          row_names_gp = 6,
          column_names_gp = 6,
          title = output_path,
          cluster_rows = cfg$rows,
          cluster_columns = cfg$cols,
          cluster_column_slices = FALSE,
          cluster_row_slices = FALSE
        )
      }
    }
  }
  
  message("Completed processing for: ", go_path)
}
#########################onlyGo###########
for (go_path in go_paths) {
  message("Processing: ", go_path)
  
  go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1)
  go_results_tmp=go_results_tmp[go_results_tmp$Region!="iCTX",]
  go_results_tmp=go_results_tmp[go_results_tmp$source%in%c("GO:BP","GO:MF"),]
  #go_results_tmp$padj=p.adjust(go_results_tmp$p_value,method = "BH")
  go_results_tmp$nlogp=-log10(go_results_tmp$padj)*sign(go_results_tmp$nlogp)
    go_results_tmp$nlog10p=-log10(go_results_tmp$padj)
  
  go_results_tmp=go_results_tmp[go_results_tmp$padj<0.05,]
  colnames(go_results_tmp)=gsub("Region.Subclass","Region Subclass",colnames(go_results_tmp))
  go_results_tmp_unique <- go_results_tmp %>%
    distinct(`Region Subclass`, Sex, GO, .keep_all = TRUE)
  if (nrow(go_results_tmp_unique) < nrow(go_results_tmp)) {
    go_path_unique <- sub("\\.csv$", "_unique.csv", go_path)
    write.csv(go_results_tmp_unique, go_path_unique, row.names = FALSE)
  }
  
  # go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1) %>%
  #   filter(Sex %in% sexes)
  
  #go_results <- switch(target_data,
  #                     "SUS" = go_results_tmp_unique
  #)
  go_results=go_results_tmp_unique
  
  output_dir_path <- dirname(go_path)
  file_name <- tools::file_path_sans_ext(basename(go_path))
  #output_dir_path <- file.path(output_dir_path, file_name)
  output_dir_path <- file.path(output_dir_path, paste0(file_name,"_onlyGOadj"))
  
  if (!dir.exists(output_dir_path)) {
    dir.create(output_dir_path, recursive = TRUE)
  }
  
  subsets <- list(
    list(tag = "N", filter = function(df) df %>% filter(Neurotransmitter != "NN")),
    list(tag = "NN",    filter = function(df) df %>% filter(Neurotransmitter == "NN"))
  )
  
  group_col <- 'predicted_group' # 'Module.P','predicted_group' or 'final_group'
  use_row_category <-TRUE
  class_target <-  NULL # c("Microglia")
  # go_results$Module.P <- gsub("^G", "M", go_results[[group_col]])
  go_results <- go_results %>%
    mutate(!!group_col := ifelse(.data[[group_col]] == "others", "other", .data[[group_col]]))
  module_levels <- c("G1","G2-1","G2-2","G3-1","G3-2","G3-3","G4","G5","G6","G7","G7-1","G7-2","G8","G8-1","G8-2","G9","G10","G11","G12","other") # c("innate", "overlapped", "adaptive", "inflammatory")
  
  if (use_row_category) {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter,all_of(group_col)) %>%
      rename(module = all_of(group_col)) %>%
      {
        if (!is.null(module_levels)) {
          mutate(., module = factor(module, levels = module_levels))
        } else {
          mutate(., module = factor(module))
        }
      } %>%
      arrange(module)
  } else {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter)
  }
  
  
  # Generate heatmaps for each subset and clustering configuration
  # split N and NN for separate heatmaps
  # 生成不同子集与聚类配置的热图（包含：不过滤GO/列的unfiltered版本）
  # 提示：确保 dplyr 在最后加载，或统一使用 dplyr:: 前缀（本代码已使用前缀）
  # library(dplyr)  # 如果你不想每次都写 dplyr::，可把这句放到所有包加载的最后
  
  # split N and NN for separate heatmaps
  for (subset_cfg in subsets) {
    subset_tag <- subset_cfg$tag
    message("Generating subset: ", subset_tag)
    
    # ---- 构建主数据（列数据） ----
    go_results_subset <- subset_cfg$filter(go_results) %>%
      { if (!is.null(class_target)) dplyr::filter(., Class %in% class_target) else . } %>%
      dplyr::mutate(
        sample = paste(`Region Subclass`, Sex, sep = "_"),
        type   = ifelse(Neurotransmitter == "NN", "NN", "N")
      ) %>%
      dplyr::select(
        sample = sample,
        term_name = GO,
        value = nlogp,
        region = Region,
        sex = Sex,
        class = Class,
        type
      )
    
    # ---- 构建行分组（行数据） ----
    go_group_subset <- subset_cfg$filter(go_group)
    
    # ---- 预置三套过滤集合 ----
    filtered_sets <- list(
      common_main = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 3, row_col = "GO"), # >10
      common_sup  = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 1, row_col = "GO"),
      specific_sup = NULL
    )
    
    # ---- 新增：unfiltered（不过滤 GO term 或列）----
    filtered_sets <- c(
      list(unfiltered = list(
        main_data = go_results_subset,
        row_data  = go_group_subset
      )),
      filtered_sets
    )
    
    # ---- 基于 common_sup 计算 specific_sup ----
    common_sup_genes <- unique(filtered_sets$common_sup$row_data$GO)
    specific_sup_main <- go_results_subset %>%
      dplyr::filter(!term_name %in% common_sup_genes)
    specific_sup_row <- go_group_subset %>%
      dplyr::filter(!GO %in% common_sup_genes)
    filtered_sets$specific_sup <- list(
      main_data = specific_sup_main,
      row_data  = specific_sup_row
    )
    
    # ---- 聚类配置 ----
    cluster_cfgs <- list(
      list(rows = FALSE, cols = FALSE, tag = "nocluster"),
      list(rows = FALSE, cols = TRUE,  tag = "cluster_col"),
      list(rows = TRUE,  cols = FALSE, tag = "cluster_row"),
      list(rows = TRUE,  cols = TRUE,  tag = "cluster_row_col")
    )
    
    # ---- 逐集合、逐聚类配置画图 ----
    for (set_name in names(filtered_sets)) {
      data_pair <- filtered_sets[[set_name]]
      
      # 兼容列表结构或 data.frame 结构
      main_data <- if (is.list(data_pair) && !is.null(data_pair$main_data)) data_pair$main_data else data_pair
      row_data  <- if (is.list(data_pair) && !is.null(data_pair$row_data))  data_pair$row_data  else go_group_subset
      main_data[main_data$type=="N","class"]=gsub("(.* )|(_.*)","",main_data[main_data$type=="N","sample"])
      
      if (is.null(main_data) || nrow(main_data) == 0) {
        message("Skipping empty dataset: ", set_name)
        next
      }
      
      base_output <- file.path(
        output_dir_path,
        paste0("heatmap_", file_name, "_", subset_tag, "_", set_name)
      )
      
      for (cfg in cluster_cfgs) {
        output_path <- paste0(base_output, "_", cfg$tag, ".pdf")
        message("Generating heatmap with config: ", subset_tag, " + ", set_name, " + ", cfg$tag)
        
        results <- heatmap_final2(
          go = main_data,
          row = row_data,
          sample_col = "sample",
          term_col = "term_name",
          value_col = "value",
          value_name = "-log10(p_val_adj)",
          region_col = "region",
          sex_col = if ("sex" %in% names(main_data)) "sex" else NULL,
          class_col = if ("class" %in% names(main_data)) "class" else NULL,
          neurotransmitter_col = NULL,
          row_go_col = "GO",
          row_module_col = if (use_row_category) "module" else NULL,
          row_category_col = NULL,
          col_order_priority = c("region", "sex", "class"),
          col_split_col = "region",
          cutoff = 10,
          distance = "euclidean",
          method = "complete", # "ward.D2" 也可
          logfc_cutoff = NULL,
          column_gap = grid::unit(3, "mm"),
          row_gap = grid::unit(3, "mm"),
          row_names_gp = 6,
          column_names_gp = 6,
          title = output_path,
          cluster_rows = cfg$rows,
          cluster_columns = cfg$cols,
          cluster_column_slices = FALSE,
          cluster_row_slices = FALSE
        )
      }
    }
  }
  
  message("Completed processing for: ", go_path)
}

#########################alldatabase###########
for (go_path in go_paths) {
  message("Processing: ", go_path)
  
  go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1)
  go_results_tmp=go_results_tmp[go_results_tmp$Region!="iCTX",]
  go_results_tmp=go_results_tmp[go_results_tmp$p_value<0.01,]
  #go_results_tmp=go_results_tmp[go_results_tmp$source%in%c("GO:BP","GO:MF"),]
  colnames(go_results_tmp)=gsub("Region.Subclass","Region Subclass",colnames(go_results_tmp))
  go_results_tmp_unique <- go_results_tmp %>%
    distinct(`Region Subclass`, Sex, GO, .keep_all = TRUE)
  if (nrow(go_results_tmp_unique) < nrow(go_results_tmp)) {
    go_path_unique <- sub("\\.csv$", "_unique.csv", go_path)
    write.csv(go_results_tmp_unique, go_path_unique, row.names = FALSE)
  }
  
  # go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1) %>%
  #   filter(Sex %in% sexes)
  
  #go_results <- switch(target_data,
  #                     "SUS" = go_results_tmp_unique
  #)
  go_results=go_results_tmp_unique
  
  output_dir_path <- dirname(go_path)
  file_name <- tools::file_path_sans_ext(basename(go_path))
  #output_dir_path <- file.path(output_dir_path, file_name)
  output_dir_path <- file.path(output_dir_path, paste0(file_name,"_alldatabase"))
  
  if (!dir.exists(output_dir_path)) {
    dir.create(output_dir_path, recursive = TRUE)
  }
  
  subsets <- list(
    list(tag = "N", filter = function(df) df %>% filter(Neurotransmitter != "NN")),
    list(tag = "NN",    filter = function(df) df %>% filter(Neurotransmitter == "NN"))
  )
  
  group_col <- 'predicted_group' # 'Module.P','predicted_group' or 'final_group'
  use_row_category <-TRUE
  class_target <-  NULL # c("Microglia")
  # go_results$Module.P <- gsub("^G", "M", go_results[[group_col]])
  go_results <- go_results %>%
    mutate(!!group_col := ifelse(.data[[group_col]] == "others", "other", .data[[group_col]]))
  module_levels <- c("G1","G2-1","G2-2","G3-1","G3-2","G3-3","G4","G5","G6","G7","G7-1","G7-2","G8","G8-1","G8-2","G9","G10","G11","G12","other") # c("innate", "overlapped", "adaptive", "inflammatory")
  
  if (use_row_category) {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter,all_of(group_col)) %>%
      rename(module = all_of(group_col)) %>%
      {
        if (!is.null(module_levels)) {
          mutate(., module = factor(module, levels = module_levels))
        } else {
          mutate(., module = factor(module))
        }
      } %>%
      arrange(module)
  } else {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter)
  }
  
  
  # Generate heatmaps for each subset and clustering configuration
  # split N and NN for separate heatmaps
  # 生成不同子集与聚类配置的热图（包含：不过滤GO/列的unfiltered版本）
  # 提示：确保 dplyr 在最后加载，或统一使用 dplyr:: 前缀（本代码已使用前缀）
  # library(dplyr)  # 如果你不想每次都写 dplyr::，可把这句放到所有包加载的最后
  
  # split N and NN for separate heatmaps
  for (subset_cfg in subsets) {
    subset_tag <- subset_cfg$tag
    message("Generating subset: ", subset_tag)
    
    # ---- 构建主数据（列数据） ----
    go_results_subset <- subset_cfg$filter(go_results) %>%
      { if (!is.null(class_target)) dplyr::filter(., Class %in% class_target) else . } %>%
      dplyr::mutate(
        sample = paste(`Region Subclass`, Sex, sep = "_"),
        type   = ifelse(Neurotransmitter == "NN", "NN", "N")
      ) %>%
      dplyr::select(
        sample = sample,
        term_name = GO,
        value = nlogp,
        region = Region,
        sex = Sex,
        class = Class,
        type
      )
    
    # ---- 构建行分组（行数据） ----
    go_group_subset <- subset_cfg$filter(go_group)
    
    # ---- 预置三套过滤集合 ----
    filtered_sets <- list(
      common_main = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 3, row_col = "GO"), # >10
      common_sup  = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 1, row_col = "GO"),
      specific_sup = NULL
    )
    
    # ---- 新增：unfiltered（不过滤 GO term 或列）----
    filtered_sets <- c(
      list(unfiltered = list(
        main_data = go_results_subset,
        row_data  = go_group_subset
      )),
      filtered_sets
    )
    
    # ---- 基于 common_sup 计算 specific_sup ----
    common_sup_genes <- unique(filtered_sets$common_sup$row_data$GO)
    specific_sup_main <- go_results_subset %>%
      dplyr::filter(!term_name %in% common_sup_genes)
    specific_sup_row <- go_group_subset %>%
      dplyr::filter(!GO %in% common_sup_genes)
    filtered_sets$specific_sup <- list(
      main_data = specific_sup_main,
      row_data  = specific_sup_row
    )
    
    # ---- 聚类配置 ----
    cluster_cfgs <- list(
      list(rows = FALSE, cols = FALSE, tag = "nocluster"),
      list(rows = FALSE, cols = TRUE,  tag = "cluster_col"),
      list(rows = TRUE,  cols = FALSE, tag = "cluster_row"),
      list(rows = TRUE,  cols = TRUE,  tag = "cluster_row_col")
    )
    
    # ---- 逐集合、逐聚类配置画图 ----
    for (set_name in names(filtered_sets)) {
      data_pair <- filtered_sets[[set_name]]
      
      # 兼容列表结构或 data.frame 结构
      main_data <- if (is.list(data_pair) && !is.null(data_pair$main_data)) data_pair$main_data else data_pair
      row_data  <- if (is.list(data_pair) && !is.null(data_pair$row_data))  data_pair$row_data  else go_group_subset
      main_data[main_data$type=="N","class"]=gsub("(.* )|(_.*)","",main_data[main_data$type=="N","sample"])
      
      if (is.null(main_data) || nrow(main_data) == 0) {
        message("Skipping empty dataset: ", set_name)
        next
      }
      
      base_output <- file.path(
        output_dir_path,
        paste0("heatmap_", file_name, "_", subset_tag, "_", set_name)
      )
      
      for (cfg in cluster_cfgs) {
        output_path <- paste0(base_output, "_", cfg$tag, ".pdf")
        message("Generating heatmap with config: ", subset_tag, " + ", set_name, " + ", cfg$tag)
        
        results <- heatmap_final2(
          go = main_data,
          row = row_data,
          sample_col = "sample",
          term_col = "term_name",
          value_col = "value",
          value_name = "-log10(p_val_adj)",
          region_col = "region",
          sex_col = if ("sex" %in% names(main_data)) "sex" else NULL,
          class_col = if ("class" %in% names(main_data)) "class" else NULL,
          neurotransmitter_col = NULL,
          row_go_col = "GO",
          row_module_col = if (use_row_category) "module" else NULL,
          row_category_col = NULL,
          col_order_priority = c("region", "sex", "class"),
          col_split_col = "region",
          cutoff = 10,
          distance = "euclidean",
          method = "complete", # "ward.D2" 也可
          logfc_cutoff = NULL,
          column_gap = grid::unit(3, "mm"),
          row_gap = grid::unit(3, "mm"),
          row_names_gp = 6,
          column_names_gp = 6,
          title = output_path,
          cluster_rows = cfg$rows,
          cluster_columns = cfg$cols,
          cluster_column_slices = FALSE,
          cluster_row_slices = FALSE
        )
      }
    }
  }
  
  message("Completed processing for: ", go_path)
}

##############################################alldatabase_small
for (go_path in go_paths) {
  message("Processing: ", go_path)
  
  go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1)
  go_results_tmp[go_results_tmp$Region!="iCTX",]
  colnames(go_results_tmp)=gsub("Region.Subclass","Region Subclass",colnames(go_results_tmp))
  
  go_results_tmp_unique <- go_results_tmp %>%
    distinct(`Region Subclass`, Sex, GO, .keep_all = TRUE)
  if (nrow(go_results_tmp_unique) < nrow(go_results_tmp)) {
    go_path_unique <- sub("\\.csv$", "_unique.csv", go_path)
    write.csv(go_results_tmp_unique, go_path_unique, row.names = FALSE)
  }
  
  # go_results_tmp <- read_GO_DEG(go_path, sex_col = "Sex", sheet = 1) %>%
  #   filter(Sex %in% sexes)
  
  #go_results <- switch(target_data,
  #                     "SUS" = go_results_tmp_unique
  #)
  go_results=go_results_tmp_unique
  
  output_dir_path <- dirname(go_path)
  file_name <- tools::file_path_sans_ext(basename(go_path))
  output_dir_path <- file.path(output_dir_path, file_name)
  output_dir_path1 <- file.path(output_dir_path, paste0(file_name,"_alldatabasesmall"))
  
  if (!dir.exists(output_dir_path)) {
    dir.create(output_dir_path, recursive = TRUE)
  }
  
  subsets <- list(
    list(tag = "N", filter = function(df) df %>% filter(Neurotransmitter != "NN")),
    list(tag = "NN",    filter = function(df) df %>% filter(Neurotransmitter == "NN"))
  )
  
  group_col <- 'predicted_group' # 'Module.P','predicted_group' or 'final_group'
  use_row_category <-TRUE
  class_target <-  NULL # c("Microglia")
  # go_results$Module.P <- gsub("^G", "M", go_results[[group_col]])
  go_results <- go_results %>%
    mutate(!!group_col := ifelse(.data[[group_col]] == "others", "other", .data[[group_col]]))
  module_levels <- c("G1","G2-1","G2-2","G3-1","G3-2","G3-3","G4","G5","G6","G7","G7-1","G7-2","G8","G8-1","G8-2","G9","G10","G11","G12","other") # c("innate", "overlapped", "adaptive", "inflammatory")
  
  if (use_row_category) {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter,all_of(group_col)) %>%
      rename(module = all_of(group_col)) %>%
      {
        if (!is.null(module_levels)) {
          mutate(., module = factor(module, levels = module_levels))
        } else {
          mutate(., module = factor(module))
        }
      } %>%
      arrange(module)
  } else {
    go_group <- go_results %>%
      dplyr::select(GO,Neurotransmitter)
  }
  
  
  # Generate heatmaps for each subset and clustering configuration
  # split N and NN for separate heatmaps
  # 生成不同子集与聚类配置的热图（包含：不过滤GO/列的unfiltered版本）
  # 提示：确保 dplyr 在最后加载，或统一使用 dplyr:: 前缀（本代码已使用前缀）
  # library(dplyr)  # 如果你不想每次都写 dplyr::，可把这句放到所有包加载的最后
  
  # split N and NN for separate heatmaps
  for (subset_cfg in subsets) {
    subset_tag <- subset_cfg$tag
    message("Generating subset: ", subset_tag)
    
    # ---- 构建主数据（列数据） ----
    go_results_subset <- subset_cfg$filter(go_results) %>%
      { if (!is.null(class_target)) dplyr::filter(., Class %in% class_target) else . } %>%
      dplyr::mutate(
        sample = paste(`Region Subclass`, Sex, sep = "_"),
        type   = ifelse(Neurotransmitter == "NN", "NN", "N")
      ) %>%
      dplyr::select(
        sample = sample,
        term_name = GO,
        value = nlogp,
        region = Region,
        sex = Sex,
        class = Class,
        type
      )
    
    # ---- 构建行分组（行数据） ----
    go_group_subset <- subset_cfg$filter(go_group)
    
    # ---- 预置三套过滤集合 ----
    filtered_sets <- list(
      common_main = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 3, row_col = "GO"), # >10
      common_sup  = filter_col_and_row(go_results_subset, go_group_subset,
                                       n_cells = 3, n_genes = 1, row_col = "GO"),
      specific_sup = NULL
    )
    
    # ---- 新增：unfiltered（不过滤 GO term 或列）----
    filtered_sets <- c(
      list(unfiltered = list(
        main_data = go_results_subset,
        row_data  = go_group_subset
      )),
      filtered_sets
    )
    
    # ---- 基于 common_sup 计算 specific_sup ----
    common_sup_genes <- unique(filtered_sets$common_sup$row_data$GO)
    specific_sup_main <- go_results_subset %>%
      dplyr::filter(!term_name %in% common_sup_genes)
    specific_sup_row <- go_group_subset %>%
      dplyr::filter(!GO %in% common_sup_genes)
    filtered_sets$specific_sup <- list(
      main_data = specific_sup_main,
      row_data  = specific_sup_row
    )
    
    # ---- 聚类配置 ----
    cluster_cfgs <- list(
      list(rows = FALSE, cols = FALSE, tag = "nocluster"),
      list(rows = FALSE, cols = TRUE,  tag = "cluster_col"),
      list(rows = TRUE,  cols = FALSE, tag = "cluster_row"),
      list(rows = TRUE,  cols = TRUE,  tag = "cluster_row_col")
    )
    
    # ---- 逐集合、逐聚类配置画图 ----
    for (set_name in names(filtered_sets)) {
      data_pair <- filtered_sets[[set_name]]
      
      # 兼容列表结构或 data.frame 结构
      main_data <- if (is.list(data_pair) && !is.null(data_pair$main_data)) data_pair$main_data else data_pair
      row_data  <- if (is.list(data_pair) && !is.null(data_pair$row_data))  data_pair$row_data  else go_group_subset
      main_data[main_data$type=="N","class"]=gsub("(.* )|(_.*)","",main_data[main_data$type=="N","sample"])
      
      if (is.null(main_data) || nrow(main_data) == 0) {
        message("Skipping empty dataset: ", set_name)
        next
      }
      
      base_output <- file.path(
        output_dir_path,
        paste0("heatmap_", file_name, "_", subset_tag, "_", set_name)
      )
      
      for (cfg in cluster_cfgs) {
        output_path <- paste0(base_output, "_", cfg$tag, ".pdf")
        message("Generating heatmap with config: ", subset_tag, " + ", set_name, " + ", cfg$tag)
        
        results <- heatmap_final2(
          go = main_data,
          row = row_data,
          sample_col = "sample",
          term_col = "term_name",
          value_col = "value",
          value_name = "-log10(p_val_adj)",
          region_col = "region",
          sex_col = if ("sex" %in% names(main_data)) "sex" else NULL,
          class_col = if ("class" %in% names(main_data)) "class" else NULL,
          neurotransmitter_col = NULL,
          row_go_col = "GO",
          row_module_col = if (use_row_category) "module" else NULL,
          row_category_col = NULL,
          col_order_priority = c("region", "sex", "class"),
          col_split_col = "region",
          cutoff = 10,
          distance = "euclidean",
          method = "complete", # "ward.D2" 也可
          logfc_cutoff = NULL,
          column_gap = grid::unit(3, "mm"),
          row_gap = grid::unit(3, "mm"),
          row_names_gp = 6,
          column_names_gp = 6,
          title = output_path,
          cluster_rows = cfg$rows,
          cluster_columns = cfg$cols,
          cluster_column_slices = FALSE,
          cluster_row_slices = FALSE
        )
      }
    }
  }
  
  message("Completed processing for: ", go_path)
}

