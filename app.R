#################################################
# Title: Single Cell Explorer
# Author: Emily Schueddig and Allison Makovec
# Date: 2/6/2026
# Last modified: 2/11/2026
#################################################
#Changes example

## Example text changes
# more chang4w

# app.R
options(shiny.maxRequestSize = 8000*1024^3)

suppressPackageStartupMessages({
  library(shiny)
  library(shinycssloaders)
  library(shinythemes)
  library(Seurat)
  library(ggplot2)
  library(plotly)
  library(DT)
  library(dplyr)
  library(stringr)
  library(patchwork)
  library(pals)
})

# ---------- Utilities ----------
is_seurat <- function(x) inherits(x, "Seurat")

available_reductions <- function(obj) {
  if (!is_seurat(obj)) return(character())
  reds <- Reductions(obj)
  if (is.null(reds)) character() else reds
}

available_assays <- function(obj) {
  if (!is_seurat(obj)) return(character())
  asy <- Assays(obj)
  if (is.null(asy)) character() else asy
}

ensure_basic_counts <- function(obj, assay = NULL) {
  stopifnot(is_seurat(obj))
  if (is.null(assay)) assay <- DefaultAssay(obj)
  ncount_col   <- paste0("nCount_", assay)
  nfeature_col <- paste0("nFeature_", assay)
  
  need_counts  <- !ncount_col %in% colnames(obj@meta.data)
  need_feats   <- !nfeature_col %in% colnames(obj@meta.data)
  
  if (need_counts || need_feats) {
    counts <- tryCatch(GetAssayData(obj, assay = assay, layer = "counts"), error = function(e) NULL)
    if (!is.null(counts)) {
      if (need_counts)  obj[[ncount_col]]   <- Matrix::colSums(counts)
      if (need_feats)   obj[[nfeature_col]] <- Matrix::colSums(counts > 0)
    }
  }
  obj
}

cat_meta_cols <- function(obj, exclude_too_unique = TRUE, max_unique_ratio = 0.8) {
  md <- obj@meta.data
  cand <- names(md)[!sapply(md, is.numeric)]
  if (!exclude_too_unique) return(cand)
  n <- nrow(md)
  Filter(function(nm) {
    nunique <- dplyr::n_distinct(md[[nm]])
    nunique / n < max_unique_ratio
  }, cand)
}

# Save last plot for download handlers
last_plot_env <- new.env(parent = emptyenv())
set_last_plot <- function(p, key) assign(key, p, envir = last_plot_env)
get_last_plot <- function(key) get(key, envir = last_plot_env, inherits = FALSE)

# Define UI ------------------------------------------------------------- 
ui <- fluidPage(
  theme = shinytheme("flatly"),
  # Create navigation bar -----------------------------------------------
  navbarPage(
    title = "Single-cell Explorer",
    # Create app section -----------------------------------------------
    tabPanel(
      "Explore the Data",
      titlePanel(h1("Explore the Data")),
      p(strong("Author:"), " Emily Schueddig | Department of Biostatistics & Data Science, University of Kansas Medical Center"),
      p(strong("Email:"), " eschueddig@kumc.edu"),
      sidebarLayout(
        # Side bar ------------------------------------------------------------
        sidebarPanel(
          width = 3,
          fileInput(
            "seurat_file",
            "Upload Seurat Object",
            accept = c(".rds", ".RData")
          ),
          tags$hr()
        ),
        
        mainPanel(
          width = 9,
          
          tabsetPanel(
            
            tabPanel(
              "Overview",
              fluidRow(
                column(3, wellPanel(h4("Cells"), uiOutput("nCells"))),
                column(3, wellPanel(h4("Genes"), uiOutput("nGenes"))),
                column(3, wellPanel(h4("Assays"), uiOutput("nAssays"))),
                column(3, wellPanel(h4("Reductions"), uiOutput("nReductions")))
              ),
              tags$hr(),
              h4("Metadata preview"),
              DTOutput("meta_preview") %>% withSpinner()
            ),
            
            tabPanel(
              "Reductions",
              uiOutput("dimplot_reduction_ui"),
              uiOutput("dimplot_group_ui"),
              uiOutput("dimplot_splitby_ui"),
              plotOutput("dimplot_reduction", height = "600px") %>% withSpinner(),
              downloadButton("download_dimplot", "Download DimPlot")
            ),
            
            tabPanel(
              "Feature Plot",
              uiOutput("featureplot_gene_ui"),
              uiOutput("featureplot_reduction_ui"),
              uiOutput("featureplot_splitby_ui"),
              plotOutput("feature_plot", height = "600px") %>% withSpinner(),
              downloadButton("download_feature_plot", "Download FeaturePlot")
            ),
            
            tabPanel(
              "Cluster Composition",
              uiOutput("comp_group_ui"),
              uiOutput("comp_stack_ui"),
              plotlyOutput("composition_plot") %>% withSpinner(),
              downloadButton("download_comp_plot", "Download composition plot")
            ),
            
            tabPanel(
              "Gene Violins",
              textInput("genes_violin", "Genes (comma-separated)", value = ""),
              uiOutput("violin_group_ui"),
              plotlyOutput("gene_violin_plot") %>% withSpinner(),
              downloadButton("download_gene_violin", "Download gene violins"),
              tags$hr(),
              DTOutput("gene_violin_summary") %>% withSpinner()
            ),
            
            tabPanel(
              "Filtered Gene Violins",
              uiOutput("filter_var_ui"),
              uiOutput("filter_value_ui"),
              uiOutput("filtered_violin_group_ui"),
              textInput("filtered_genes", "Genes (comma-separated)", value = ""),
              plotlyOutput("filtered_violin_plot") %>% withSpinner(),
              downloadButton("download_filtered_violins", "Download violin plot"),
              tags$hr(),
              DTOutput("filtered_violin_summary") %>% withSpinner()
            )
            
          )
        )
      )
    ),
    tabPanel(
      "Analyze the Data",
      titlePanel(h1("Analyze the Data"))
    )
  )
)

# ---------- Server ----------
server <- function(input, output, session) {
  
  seurat_obj <- reactiveVal(NULL)
  
  observeEvent(input$seurat_file, {
    req(input$seurat_file)
    withProgress(message = "Loading Seurat object...", value = 0.1, {
      obj <- readRDS(input$seurat_file$datapath)
      validate(need(is_seurat(obj), "The file you uploaded is not a Seurat object."))
      incProgress(0.3, detail = "Preparing metadata & choices")
      seurat_obj(obj)
      # Update UI choices
      updateSelectInput(session, "assay", choices = available_assays(obj), selected = DefaultAssay(obj))
      reds <- available_reductions(obj)
      # sel_red <- if ("umap" %in% reds) "umap" else if (length(reds) > 0) reds[[1]] else NULL
      updateSelectInput(session, "dimplot_reduction", choices = reds, selected = NULL)
      md <- obj@meta.data
      cat_cols <- cat_meta_cols(obj)
      default_group <- if ("seurat_clusters" %in% colnames(md)) "seurat_clusters" else ""
      updateSelectInput(session, "groupby", choices = c("(none)" = "", cat_cols), selected = default_group)
      # Genes
      updateSelectizeInput(session, "gene", choices = rownames(obj), server = TRUE)
      incProgress(0.6, detail = "Done")
    })
  })
  
  # Assay / reductions / group-by UI
  output$assay_ui <- renderUI({
    req(seurat_obj())
    selectInput("assay", "Assay", choices = available_assays(seurat_obj()), selected = DefaultAssay(seurat_obj()))
  })

  output$dimplot_reduction_ui <- renderUI({
    req(seurat_obj())
    selectInput("dimplot_reduction", "Dimensional Reduction",
                choices = available_reductions(seurat_obj()),
                selected = NULL)
  })
  
  output$dimplot_group_ui <- renderUI({
    req(seurat_obj())
    choices <- cat_meta_cols(seurat_obj())
    selectInput("dimplot_group", "Group by (color)", choices = choices, selected = "seurat_clusters")
  })
  
  output$dimplot_splitby_ui <- renderUI({
    req(seurat_obj())
    choices <- c("(none)" = "", cat_meta_cols(seurat_obj()))
    selectInput("dimplot_splitby", "Split by (optional)", choices = choices, selected = "")
  })
  
  # UI: Gene input
  output$featureplot_gene_ui <- renderUI({
    req(seurat_obj())
    selectizeInput("featureplot_gene", "Select Gene",
                   choices = rownames(seurat_obj()),
                   selected = NULL,
                   options = list(placeholder = 'Start typing a gene...', maxOptions = 2000))
  })
  
  # UI: Dimensional reduction input
  output$featureplot_reduction_ui <- renderUI({
    req(seurat_obj())
    selectInput("featureplot_reduction", "Reduction method",
                choices = available_reductions(seurat_obj()),
                selected = "umap")
  })
  
  # UI: Optional split-by input
  output$featureplot_splitby_ui <- renderUI({
    req(seurat_obj())
    choices <- c("(none)" = "", cat_meta_cols(seurat_obj()))
    selectInput("featureplot_splitby", "Split by (optional)", choices = choices, selected = "")
  })
  
  output$gene_ui <- renderUI({
    req(seurat_obj())
    selectizeInput("gene", "Gene (for FeaturePlot)", choices = NULL, selected = NULL, multiple = FALSE,
                   options = list(placeholder = 'Start typing a gene...', maxOptions = 2000))
  })
  
  output$violin_group_ui <- renderUI({
    req(seurat_obj())
    selectInput("violin_group", "Group by", choices = cat_meta_cols(seurat_obj()), selected = "seurat_clusters")
  })
  
  # Dropdown to select variable to filter by
  output$filter_var_ui <- renderUI({
    req(seurat_obj())
    choices <- cat_meta_cols(seurat_obj(), exclude_too_unique = TRUE)
    validate(need(length(choices) > 0, "No suitable categorical variables found in metadata."))
    selectInput("filter_var", "Filter by variable", choices = choices, selected = choices[[1]])
  })
  
  # Dropdown to select **one** value from the selected variable
  output$filter_value_ui <- renderUI({
    req(seurat_obj(), input$filter_var)
    values <- seurat_obj()@meta.data[[input$filter_var]]
    
    # Ensure it's character for display
    if (is.factor(values)) values <- as.character(values)
    
    lvls <- sort(unique(values))
    selectInput("filter_value", "Value to keep", choices = lvls, selected = lvls[[1]])
  })
  
  # Grouping variable UI
  output$filtered_violin_group_ui <- renderUI({
    req(seurat_obj())
    selectInput("filtered_violin_group", "Group by", choices = cat_meta_cols(seurat_obj()), selected = "seurat_clusters")
  })
  
  # Overview tiles
  output$nCells <- renderUI({ req(seurat_obj()); tags$h3(format(ncol(seurat_obj()), big.mark=",")) })
  output$nGenes <- renderUI({ req(seurat_obj()); tags$h3(format(nrow(seurat_obj()), big.mark=",")) })
  output$nAssays <- renderUI({ req(seurat_obj()); tags$h3(length(available_assays(seurat_obj()))) })
  output$nReductions <- renderUI({ req(seurat_obj()); tags$h3(length(available_reductions(seurat_obj()))) })
  
  # Metadata preview
  output$meta_preview <- renderDT({
    req(seurat_obj())
    datatable(head(seurat_obj()@meta.data, 20), extensions = "Buttons",
              options = list(scrollX = TRUE, dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')))
  })
  
  groupby_or_ident <- reactive({
    if (!is.null(input$groupby) && nzchar(input$groupby)) input$groupby else NULL
  })
  
  assay_metrics <- reactive({
    req(seurat_obj())
    assay <- if (!is.null(input$assay) && nzchar(input$assay)) input$assay else DefaultAssay(seurat_obj())
    ncount_col   <- paste0("nCount_", assay)
    nfeature_col <- paste0("nFeature_", assay)
    list(assay = assay, ncount = ncount_col, nfeature = nfeature_col)
  })
  
  # ---- Reduction ----
  output$dimplot_reduction <- renderPlot({
    req(seurat_obj(), input$dimplot_reduction, input$dimplot_group)

    split_var <- if (nzchar(input$dimplot_splitby)) input$dimplot_splitby else NULL
    
    p <- DimPlot(seurat_obj(),
                 reduction = input$dimplot_reduction,
                 group.by = input$dimplot_group,
                 split.by = split_var,
                 pt.size = 1, ncol=2, label = T)
    
    set_last_plot(p, "dimplot_raw")
    print(p)
  })
  
  output$download_dimplot <- downloadHandler(
    filename = function() sprintf("DimPlot_%s.png", Sys.Date()),
    content = function(file) {
      req(exists("dimplot_raw", envir = last_plot_env))
      ggsave(file, plot = get_last_plot("dimplot_raw"), width = 8, height = 6, dpi = 300)
    }
  )
  
  # ----Feature Plot ----
  output$feature_plot <- renderPlot({
    req(seurat_obj(), input$featureplot_gene, input$featureplot_reduction)
    
    split_var <- if (nzchar(input$featureplot_splitby)) input$featureplot_splitby else NULL
    
    p_list <- FeaturePlot(seurat_obj(),
                     features = input$featureplot_gene,
                     reduction = input$featureplot_reduction,
                     split.by = split_var,
                     pt.size = 0.5, order = T, ncol = 2) 
    
    
    p_list <- lapply(p_list, function(p) {
      p + theme_minimal() + scale_color_gradientn(colors = c("#D1E5F0", "#FEE090", "#B2182B", "#67001F"))
    })
    
    p <- patchwork::wrap_plots(p_list, ncol = 2)
    print(p)

    print(p)
  })
  
  output$download_feature_plot <- downloadHandler(
    filename = function() sprintf("FeaturePlot_%s.png", Sys.Date()),
    content = function(file) {
      req(exists("featureplot_raw", envir = last_plot_env))
      ggsave(file, plot = get_last_plot("featureplot_raw"), width = 8, height = 6, dpi = 300)
    }
  )
  
  # ---- Composition ----
  output$comp_group_ui <- renderUI({
    req(seurat_obj())
    selectInput("comp_group", "Primary grouping", choices = cat_meta_cols(seurat_obj()), selected = "seurat_clusters")
  })
  
  output$comp_stack_ui <- renderUI({
    req(seurat_obj())
    choices <- c("(none)" = "", cat_meta_cols(seurat_obj()))
    selectInput("comp_stack", "Stack by (optional)", choices = choices, selected = "")
  })
  
  output$composition_plot <- renderPlotly({
    req(seurat_obj(), input$comp_group)
    md <- seurat_obj()@meta.data
    grp <- input$comp_group
    stk <- if (nzchar(input$comp_stack)) input$comp_stack else NULL
    df <- md %>%
      mutate(.g = .data[[grp]], .s = if (!is.null(stk)) .data[[stk]] else factor("all")) %>%
      count(.g, .s) %>%
      group_by(.g) %>%
      mutate(pct = n / sum(n)) %>%
      ungroup()
    
    fill_levels <- unique(df$.s)
    
    colors <- setNames(
      pals::glasbey(length(fill_levels)),
      fill_levels
    )
    
    p <- ggplot(df, aes(x = .g, y = pct, fill = .s)) +
      geom_bar(stat = "identity") +
      ylab("Proportion") + xlab(grp) +
      theme_minimal() +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(fill = ifelse(is.null(stk), "", stk)) +
      scale_fill_manual(values = colors)
    set_last_plot(p, "comp_raw")
    ggplotly(p)
  })
  
  output$download_comp_plot <- downloadHandler(
    filename = function() sprintf("composition_plot_%s.png", Sys.Date()),
    content = function(file) {
      req(exists("comp_raw", envir = last_plot_env))
      ggplot2::ggsave(file, plot = get_last_plot("comp_raw"), width = 7, height = 5, dpi = 300)
    }
  )
  
  # ---- Gene violins ----
  output$gene_violin_plot <- renderPlotly({
    req(seurat_obj(), nzchar(input$genes_violin), input$violin_group)
    obj <- seurat_obj()
    genes <- unlist(strsplit(input$genes_violin, split = ","))
    genes <- trimws(genes)
    genes <- genes[genes %in% rownames(obj)]
    validate(need(length(genes) > 0, "No valid genes found"))
    
    group_var <- input$violin_group
    validate(need(group_var %in% colnames(obj@meta.data), "Invalid group variable"))
    
    md <- obj@meta.data
    df <- FetchData(obj, vars = genes)
    df[[group_var]] <- md[[group_var]]
    
    df_long <- tidyr::pivot_longer(df, cols = all_of(genes), names_to = "gene", values_to = "expr")
    
    p <- ggplot(df_long, aes(x = .data[[group_var]], y = expr, fill = .data[[group_var]])) +
      geom_violin(scale = "width") +
      facet_wrap(~ gene, scales = "free_y") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      guides(fill = "none")
    
    set_last_plot(p, "gene_violin_raw")
    ggplotly(p)
  })
  
  # Summary table for "Gene Violins" tab (grouped by input$violin_group)
  output$gene_violin_summary <- renderDT({
    req(seurat_obj(), nzchar(input$genes_violin), input$violin_group)
    
    genes <- trimws(unlist(strsplit(input$genes_violin, ",")))
    genes <- genes[genes %in% rownames(seurat_obj())]
    validate(need(length(genes) > 0, "No valid genes entered"))
    
    group_var <- input$violin_group
    validate(need(group_var %in% colnames(seurat_obj()@meta.data), "Invalid group variable"))
    
    df <- FetchData(seurat_obj(), vars = c(genes, group_var))
    
    summary_df <- df %>%
      tidyr::pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expr") %>%
      group_by(.data[[group_var]], gene) %>%
      summarise(mean_expr = round(mean(expr, na.rm = TRUE), 3), .groups = "drop") %>%
      rename(group = !!group_var)
    
    datatable(summary_df, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # ---- Filtered Gene violins ----
  filtered_metadata <- reactive({
    req(seurat_obj(), input$filter_var, input$filter_value)
    md <- seurat_obj()@meta.data
    md[md[[input$filter_var]] == input$filter_value, , drop = FALSE]
  })
  
  output$filtered_violin_plot <- renderPlotly({
    req(seurat_obj(), input$filtered_genes, input$filtered_violin_group)
    
    genes <- trimws(unlist(strsplit(input$filtered_genes, ",")))
    genes <- genes[genes %in% rownames(seurat_obj())]
    validate(need(length(genes) > 0, "No valid genes entered"))
    
    cells_to_keep <- rownames(filtered_metadata())
    obj_filtered <- subset(seurat_obj(), cells = cells_to_keep)
    
    df <- FetchData(obj_filtered, vars = genes)
    df[[input$filtered_violin_group]] <- obj_filtered@meta.data[[input$filtered_violin_group]]
    
    df_long <- tidyr::pivot_longer(df, cols = all_of(genes), names_to = "gene", values_to = "expr")
    
    p <- ggplot(df_long, aes(x = .data[[input$filtered_violin_group]], y = expr, fill = .data[[input$filtered_violin_group]])) +
      geom_violin(scale = "width") +
      facet_wrap(~ gene, scales = "free_y") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      guides(fill = "none")
    
    set_last_plot(p, "filtered_violin_raw")
    ggplotly(p)
  })
  
  output$filtered_violin_summary <- renderDT({
    req(seurat_obj(), nzchar(input$filtered_genes), input$filtered_violin_group, input$filter_var, input$filter_value)
    
    genes <- trimws(unlist(strsplit(input$filtered_genes, ",")))
    genes <- genes[genes %in% rownames(seurat_obj())]
    validate(need(length(genes) > 0, "No valid genes entered"))
    
    group_var <- input$filtered_violin_group
    validate(need(group_var %in% colnames(seurat_obj()@meta.data), "Invalid group variable"))
    
    # Filtered cells only
    cells_to_keep <- rownames(filtered_metadata())
    obj_filtered <- subset(seurat_obj(), cells = cells_to_keep)
    
    df <- FetchData(obj_filtered, vars = c(genes, group_var))
    
    summary_df <- df %>%
      tidyr::pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expr") %>%
      group_by(.data[[group_var]], gene) %>%
      summarise(mean_expr = round(mean(expr, na.rm = TRUE), 3), .groups = "drop") %>%
      rename(group = !!group_var)
    
    datatable(summary_df, options = list(pageLength = 10, scrollX = TRUE))
  })
}

# ---------- Run App ----------
shinyApp(ui, server)
