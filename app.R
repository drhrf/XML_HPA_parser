# HPA XML Parser - Shiny app
#
# Graphical front end for the hpa_xml_parser Python package. Parsing,
# downloading and exporting are delegated to Python through reticulate
# (see R/hpa_bridge.R), so the app and the command-line tool produce the same
# columns and the same files.
#
# Run with: shiny::runApp()   (from the repository root)
# Requires: R packages shiny, bslib, reticulate, DT, ggplot2 and the Python
# package hpa-xml-parser in the Python used by reticulate.

library(shiny)
library(bslib)
library(DT)
library(ggplot2)

source(file.path("R", "hpa_bridge.R"), local = TRUE)

options(shiny.maxRequestSize = 50 * 1024^2)

ui <- page_navbar(
  title = "HPA XML Parser",
  theme = bs_theme(version = 5),
  sidebar = sidebar(
    width = 320,
    radioButtons("input_mode", "Input",
                 c("Gene ID or symbol" = "fetch", "Upload XML file(s)" = "upload")),
    conditionalPanel(
      "input.input_mode == 'fetch'",
      textAreaInput("gene_ids", "Ensembl IDs or gene symbols",
                    placeholder = "ACADM, GFAP, ENSG00000075624", rows = 3),
      helpText("Downloaded from www.proteinatlas.org and cached locally.")
    ),
    conditionalPanel(
      "input.input_mode == 'upload'",
      fileInput("xml_files", "HPA per-gene XML (.xml or .xml.gz)",
                multiple = TRUE, accept = c(".xml", ".gz"))
    ),
    actionButton("run", "Load data", class = "btn-primary"),
    hr(),
    uiOutput("version_box"),
    hr(),
    checkboxGroupInput("assays", "Assays to show and export", choices = NULL),
    radioButtons("export_format", "Export format",
                 c("Excel (.xlsx)" = "xlsx", "CSV (.zip)" = "csv", "TSV (.zip)" = "tsv")),
    downloadButton("download", "Download")
  ),
  nav_panel(
    "Samples",
    card(DTOutput("samples_table"))
  ),
  nav_panel(
    "Aggregates",
    card(
      card_header("Values aggregated by HPA per tissue / brain region (nTPM, pTPM, TPM)"),
      DTOutput("aggregates_table")
    )
  ),
  nav_panel(
    "Plots",
    layout_columns(
      col_widths = c(3, 3, 3, 3),
      selectInput("plot_gene", "Gene", choices = NULL),
      selectInput("plot_assay", "Assay", choices = NULL),
      selectInput("plot_group", "Group by",
                  c("Region group / organ" = "organ", "Region / tissue / subregion" = "region",
                    "Sex" = "sex", "Age group" = "age_group")),
      checkboxInput("plot_log", "log10(nTPM + 1) scale", value = TRUE)
    ),
    card(full_screen = TRUE, plotOutput("boxplot", height = "auto")),
    helpText(
      "Each point is one sample. Age groups: <1, 1-19, 20-39, 40-59, 60-69, 70-79, 80-89, 90+ ",
      "(years as reported by HPA). Brain samples are post mortem. HPA gives no donor ID, so ",
      "samples from the same donor cannot be identified: these plots are descriptive and ",
      "points should not be treated as independent observations."
    )
  ),
  nav_panel(
    "Metadata",
    card(
      card_header("Provenance of each gene (written to the 'metadata' sheet on export)"),
      tableOutput("metadata_table")
    )
  )
)

server <- function(input, output, session) {
  state <- reactiveValues(genes = NULL, samples = NULL, aggregates = NULL, metadata = NULL)

  observeEvent(input$run, {
    genes <- tryCatch(
      withProgress(message = "Loading HPA data...", {
        if (input$input_mode == "fetch") {
          ids <- hpa_split_ids(input$gene_ids)
          validate(need(length(ids) > 0, "Enter at least one gene."))
          hpa_fetch(ids)
        } else {
          req(input$xml_files)
          hpa_parse_files(input$xml_files$datapath)
        }
      }),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = NULL)
        NULL
      }
    )
    req(genes)
    if (reticulate::py_len(genes) == 0) {
      showNotification("No gene could be loaded.", type = "error")
      return()
    }
    state$genes <- genes
    state$samples <- hpa_samples(genes)
    state$aggregates <- hpa_aggregates(genes)
    state$metadata <- hpa_metadata(genes)

    available <- intersect(hpa_sample_assays(), unique(state$samples$assay_type))
    updateCheckboxGroupInput(session, "assays", choices = available, selected = available)
    updateSelectInput(session, "plot_gene", choices = unique(state$samples$gene_name))
    updateSelectInput(session, "plot_assay", choices = available,
                      selected = if ("humanBrain" %in% available) "humanBrain" else available[1])
  })

  selected_samples <- reactive({
    req(state$samples)
    state$samples[state$samples$assay_type %in% input$assays, , drop = FALSE]
  })

  output$version_box <- renderUI({
    meta <- state$metadata
    if (is.null(meta)) return(helpText("No data loaded."))
    versions <- unique(meta$hpa_version)
    tagList(
      strong("HPA release: ", paste(versions, collapse = ", ")),
      if (length(versions) > 1) div(class = "text-danger", "Genes come from different HPA releases!"),
      br(),
      sprintf("Ensembl %s, %s", paste(unique(meta$ensembl_version), collapse = ", "),
              paste(unique(meta$assembly), collapse = ", ")),
      br(),
      sprintf("%d gene(s), %d samples", nrow(meta), nrow(state$samples))
    )
  })

  output$samples_table <- renderDT({
    datatable(selected_samples(), rownames = FALSE,
              options = list(pageLength = 25, scrollX = TRUE))
  })

  output$aggregates_table <- renderDT({
    req(state$aggregates)
    datatable(state$aggregates, rownames = FALSE,
              options = list(pageLength = 25, scrollX = TRUE))
  })

  output$metadata_table <- renderTable(state$metadata)

  plot_data <- reactive({
    req(state$samples, input$plot_gene, input$plot_assay)
    df <- state$samples[state$samples$gene_name == input$plot_gene &
                          state$samples$assay_type == input$plot_assay, , drop = FALSE]
    df$age_group <- hpa_age_group(df$age)
    df
  })

  plot_height <- function() {
    df <- tryCatch(plot_data(), error = function(e) NULL)
    n <- if (is.null(df) || nrow(df) == 0) 1 else length(unique(df[[input$plot_group]]))
    max(400, 60 + 16 * n)
  }

  output$boxplot <- renderPlot(height = plot_height, {
    df <- plot_data()
    req(nrow(df) > 0)
    group <- input$plot_group
    df$group <- df[[group]]
    if (group %in% c("region", "organ")) {
      df$group <- stats::reorder(factor(df$group), df$value, FUN = stats::median)
    }
    y_lab <- "nTPM"
    if (isTRUE(input$plot_log)) {
      df$value <- log10(df$value + 1)
      y_lab <- "log10(nTPM + 1)"
    }
    ggplot(df, aes(x = group, y = value)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
      coord_flip() +
      labs(x = NULL, y = y_lab,
           title = sprintf("%s - %s (HPA v%s)", input$plot_gene, input$plot_assay,
                           paste(unique(df$hpa_version), collapse = ", "))) +
      theme_minimal(base_size = 12)
  })

  output$download <- downloadHandler(
    filename = function() {
      stem <- paste0("hpa_", paste(unique(state$samples$gene_name), collapse = "_"))
      if (input$export_format == "xlsx") paste0(stem, ".xlsx") else paste0(stem, ".zip")
    },
    content = function(file) {
      req(state$genes)
      assays <- input$assays
      if (length(assays) == 0) assays <- NULL
      if (input$export_format == "xlsx") {
        tmp <- tempfile(fileext = ".xlsx")
        hpa_export(state$genes, tmp, fmt = "xlsx", assays = assays)
        file.copy(tmp, file, overwrite = TRUE)
      } else {
        dir <- tempfile("hpa_export_")
        dir.create(dir)
        target <- file.path(dir, paste0("hpa_samples.", input$export_format))
        paths <- hpa_export(state$genes, target, fmt = input$export_format, assays = assays)
        utils::zip(file, paths, flags = "-j -q")
      }
    }
  )
}

shinyApp(ui, server)
