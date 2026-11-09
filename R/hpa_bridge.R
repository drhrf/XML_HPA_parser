# Bridge between the Shiny app and the hpa_xml_parser Python package.
#
# All parsing, downloading and exporting is done by the Python package through
# reticulate, so the app and the command-line tool share one implementation and
# produce identical files.
#
# Python selection: set the environment variable HPA_XML_PARSER_PYTHON to the
# Python executable (or virtualenv) that has hpa-xml-parser installed, or let
# reticulate pick its default.

HPA_INSTALL_HINT <- paste(
  "Install the Python package in the Python used by reticulate, e.g.:",
  "reticulate::py_install(\"hpa-xml-parser @ git+https://github.com/drhrf/XML_HPA_parser\", pip = TRUE)"
)

#' Import the hpa_xml_parser Python module (cached).
hpa_module <- local({
  mod <- NULL
  function() {
    if (is.null(mod)) {
      py <- Sys.getenv("HPA_XML_PARSER_PYTHON", "")
      if (nzchar(py)) reticulate::use_python(py, required = TRUE)
      if (!reticulate::py_module_available("hpa_xml_parser")) {
        stop("Python package 'hpa_xml_parser' not found. ", HPA_INSTALL_HINT, call. = FALSE)
      }
      mod <<- reticulate::import("hpa_xml_parser", convert = FALSE)
    }
    mod
  }
})

#' Parse local HPA XML files (.xml or .xml.gz). Returns a Python list of ParsedGene.
hpa_parse_files <- function(paths, retrieved_at = NULL) {
  hpa <- hpa_module()
  genes <- lapply(paths, function(p) hpa$parse_xml(p, retrieved_at = retrieved_at))
  reticulate::r_to_py(genes)
}

#' Split free text ("ACADM, GFAP ENSG00000075624") into gene identifiers.
hpa_split_ids <- function(text) {
  ids <- unlist(strsplit(text %||% "", "[[:space:],;]+"))
  unique(ids[nzchar(ids)])
}

#' Download genes by Ensembl ID or symbol. Returns a Python list of ParsedGene.
hpa_fetch <- function(ids, ...) {
  hpa_module()$fetch_genes(as.list(ids), on_error = "warn", ...)
}

# Convert a pandas DataFrame to an R data.frame. Columns are handed over as
# Python lists with missing values (None, NaN, pd.NA) mapped to None, so the
# result does not depend on how a reticulate version treats pandas' nullable
# dtypes. Numeric pandas columns become numeric R columns; all others character.
.pandas_to_r <- function(df) {
  to_lists <- reticulate::py_run_string(paste(
    "import pandas as pd",
    "def _hpa_to_lists(d):",
    "    return {c: (pd.api.types.is_numeric_dtype(d[c]),",
    "                [None if pd.isna(v) else v for v in d[c].astype(object).tolist()])",
    "            for c in d.columns}",
    sep = "\n"
  ), convert = FALSE)$`_hpa_to_lists`
  cols <- reticulate::py_to_r(to_lists(df))
  if (length(cols) == 0) return(data.frame())
  out <- lapply(cols, function(col) {
    is_num <- isTRUE(col[[1]])
    vals <- col[[2]]
    if (is_num) {
      vapply(vals, function(v) if (is.null(v)) NA_real_ else as.numeric(v), 0)
    } else {
      vapply(vals, function(v) if (is.null(v)) NA_character_ else as.character(v), "")
    }
  })
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

#' Per-sample long table (same columns as the Python package) as an R data.frame.
hpa_samples <- function(genes, assays = NULL) {
  assays_py <- if (is.null(assays)) NULL else as.list(assays)
  .pandas_to_r(hpa_module()$to_frame(genes, assays = assays_py, kind = "samples"))
}

#' Aggregated level values as an R data.frame.
hpa_aggregates <- function(genes, assays = NULL) {
  assays_py <- if (is.null(assays)) NULL else as.list(assays)
  .pandas_to_r(hpa_module()$to_frame(genes, assays = assays_py, kind = "aggregates"))
}

#' Provenance table (HPA/Ensembl versions, retrieval time) as an R data.frame.
hpa_metadata <- function(genes) {
  .pandas_to_r(hpa_module()$metadata_frame(genes))
}

#' Write files with the Python exporter. Returns the written paths.
hpa_export <- function(genes, out, fmt = NULL, assays = NULL, include_aggregates = TRUE) {
  assays_py <- if (is.null(assays)) NULL else as.list(assays)
  paths <- hpa_module()$export(genes, out, fmt = fmt, assays = assays_py,
                               include_aggregates = include_aggregates)
  vapply(reticulate::py_to_r(paths), function(p) as.character(reticulate::py_str(p)), "")
}

#' Assays with per-sample data, as defined by the Python package.
hpa_sample_assays <- function() {
  unlist(reticulate::py_to_r(hpa_module()$SAMPLE_ASSAYS))
}

#' Bin ages into groups for plotting. Breaks follow decades; NA stays NA.
hpa_age_group <- function(age, breaks = c(-Inf, 1, 20, 40, 60, 70, 80, 90, Inf)) {
  labels <- c("<1", "1-19", "20-39", "40-59", "60-69", "70-79", "80-89", "90+")
  cut(age, breaks = breaks, labels = labels, right = FALSE)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
