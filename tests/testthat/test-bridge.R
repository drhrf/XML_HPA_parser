# testthat tests for the Shiny app's bridge to the Python package.
# Run from the repository root:
#   Rscript -e 'testthat::test_dir("tests/testthat")'

root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
source(file.path(root, "R", "hpa_bridge.R"))
fixture <- file.path(root, "tests", "fixtures", "ENSG00000117054_ACADM_v25.xml.gz")

skip_without_python <- function() {
  py <- Sys.getenv("HPA_XML_PARSER_PYTHON", "")
  if (nzchar(py)) reticulate::use_python(py, required = TRUE)
  if (!reticulate::py_module_available("hpa_xml_parser")) {
    testthat::skip("Python package hpa_xml_parser not available")
  }
}

test_that("gene identifiers are split on commas, spaces and new lines", {
  expect_equal(hpa_split_ids("ACADM, GFAP\nENSG00000075624;  ACADM"),
               c("ACADM", "GFAP", "ENSG00000075624"))
  expect_equal(hpa_split_ids(""), character(0))
})

test_that("age groups", {
  expect_equal(as.character(hpa_age_group(c(0.2, 1, 18, 42, 61, 75, 85, 94, NA))),
               c("<1", "1-19", "1-19", "40-59", "60-69", "70-79", "80-89", "90+", NA))
})

test_that("R receives the same table as Python (columns, rows, types)", {
  skip_without_python()
  genes <- hpa_parse_files(fixture, retrieved_at = "2026-11-09T16:23:00Z")
  s <- hpa_samples(genes)
  core <- unlist(reticulate::py_to_r(hpa_module()$SAMPLE_COLUMNS))
  expect_equal(names(s)[seq_along(core)], core)
  expect_equal(nrow(s), 1447)
  expect_equal(as.list(table(s$assay_type)),
               list(humanBrain = 966L, immuneCell = 109L, mouseBrain = 67L,
                    pigBrain = 119L, tissue = 186L))
  expect_type(s$value, "double")
  expect_type(s$age, "double")
  expect_type(s$sample_id, "character")
  expect_true(all(s$hpa_version == "25"))

  amy <- s[s$assay_type == "humanBrain" & s$sample_id == "1987", ]
  expect_equal(amy$region, "Amygdala")
  expect_equal(amy$organ, "AMY")
  expect_equal(amy$value, 2.4)
  expect_equal(amy$sex, "Female")
  expect_equal(amy$age, 91)

  bm <- s[s$assay_type == "tissue" & s$sample_id == "97", ]
  expect_equal(bm$sex, "Female")
  expect_true(is.na(bm$age))
  imm <- s[s$assay_type == "immuneCell", ]
  expect_true(all(is.na(imm$age)))
})

test_that("concordance holds on the R side too", {
  skip_without_python()
  genes <- hpa_parse_files(fixture)
  s <- hpa_samples(genes)
  a <- hpa_aggregates(genes)
  a <- a[a$level_type == "normalizedRNAExpression" & a$n_samples > 0, ]
  m <- aggregate(value ~ assay_type + organ + region, data = s, FUN = mean)
  j <- merge(m, a, by = c("assay_type", "organ", "region"), suffixes = c("_mean", "_hpa"))
  expect_equal(nrow(j), nrow(a))
  expect_lte(max(abs(j$value_mean - j$value_hpa)), 0.1)
})

test_that("metadata and export match the Python exporter", {
  skip_without_python()
  genes <- hpa_parse_files(fixture)
  meta <- hpa_metadata(genes)
  expect_equal(meta$hpa_version, "25")
  expect_equal(meta$ensembl_version, "109")

  xlsx <- tempfile(fileext = ".xlsx")
  expect_equal(hpa_export(genes, xlsx), xlsx)
  expect_equal(openxlsx::getSheetNames(xlsx), c(hpa_sample_assays(), "aggregates", "metadata"))

  csv <- tempfile(fileext = ".csv")
  paths <- hpa_export(genes, csv, assays = "humanBrain")
  expect_length(paths, 3)
  back <- utils::read.csv(csv, colClasses = c(sample_id = "character"))
  expect_equal(nrow(back), 966)
})
