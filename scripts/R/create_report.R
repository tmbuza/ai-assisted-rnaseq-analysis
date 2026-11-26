#!/usr/bin/env Rscript
# Base-R helper; preserve existing learner writing.
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 1L || !a[1] %in% c("synthetic", "airway", "my-study"))
  stop("Usage: Rscript scripts/R/create_report.R synthetic|airway|my-study", call. = FALSE)
if (!file.exists("_quarto.yml") || !file.exists("templates/report.qmd"))
  stop("Run from the main project folder containing _quarto.yml and templates/report.qmd.", call. = FALSE)
dataset <- a[1]
labels <- c(synthetic = "Synthetic RNA-seq practice", airway = "Airway RNA-seq case study", "my-study" = "My RNA-seq study")
dir.create("reports", showWarnings = FALSE)
if (!file.exists("reports/_quarto.yml"))
  writeLines(c("project:", "  type: default", "  output-dir: ."), "reports/_quarto.yml")
path <- file.path("reports", paste0(dataset, "-report.qmd"))
if (file.exists(path)) {
  cat("Existing report kept:", path, "\n")
} else {
  s <- readLines("templates/report.qmd", warn = FALSE)
  s <- gsub("__DATASET__", dataset, s, fixed = TRUE)
  s <- gsub("__LABEL__", labels[[dataset]], s, fixed = TRUE)
  writeLines(s, path)
  cat("Created", path, "\n")
}
cat("Edit the report; enter reports/ to render it. Return with cd .. afterward.\n")
