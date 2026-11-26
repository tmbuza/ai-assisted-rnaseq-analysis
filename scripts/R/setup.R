#!/usr/bin/env Rscript
# Restore the supplied, maintainer-tested environment; do not snapshot here.
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!file.exists("_quarto.yml") || !dir.exists("scripts/R"))
  stop("Run setup from the main ai-assisted-rnaseq-analysis folder.", call. = FALSE)
if (!file.exists("renv.lock"))
  stop("renv.lock is missing. Obtain the complete project with its tested lockfile before learner setup. No packages were changed.", call. = FALSE)
if (!file.exists(".Rprofile") || !file.exists("renv/activate.R"))
  stop("renv activation files are missing. Obtain .Rprofile and renv/activate.R with the project.", call. = FALSE)
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::load(project = getwd())
tryCatch(renv::restore(prompt = FALSE), error = function(e)
  stop("Environment restore failed: ", conditionMessage(e),
       "\nStop here and inspect the first installation error. Do not snapshot a failed installation.", call. = FALSE))
required <- c("DESeq2", "ggplot2", "airway", "SummarizedExperiment")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Required packages cannot load: ", paste(missing, collapse = ", "), call. = FALSE)
s <- renv::status()
if (!isTRUE(s$synchronized)) stop("Project is not synchronized. Inspect renv::status() before continuing.", call. = FALSE)
cat("Setup complete. Continue with synthetic input checking.\n")
