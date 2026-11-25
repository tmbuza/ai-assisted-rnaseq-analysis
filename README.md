# AI-Assisted RNA-Seq Analysis with R

A hands-on Quarto book: fictional counts → Bioconductor airway experiment → your study. Each example follows Define → Preview → Meet → Check → Run → Verify and Explain → Write.

Open the folder containing `_quarto.yml`. Follow `00-setup.qmd` to initialize `renv` and install the R/Bioconductor packages; begin with `Rscript scripts/R/run_workflow.R data/synthetic results/synthetic --check-only`. Then run the synthetic and airway analyses as the chapters direct. The book itself has execution disabled and does not create analysis outputs during rendering.

**Provenance:** Airway data are loaded from the Bioconductor `airway` package, based on Himes et al. 2014, DOI 10.1371/journal.pone.0099625. Source: https://bioconductor.org/packages/airway/ . Statistical design guidance: https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html .

This guide starts from count matrices, not sequencing reads. Synthetic data are fictional and have no biological interpretation. For your study, review the design and contrast before running or interpreting results.

## Reproducible R environment

This starter has no tested `renv.lock`. In the main folder, initialize `renv` and install packages as described in `00-setup.qmd`. After the workflow runs successfully, call `renv::snapshot()` and retain `renv.lock`, `.Rprofile` and `renv/activate.R` in the repository. A learner receiving that validated lockfile can use `renv::restore()`. Do not commit `renv/library/`.
