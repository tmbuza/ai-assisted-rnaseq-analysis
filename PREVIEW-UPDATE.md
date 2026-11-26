# RNA-Seq preview update

Copy this package's contents into the main project folder, preserving paths. It updates the supplied lessons and scripts and includes matching synthetic inputs and results. It contains no CSS or book configuration, so your current styling remains in place.

The standard entry point remains:

```bash
Rscript scripts/R/run_workflow.R data/synthetic results/synthetic
Rscript scripts/R/prepare_airway.R
Rscript scripts/R/run_workflow.R data/airway results/airway
```

`02-analyse-example.R` is an equivalent entry point forwarding to the same workflow. It is not a second analysis implementation.

For guide maintainers, rebuild both prepared result sets before distributing:

```bash
Rscript scripts/R/build_preview_results.R
quarto render
```

The builder regenerates fictional paired counts (six donors, twelve samples), exports the real airway package and runs the same DESeq2 workflow for each. It stops on failure. Learners should start by viewing the bundled figures, then reproduce them later in the lesson.

The synthetic inputs replace the earlier independent-group example. Old A01 totals are no longer applicable. The own-study workflow still supports independent A/B samples without a donor column, or complete paired A/B samples with donor. Adapt the design deliberately for studies outside these examples.

Six figures: library sizes, PCA, sample distances, MA, volcano and selected-gene normalized counts. Figures have supporting TSVs. DE results are unshrunken; volcano uses adjusted p-values. Selected gene is chosen after fitting for illustration, not confirmation. Normalized counts and VST do not remove donor effects.

## Validation of this update

Executed with R 4.5.3, DESeq2 1.50.2, airway 1.30.0 and ggplot2 4.0.4.
Both workflows generated all six figures and supporting tables.
Independent checks passed for all raw sample summaries, all retained normalized counts, selected-gene values and sample distances.
A01: 10,456 counts and 120 detected genes. Airway SRR1039508: 20,637,971 counts and 24,633 detected genes.
An incomplete donor pair was rejected by input validation.
The complete 24-page Quarto book rendered successfully with Quarto 1.11.9.
The package excludes CSS, `_quarto.yml` and rendered `docs/`; render within your existing project to retain its styling and publishing configuration.
