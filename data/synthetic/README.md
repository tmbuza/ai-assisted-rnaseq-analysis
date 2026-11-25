# Fictional paired RNA-seq experiment, version 2

Regenerate with `Rscript scripts/R/prepare_synthetic.R` from the main project folder.
Seed 20261126; 2,000 fictional genes; six donors measured once in A and once in B.
Model: ~ donor + group; contrast: B versus A. Data are not from a real organism.

Baseline expression spans a broad range; gene-wise NB dispersion varies. Most condition effects are zero, and nonzero effects vary in size and direction. Donor effects persist across pairs, while treatment responses also vary among donors. Both conditions have the same set of generating library multipliers in a different order. Twenty genes have zero generating mean to illustrate filtering.

`simulation-truth.tsv` records baseline means, dispersions, central log2 effects and donor-specific effects. These are generating parameters, not fitted DESeq2 estimates or a list of statistically significant genes. `simulation-samples.tsv` records library multipliers, which differ from DESeq2 estimated size factors. `simulation-record.txt` records the version and seed. The rules and seed were fixed before model fitting; there is no target number of discoveries.

`plot-examples.tsv` prespecifies a null-effect gene and a variable-response gene for the three-panel expression display. Neither is chosen by its fitted p-value. A third panel shows a gene selected after analysis for illustration.

This replaces the earlier 120-gene simulation. Reproduce both inputs and results together; previous sample totals and gene-level values no longer apply.

## Regenerate this example

From the main project folder:

```bash
Rscript scripts/R/rebuild_synthetic_preview.R
quarto render
```

The rebuild script prepares synthetic inputs, runs the common workflow and verifies A01. It does not prepare or analyse airway. All prepared synthetic figures are already bundled for previewing before practice.

## Checked results

This version was executed with R 4.5.3 and DESeq2 1.50.2. The run retained 1,839 genes from 2,000 input rows; A01 has 629,417 counts and 1,935 detected genes. These values replace the earlier 120-gene simulation's results. Count summaries, size-factor normalization, distances, PCA and all three gene-expression panels passed independent numerical checks. Fixed-seed input regeneration was byte-for-byte reproducible. The 24-page guide rendered with Quarto 1.11.9.

The prepared update contains only synthetic inputs/results, the common workflow script, the revised simulation and affected lessons/setup/report template. It excludes airway inputs/results, CSS and book configuration. Copy its folders into the existing project while preserving paths. Re-render the guide to update your published pages. Any personal report edits should be merged with the supplied template.
