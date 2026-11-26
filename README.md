# AI-Assisted RNA-Seq Analysis with R

A Quarto guide from gene-level count data to results you can explain and defend: **Define → Preview → Meet → Check → Run → Verify and Explain → Write**.

Begin with a fictional paired experiment, repeat the process with the real Bioconductor airway data, then adapt the supported workflow to your study. This practical starts from count matrices, not FASTQ reads.

## Start here

Download the complete repository using Code → Download ZIP, extract it and open the inner project folder. Follow [00-setup.qmd](00-setup.qmd). Learners restore a supplied, tested renv.lock and activation files using:

```bash
Rscript scripts/R/setup.R
Rscript scripts/R/run_workflow.R data/synthetic results/synthetic --check-only
```

If the project environment has not been supplied, the setup helper stops clearly. It does not create an untested lockfile for a learner. Maintainer-only environment creation is explained in the appendix.

## The worked paths

- Synthetic v2: 2,000 genes, six paired donors and 12 samples; seed 20261126; formula `~ donor + group`, B relative to A.
- Airway: eight samples from four paired cell lines; formula `~ cell + dex`, trt relative to untrt. Source: airway package and Himes et al. (2014), DOI 10.1371/journal.pone.0099625.
- My study: review the supported design interface, source and contrast before inference. The script does not automatically fit arbitrary covariates or studies.

## Outputs and reporting

Six figures and their supporting TSVs are written to `results/DATASET/figures/` and `results/DATASET/tables/`. Run logs and session information identify the analysis. See [01-run-the-workflow.qmd](01-run-the-workflow.qmd) for the output map.

Use `Rscript scripts/R/create_report.R synthetic` to create a report starter; use airway or my-study for other paths. The helper preserves existing reports and isolates rendering inside ignored reports/. Templates are tracked under templates/. See [02-write-the-report.qmd](02-write-the-report.qmd).

Quarto rendering displays saved results; it does not fit DESeq2. Keep scripts, inputs, environment and run provenance aligned. A successful raw sample verification does not validate the full statistical model.
