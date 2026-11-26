# R scripts

- `prepare_synthetic.R`: regenerate the fictional paired dataset (seed 20261126).
- `prepare_airway.R`: export original airway counts and paired metadata.
- `run_workflow.R`: validate inputs, fit DESeq2 and save six figures with supporting tables.
- `02-analyse-example.R`: alternate entry point forwarding to `run_workflow.R`.
- `verify_sample.R`: independently check a sample total and detected genes.
- `build_preview_results.R`: maintainer command regenerating both datasets and all bundled previews.

Run commands from the main project folder. Learners view bundled results first, then reproduce them using `run_workflow.R`.

- `setup.R`: restore the supplied tested lockfile; no learner snapshot or package-version upgrades.
- `create_report.R`: create an ignored report starter and its separate Quarto project, preserving existing writing.

The fitted workflow and sample verifier retain their existing interfaces. The verifier checks raw totals and positive-count gene counts only; input/script hashes require a separate review. Do not use maintainer preview rebuilding as a routine learner step.
