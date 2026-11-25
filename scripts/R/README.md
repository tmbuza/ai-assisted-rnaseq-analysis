# R scripts

- `prepare_synthetic.R`: regenerate the fictional paired dataset (seed 20261125).
- `prepare_airway.R`: export original airway counts and paired metadata.
- `run_workflow.R`: validate inputs, fit DESeq2 and save six figures with supporting tables.
- `02-analyse-example.R`: alternate entry point forwarding to `run_workflow.R`.
- `verify_sample.R`: independently check a sample total and detected genes.
- `build_preview_results.R`: maintainer command regenerating both datasets and all bundled previews.

Run commands from the main project folder. Learners view bundled results first, then reproduce them using `run_workflow.R`.
