#!/usr/bin/env Rscript
# Regenerate only the synthetic v2 inputs and outputs, then verify one sample.
rscript <- file.path(R.home("bin"),"Rscript")
run <- function(args) {
 status <- system2(rscript,args)
 if(status!=0) stop("Synthetic preview rebuild failed: ",paste(args,collapse=" "))
}
run("scripts/R/prepare_synthetic.R")
run(c("scripts/R/run_workflow.R","data/synthetic","results/synthetic"))
run(c("scripts/R/verify_sample.R","data/synthetic","results/synthetic","A01"))
