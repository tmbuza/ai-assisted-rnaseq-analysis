#!/usr/bin/env Rscript
# Maintainer step: run once before rendering or distributing the guide.
# Learners see bundled outputs first and reproduce them later.
rscript <- file.path(R.home("bin"),"Rscript")
run <- function(args) {
 status <- system2(rscript,args)
 if(status!=0) stop("Preview preparation failed: ",paste(args,collapse=" "))
}
run("scripts/R/prepare_synthetic.R")
run(c("scripts/R/run_workflow.R","data/synthetic","results/synthetic"))
run("scripts/R/prepare_airway.R")
run(c("scripts/R/run_workflow.R","data/airway","results/airway"))
message("Both sets of preview results are ready. Render the guide with quarto render.")
