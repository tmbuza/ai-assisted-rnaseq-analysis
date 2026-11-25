#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3) stop("Usage: Rscript scripts/R/verify_sample.R data_dir results_dir sample_id")
x <- read.delim(file.path(args[1],"gene-counts.tsv"),check.names=FALSE)
qc <- read.delim(file.path(args[2],"tables/sample-qc.tsv"))
id <- args[3]
if(!id %in% names(x) || !id %in% qc$sample_id) stop("Sample missing from input or output")
v <- x[[id]]; saved <- qc[match(id,qc$sample_id),]
stopifnot(all(!is.na(v)),all(v>=0),sum(v)==saved$total_counts,sum(v>0)==saved$detected_genes)
cat("sample_id:",id,"\ntotal_counts:",sum(v),"\ndetected_genes:",sum(v>0),"\nSample verification passed.\n")
