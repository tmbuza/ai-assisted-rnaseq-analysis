#!/usr/bin/env Rscript
if (!requireNamespace("airway", quietly=TRUE)) stop("Install airway with BiocManager::install('airway')")
if (!requireNamespace("SummarizedExperiment", quietly=TRUE)) stop("Install SummarizedExperiment")
data("airway", package="airway")
se <- airway
counts <- SummarizedExperiment::assay(se)
meta <- as.data.frame(SummarizedExperiment::colData(se))
stopifnot(identical(colnames(counts), rownames(meta)))
dir.create("data/airway", recursive=TRUE, showWarnings=FALSE)
write.table(data.frame(gene_id=rownames(counts), counts, check.names=FALSE), "data/airway/gene-counts.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(data.frame(sample_id=rownames(meta), cell=meta$cell, dex=meta$dex), "data/airway/sample-metadata.tsv", sep="\t", quote=FALSE, row.names=FALSE)
message("Exported airway counts and metadata for ", ncol(counts), " samples")
