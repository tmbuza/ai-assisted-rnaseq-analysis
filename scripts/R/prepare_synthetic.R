#!/usr/bin/env Rscript
# Fictional paired RNA-seq experiment, version 2. Rules fixed before analysis.
# Run from the main project folder. This overwrites synthetic inputs only.
set.seed(20261126)
folder <- "data/synthetic"
dir.create(folder,recursive=TRUE,showWarnings=FALSE)
n_genes <- 2000L; n_donors <- 6L
genes <- sprintf("Gene_%04d",seq_len(n_genes))
ids <- c(sprintf("A%02d",1:6),sprintf("B%02d",1:6))
meta <- data.frame(sample_id=ids,group=rep(c("A","B"),each=6),donor=rep(sprintf("D%02d",1:6),2))
# Broad expression range, with explicitly unexpressed genes for filter practice.
baseline <- exp(rnorm(n_genes,4.2,1.8))
baseline[1:3] <- c(800,1200,600)
baseline[1981:2000] <- 0
# Most genes have no simulated condition effect. Effects vary continuously.
log2_effect <- numeric(n_genes)
changed <- sample(4:1980,300,replace=FALSE)
log2_effect[changed[1:150]] <- runif(150,.35,2.3)
log2_effect[changed[151:300]] <- -runif(150,.35,2.3)
log2_effect[1] <- 1.8; log2_effect[3] <- .7
# Persistent, gene-specific donor effects occur in both paired samples.
donor_effect <- matrix(rnorm(n_genes*n_donors,0,.5),n_genes,n_donors)
# Treatment response varies among donors, including a prespecified example.
response_deviation <- matrix(0,n_genes,n_donors)
nonzero <- which(log2_effect!=0)
response_deviation[nonzero,] <- matrix(rnorm(length(nonzero)*n_donors,0,.3),length(nonzero),n_donors)
response_deviation[1,] <- c(-.2,.15,-.1,.25,-.15,.05)
response_deviation[3,] <- c(-1.2,.8,-.6,.9,-.7,.8)
# Equal sets of library multipliers in A and B; not estimated size factors.
library_multiplier <- c(.7,1.05,1.35,.85,1.2,.9,1.2,.85,.7,1.35,.9,1.05)
# Gene-wise NB dispersion: more relative noise at lower expression, with variation.
dispersion <- (.035 + .65 / sqrt(pmax(baseline,1))) * exp(rnorm(n_genes,0,.35))
mu <- sapply(1:12,function(j) {
 donor <- ((j-1) %% n_donors)+1
 response <- if(j>6) log2_effect+response_deviation[,donor] else 0
 baseline * exp(donor_effect[,donor]) * 2^response * library_multiplier[j]
})
cnt <- matrix(rnbinom(length(mu),mu=as.vector(mu),size=rep(1/dispersion,12)),n_genes,12,dimnames=list(genes,ids))
write.table(data.frame(gene_id=genes,cnt,check.names=FALSE),file.path(folder,"gene-counts.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(meta,file.path(folder,"sample-metadata.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
truth <- data.frame(gene_id=genes,baseline_mean=baseline,simulated_log2_effect=log2_effect,dispersion=dispersion)
for(k in 1:6) truth[[paste0("log2_effect_D",sprintf("%02d",k))]] <- log2_effect+response_deviation[,k]
write.table(truth,file.path(folder,"simulation-truth.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(data.frame(gene_id=genes[c(2,3)],example=c("No simulated condition effect","Variable donor response")),file.path(folder,"plot-examples.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(data.frame(sample_id=ids,library_multiplier=library_multiplier),file.path(folder,"simulation-samples.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(c("Simulation version: 2","Seed: 20261126","2000 genes; 6 paired donors; 12 samples","20 genes have zero generating mean","Rules and seed fixed before fitting DESeq2; no target discovery count","Truth tables describe generating parameters, not fitted estimates","Library multipliers are not DESeq2 estimated size factors"),file.path(folder,"simulation-record.txt"))
message("Prepared synthetic v2: 2000 genes, 6 donors, 12 samples")
