#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2 || length(args) > 3) stop("Usage: Rscript scripts/R/run_workflow.R data_dir results_dir [--check-only]")
check_only <- length(args)==3 && identical(args[3], "--check-only")
if (length(args)==3 && !check_only) stop("Third argument must be --check-only")
data_dir <- args[1]; output_dir <- args[2]
read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  read.delim(path, check.names=FALSE, stringsAsFactors=FALSE)
}
x <- read_tsv(file.path(data_dir,"gene-counts.tsv"))
m <- read_tsv(file.path(data_dir,"sample-metadata.tsv"))
if (!identical(names(x)[1],"gene_id") || !identical(names(m)[1],"sample_id")) stop("Expected gene_id and sample_id as first columns")
if (anyDuplicated(names(x)) || anyDuplicated(x$gene_id) || anyDuplicated(m$sample_id) || anyNA(x) || anyNA(m)) stop("Duplicate identifiers or missing values")
if (ncol(x)<3 || nrow(x)<2 || nrow(m)<2) stop("Too few genes or samples")
if (!setequal(names(x)[-1],m$sample_id)) stop("Count and metadata sample IDs differ")
m <- m[match(names(x)[-1], m$sample_id), , drop=FALSE]
if (!identical(m$sample_id,names(x)[-1])) stop("Sample alignment failed")
cnt <- as.matrix(x[-1]); suppressWarnings(storage.mode(cnt)<-"numeric")
if (anyNA(cnt) || any(!is.finite(cnt)) || any(cnt<0) || any(cnt!=floor(cnt))) stop("Counts must be nonnegative finite integers")
rownames(cnt) <- x$gene_id
mode <- if ("dex" %in% names(m) && "cell" %in% names(m)) "airway" else if ("group" %in% names(m)) "synthetic" else stop("Need group or both cell and dex columns")
if (mode=="synthetic") {
 if (!all(c("A","B") %in% m$group) || any(!m$group %in% c("A","B")) || min(table(m$group))<2) stop("Synthetic example needs replicated A and B groups")
 m$group <- factor(m$group,levels=c("A","B")); design <- ~ group
} else {
 if (!all(c("untrt","trt") %in% m$dex) || any(!m$dex %in% c("untrt","trt"))) stop("Expected airway dex values untrt and trt")
 if (any(table(m$cell,m$dex)!=1)) stop("Each airway cell line must have one sample per condition")
 m$cell <- factor(m$cell); m$dex <- factor(m$dex,levels=c("untrt","trt")); design <- ~ cell + dex
}
message("Input checks passed: ", ncol(cnt)," samples; ",nrow(cnt)," genes; design ",deparse(design))
if (check_only) quit(save="no")
if (!requireNamespace("DESeq2",quietly=TRUE)) stop("Install DESeq2 with BiocManager::install('DESeq2')")
if (!requireNamespace("ggplot2",quietly=TRUE)) stop("Install ggplot2")
dir.create(file.path(output_dir,"tables"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(output_dir,"figures"), recursive=TRUE, showWarnings=FALSE)
qc <- data.frame(sample_id=colnames(cnt), total_counts=colSums(cnt), detected_genes=colSums(cnt>0))
write.table(qc,file.path(output_dir,"tables/sample-qc.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
dds <- DESeq2::DESeqDataSetFromMatrix(countData=cnt,colData=m,design=design)
keep <- rowSums(DESeq2::counts(dds)>=10)>=2
if (sum(keep)<2) stop("Too few genes pass the documented filter")
dds <- dds[keep,]; dds <- DESeq2::DESeq(dds,quiet=TRUE)
contrast <- if (mode=="synthetic") c("group","B","A") else c("dex","trt","untrt")
res <- DESeq2::results(dds,contrast=contrast,alpha=.05)
tab <- data.frame(gene_id=rownames(res),as.data.frame(res),check.names=FALSE)
tab <- tab[order(tab$padj,na.last=TRUE),]
write.table(tab,file.path(output_dir,"tables/differential-expression.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
# vsd <- DESeq2::vst(dds,blind=FALSE)
vsd <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)
pca <- DESeq2::plotPCA(vsd,intgroup=if(mode=="airway") c("dex","cell") else "group",returnData=TRUE)
variance <- round(100*attr(pca,"percentVar"))
pca$label <- if(mode=="airway") pca$dex else pca$group
p <- ggplot2::ggplot(pca,ggplot2::aes(x=PC1,y=PC2,color=label,label=name))+ggplot2::geom_point(size=3)+ggplot2::geom_text(vjust=-.7,show.legend=FALSE)+ggplot2::labs(x=paste0("PC1: ",variance[1],"%"),y=paste0("PC2: ",variance[2],"%"),color=if(mode=="airway") "Treatment" else "Group")+ggplot2::theme_minimal()
ggplot2::ggsave(file.path(output_dir,"figures/sample-pca.png"),p,width=8,height=5,dpi=160)
writeLines(c(paste("input",normalizePath(data_dir)),paste("design",deparse(design)),paste("contrast",paste(contrast,collapse=" ")),paste("genes_tested",nrow(tab)),paste("padj_below_0.05",sum(tab$padj<.05,na.rm=TRUE)),paste("R",R.version.string),paste("DESeq2",as.character(utils::packageVersion("DESeq2"))),paste("filter","counts >= 10 in >= 2 samples")),file.path(output_dir,"run-log.txt"))
message("Saved results in ",output_dir)
