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
if (anyDuplicated(names(x)) || anyDuplicated(names(m)) || anyDuplicated(x$gene_id) || anyDuplicated(m$sample_id) || anyNA(x) || anyNA(m)) stop("Duplicate identifiers or missing values")
if (ncol(x)<3 || nrow(x)<2 || nrow(m)<2) stop("Too few genes or samples")
if (!setequal(names(x)[-1],m$sample_id)) stop("Count and metadata sample IDs differ")
m <- m[match(names(x)[-1], m$sample_id), , drop=FALSE]
if (!identical(m$sample_id,names(x)[-1])) stop("Sample alignment failed")
cnt <- as.matrix(x[-1]); suppressWarnings(storage.mode(cnt)<-"numeric")
if (anyNA(cnt) || any(!is.finite(cnt)) || any(cnt<0) || any(cnt!=floor(cnt))) stop("Counts must be nonnegative finite integers")
if (any(!nzchar(trimws(x$gene_id))) || any(!nzchar(trimws(m$sample_id)))) stop("Blank identifiers")
if (any(colSums(cnt)==0)) stop("All-zero sample")
rownames(cnt) <- x$gene_id
mode <- if ("dex" %in% names(m) && "cell" %in% names(m)) "airway" else if ("group" %in% names(m)) "synthetic" else stop("Need group or both cell and dex columns")
if (mode=="synthetic") {
 if (!all(c("A","B") %in% m$group) || any(!m$group %in% c("A","B")) || min(table(m$group))<2) stop("Synthetic example needs replicated A and B groups")
 m$group <- factor(m$group,levels=c("A","B"))
 if ("donor" %in% names(m)) {
   m$donor <- factor(m$donor)
   if (nlevels(m$donor)<2 || any(table(m$donor,m$group)!=1)) stop("Each donor needs one A and one B sample")
   design <- ~ donor + group
 } else design <- ~ group
} else {
 if (!all(c("untrt","trt") %in% m$dex) || any(!m$dex %in% c("untrt","trt"))) stop("Expected airway dex values untrt and trt")
 if (any(table(m$cell,m$dex)!=1)) stop("Each airway cell line must have one sample per condition")
 m$cell <- factor(m$cell); m$dex <- factor(m$dex,levels=c("untrt","trt")); design <- ~ cell + dex
}
rownames(m) <- m$sample_id
model <- model.matrix(design, m)
if (qr(model)$rank<ncol(model) || nrow(model)<=ncol(model)) stop("Design is confounded or has no residual degrees of freedom")
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
writeLines(c(paste("input",normalizePath(data_dir)),paste("design",deparse(design)),paste("contrast",paste(contrast,collapse=" ")),paste("genes_tested",nrow(tab)),paste("padj_below_0.05",sum(tab$padj<.05,na.rm=TRUE)),paste("R",R.version.string),paste("DESeq2",as.character(utils::packageVersion("DESeq2"))),paste("filter","counts >= 10 in >= 2 samples")),file.path(output_dir,"run-log.txt"))

# 6. Save quantities behind all six teaching figures -----------------------
condition <- if (mode=="airway") as.character(m$dex) else as.character(m$group)
subject <- if (mode=="airway") as.character(m$cell) else if ("donor" %in% names(m)) as.character(m$donor) else NULL
levels_condition <- if (mode=="airway") c("untrt","trt") else c("A","B")
normal <- DESeq2::counts(dds, normalized=TRUE)
write_matrix <- function(mat, filename, id="gene_id") {
 out <- data.frame(rownames(mat), mat, check.names=FALSE); names(out)[1] <- id
 write.table(out,file.path(output_dir,"tables",filename),sep="\t",quote=FALSE,row.names=FALSE)
}
write_matrix(normal,"normalized-counts.tsv")
write_matrix(SummarizedExperiment::assay(vsd),"transformed-counts.tsv")
sf <- data.frame(sample_id=colnames(dds),size_factor=DESeq2::sizeFactors(dds))
write.table(sf,file.path(output_dir,"tables/size-factors.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(pca,file.path(output_dir,"tables/pca-scores.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(data.frame(component=c("PC1","PC2"),percent_variance=100*attr(pca,"percentVar")),file.path(output_dir,"tables/pca-variance.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
qc$condition <- factor(condition,levels=levels_condition)
qc$subject <- if(is.null(subject)) qc$sample_id else subject
save_plot <- function(plot,name,width=8,height=5) ggplot2::ggsave(file.path(output_dir,"figures",name),plot,width=width,height=height,dpi=180)
p_lib <- ggplot2::ggplot(qc,ggplot2::aes(x=sample_id,y=total_counts,fill=condition)) +
 ggplot2::geom_col() + ggplot2::labs(x="Sample",y="Sum of supplied raw gene counts",fill="Condition") +
 ggplot2::theme_minimal(base_size=12) + ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))
save_plot(p_lib,"library-sizes.png")
# PCA uses up to 500 most variable transformed genes. Donor effects remain.
pca$condition <- factor(condition,levels=levels_condition)
pca$subject <- if(is.null(subject)) "Unpaired" else subject
p_pca <- ggplot2::ggplot(pca,ggplot2::aes(PC1,PC2,color=condition))
if(!is.null(subject)) {
 p_pca <- p_pca + ggplot2::geom_line(ggplot2::aes(group=subject),color="grey70") + ggplot2::geom_point(ggplot2::aes(shape=subject),size=3)
} else {
 p_pca <- p_pca + ggplot2::geom_point(size=3)
}
p_pca <- p_pca + ggplot2::geom_text(ggplot2::aes(label=name),vjust=-.7,size=3,show.legend=FALSE) +
 ggplot2::labs(x=paste0("PC1: ",variance[1],"%"),y=paste0("PC2: ",variance[2],"%"),color="Condition",shape="Donor / cell line",subtitle="Variance-stabilized counts; grey lines connect paired samples") + ggplot2::theme_minimal(base_size=12)
p_pca <- p_pca + ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=.18))
save_plot(p_pca,"sample-pca.png",9,6)
distance <- as.matrix(dist(t(SummarizedExperiment::assay(vsd))))
write_matrix(distance,"sample-distances.tsv","sample_id")
ordering <- hclust(as.dist(distance))$order
distance_labels <- paste(m$sample_id,condition,if(is.null(subject)) "" else subject,sep=" | ")
dimnames(distance) <- list(distance_labels,distance_labels)
heat <- expand.grid(sample_x=distance_labels,sample_y=distance_labels,stringsAsFactors=FALSE)
heat$distance <- as.vector(distance)
heat$sample_x <- factor(heat$sample_x,levels=distance_labels[ordering])
heat$sample_y <- factor(heat$sample_y,levels=rev(distance_labels[ordering]))
p_dist <- ggplot2::ggplot(heat,ggplot2::aes(sample_x,sample_y,fill=distance))+ggplot2::geom_tile()+
 ggplot2::scale_fill_gradient(low="white",high="#174a71")+ggplot2::coord_equal()+
 ggplot2::labs(x=NULL,y=NULL,fill="Euclidean\ndistance",subtitle="All retained genes after variance stabilization; clustered sample order")+
 ggplot2::theme_minimal(base_size=10)+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=60,hjust=1),panel.grid=ggplot2::element_blank())
if (mode=="synthetic" && !is.null(subject)) {
 y_levels <- c(levels(heat$sample_y),"Condition","Donor")
 heat$sample_y <- factor(as.character(heat$sample_y),levels=y_levels)
 ann <- data.frame(sample_x=factor(rep(distance_labels,2),levels=levels(heat$sample_x)),
                   sample_y=factor(rep(c("Condition","Donor"),each=nrow(m)),levels=y_levels),
                   value=c(condition,subject))
 donor_colours <- setNames(grDevices::hcl.colors(length(unique(subject)),"Dark 3"),unique(subject))
 annotation_colours <- c(c(A="#f8766d",B="#00bfc4")[condition],donor_colours[subject])
 p_dist <- ggplot2::ggplot(heat,ggplot2::aes(sample_x,sample_y,fill=distance))+
   ggplot2::geom_tile()+ggplot2::scale_fill_gradient(low="white",high="#174a71")+
   ggplot2::geom_tile(data=ann,ggplot2::aes(sample_x,sample_y),fill=unname(annotation_colours),inherit.aes=FALSE)+
   ggplot2::geom_text(data=ann,ggplot2::aes(sample_x,sample_y,label=value),inherit.aes=FALSE,size=2.8)+
   ggplot2::coord_equal()+ggplot2::labs(x=NULL,y=NULL,fill="Euclidean\ndistance",
     subtitle="All retained variance-stabilized genes; condition and donor strips")+
   ggplot2::theme_minimal(base_size=12)+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=60,hjust=1),panel.grid=ggplot2::element_blank())
}
save_plot(p_dist,"sample-distances.png",10,9)
tab$status <- ifelse(is.na(tab$padj),"padj unavailable",ifelse(tab$padj<.05,"padj < 0.05","padj >= 0.05"))
colours <- c("padj unavailable"="#b8b8b8","padj < 0.05"="#b63848","padj >= 0.05"="#326b88")
finite_fc <- is.finite(tab$log2FoldChange)
p_ma <- ggplot2::ggplot(tab[finite_fc & tab$baseMean>0,],ggplot2::aes(baseMean,log2FoldChange,color=status))+
 ggplot2::geom_point(size=1,alpha=.6)+ggplot2::geom_hline(yintercept=0,linetype=2)+ggplot2::scale_x_log10()+
 ggplot2::scale_color_manual(values=colours)+ggplot2::labs(x="Mean normalized count (log scale)",y="Unshrunken log2 fold change",color=NULL,subtitle=paste(contrast[2],"versus",contrast[3],"; positive = higher in",contrast[2]))+ggplot2::theme_minimal(base_size=12)
save_plot(p_ma,"ma-plot.png")
v <- tab[finite_fc & !is.na(tab$padj),]
v$minus_log10_padj <- -log10(pmax(v$padj,.Machine$double.xmin))
p_volcano <- ggplot2::ggplot(v,ggplot2::aes(log2FoldChange,minus_log10_padj,color=status))+
 ggplot2::geom_point(size=1,alpha=.6)+ggplot2::geom_hline(yintercept=-log10(.05),linetype=2)+ggplot2::geom_vline(xintercept=0,linetype=2)+
 ggplot2::scale_color_manual(values=colours)+ggplot2::labs(x="Unshrunken log2 fold change",y="-log10(adjusted p-value)",color=NULL,subtitle="Missing adjusted p-values omitted; zero values floored for display")+ggplot2::theme_minimal(base_size=12)
save_plot(p_volcano,"volcano.png")
# Choose deterministically for teaching, not as independent validation.
valid <- tab$gene_id[!is.na(tab$padj) & finite_fc]
selected <- if(length(valid)) valid[1] else rownames(normal)[which.max(rowMeans(normal))]
gene_data <- data.frame(gene_id=selected,sample_id=colnames(normal),normalized_count=as.numeric(normal[selected,]),condition=factor(condition,levels=levels_condition),subject=if(is.null(subject)) m$sample_id else subject)
write.table(gene_data,file.path(output_dir,"tables/selected-gene-counts.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
p_gene <- ggplot2::ggplot(gene_data,ggplot2::aes(condition,normalized_count))
if(!is.null(subject)) p_gene <- p_gene+ggplot2::geom_line(ggplot2::aes(group=subject,color=subject))
p_gene <- p_gene+ggplot2::geom_point(ggplot2::aes(color=subject),size=3)+
 ggplot2::labs(x="Condition",y="DESeq2 normalized count",color="Donor / sample",title=selected,subtitle="Lowest available adjusted p-value; selected after analysis for illustration")+ggplot2::theme_minimal(base_size=12)
# Optional prespecified examples supplement the post-analysis selected gene.
# No example is chosen or discarded because of its fitted p-value.
example_path <- file.path(data_dir,"plot-examples.tsv")
if (mode=="synthetic" && !is.null(subject) && file.exists(example_path)) {
 examples <- read_tsv(example_path)
 if (!all(c("gene_id","example") %in% names(examples)) || anyDuplicated(examples$gene_id)) stop("Invalid plot-examples.tsv")
 if (!all(examples$gene_id %in% rownames(normal))) stop("A prespecified example did not survive filtering; inspect inputs")
 examples <- rbind(data.frame(gene_id=selected,example="Selected after analysis"),examples[examples$gene_id!=selected,c("gene_id","example")])
 panels <- do.call(rbind,lapply(seq_len(nrow(examples)),function(i) {
   out <- data.frame(gene_id=examples$gene_id[i],sample_id=colnames(normal),normalized_count=as.numeric(normal[examples$gene_id[i],]),
     condition=factor(condition,levels=levels_condition),subject=subject)
   out$example <- examples$example[i]
   out
 }))
 write.table(panels,file.path(output_dir,"tables/example-gene-counts.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
 panels$panel <- paste(panels$gene_id,panels$example,sep="\n")
 panels$panel <- factor(panels$panel,levels=unique(panels$panel))
 p_gene <- ggplot2::ggplot(panels,ggplot2::aes(condition,normalized_count,group=subject,color=subject))+
   ggplot2::geom_line()+ggplot2::geom_point(size=2.5)+ggplot2::facet_wrap(~panel,nrow=1,scales="free_y")+
   ggplot2::labs(x="Condition",y="DESeq2 normalized count",color="Donor",
     subtitle="One post-analysis selection; two examples specified in the simulation")+
   ggplot2::theme_minimal(base_size=12)+ggplot2::theme(legend.position="bottom")
 save_plot(p_gene,"selected-gene.png",12,5.5)
} else {
 save_plot(p_gene,"selected-gene.png")
}
# Inputs and script hashes tie preview files to the exact run.
capture.output(sessionInfo(),file=file.path(output_dir,"session-info.txt"))
hash_inputs <- c(file.path(data_dir,c("gene-counts.tsv","sample-metadata.tsv")),"scripts/R/run_workflow.R")
if (mode=="synthetic" && file.exists(example_path)) hash_inputs <- c(hash_inputs,example_path)
if (mode=="synthetic" && file.exists(file.path(data_dir,"simulation-record.txt"))) hash_inputs <- c(hash_inputs,"scripts/R/prepare_synthetic.R",file.path(data_dir,"simulation-record.txt"))
hashes <- tools::md5sum(hash_inputs)
write.table(data.frame(path=names(hashes),md5=unname(hashes)),file.path(output_dir,"tables/input-script-hashes.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
cat(paste0("selected_gene ",selected,"\nPCA_top_genes ",min(500,nrow(dds)),"\nVST blind=FALSE; donor effects are not removed\nlog2FC unshrunken; volcano uses adjusted p-values\n"),file=file.path(output_dir,"run-log.txt"),append=TRUE)
message("All six preview figures and supporting tables saved in ",output_dir)
