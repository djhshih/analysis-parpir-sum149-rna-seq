library(io)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(ggsci)
library(ggpubr)

source("../R/preamble.R")
source("../R/deseq.R")

# cell line: SUM149

in.fname <- as.filename("parpi-resist_deseq-res_resistance.rds");
out.fname <- in.fname;
out.fname$ext <- NULL;

res <- qread(in.fname);
pheno <- setup_pheno(qread("../sample-info_parpi-resist_stage2.tsv"), rename.clones=TRUE);
fannot <- qread("parpi-resist_fannot.rds");

####

gsets.h <- read_msigdb("h.all");

gset.emt <- gsets.h$data$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION;
gset.kras <- gsets.h$data$HALLMARK_KRAS_SIGNALING_UP;

intersect(gset.emt, gset.kras)


des <- deseq_summarize_extreme(res, fannot);

des.a <- deseq_annotate_significant(des, fdr.cut=0.01, delta.cut=0);
des.k <- filter(des.a, keep);

qdraw(
	deseq_plot_stat(filter(des.k, gene_name %in% gset.emt)) +
		ggtitle("RCs vs. P, EMT pathway")
	,
	height = 12, width = 4,
	file = insert(out.fname, c("clones-vs-parental", "genes", "emt"), ext="pdf")
)

qdraw(
	deseq_plot_stat(filter(des.k, gene_name %in% gset.kras)) +
		ggtitle("RCs vs. P, KRAS signaling")
	,
	height = 10, width = 4,
	file = insert(out.fname, c("clones-vs-parental", "genes", "kras"), ext="pdf")
)

