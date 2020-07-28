library(io)
library(GSVA)
library(dplyr)
library(reshape2)
library(parallel)
library(ComplexHeatmap)
library(RColorBrewer)
library(sva)

# cell line: SUM149
# triple negative, inflammatory breast cancer
# disease progressed through chemotherapy
# mutation: BRCA1 p.P724fs*12
# microamplification in PTEN exon 2 leads to truncated PTEN transcript
# beginning in exon 2; no RNA expression of PTEN for all subsequent exons

# TODO verify PTEN loss of function in RNAseq data
# TODO check BRCA1 mutation status!

in.fname <- as.filename("parpi-resist_deseq-stat_treatment-clone-interaction_clones-vs-parental.mtx");
out.fname <- in.fname;
out.fname$ext <- NULL;

pheno <- qread("../sample-info_parpi-resist_stage2.tsv");

x <- qread(in.fname);

mc.cores <- 4;

read_msigdb <- function(collection, version="6.2") {
	release <- gsub(".", "", version, fixed=TRUE);
	qread(sprintf("~/data/msigdb/release-%s/%s.v%s.symbols.gmt", release, collection, version))
}

clone.cols <- brewer.pal(8, "Accent");
names(clone.cols) <- levels(pheno$clone);

#ha <- HeatmapAnnotation(
#	df = select(pheno, clone, treatment, batch, lane, fold_resistance),
#	col = list(
#		treatment = c("DMSO" = "grey30", "None" = "grey60", "Talazoparib" = "royalblue"),
#		clone = clone.cols
#	)
#);


gsets.h <- read_msigdb("h.all");

#es.h <- gsva(x, gsets.h$data, parallel.sz=mc.cores, kcdf="none");
es.h <- gsva(x, gsets.h$data, parallel.sz=mc.cores, method="ssgsea");

summary(es.h)

pdf(tag(out.fname, c("gsva", "h"), ext="pdf"), width=10, height=10);
Heatmap(es.h, cluster_col = FALSE)
dev.off();
#Heatmap(es.h, top_annotation = ha, cluster_col = TRUE)


gsets.c6 <- read_msigdb("c6.all");

es.c6 <- gsva(x, gsets.c6$data, parallel.sz=mc.cores, method="ssgsea");

colf = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"));

Heatmap(es.c6, col=colf, cluster_col = FALSE, row_names_gp = gpar(fontsize=4))

dim(es.c6)
hist(es.c6)

means <- rowMeans(es.c6);
n.sub <- 80;

# upregulated in resistant clones
es.c6.sub <- es.c6[means > 0, , drop=FALSE]

ha <- NULL;

pdf(tag(out.fname, c("gsva", "c6", "parental-up"), ext="pdf"), width=10, height=15);
Heatmap(es.c6.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));
dev.off();

# downregulated in resistant clones
es.c6.sub <- es.c6[order(means)[1:n.sub], , drop=FALSE]

pdf(tag(out.fname, c("gsva", "c6", "parental-down"), ext="pdf"), width=10, height=15);
Heatmap(es.c6.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8))
dev.off();

# resistant clones have expression patterns indictative of
# - NFkappaB deactivation
# - downregulation of genes associated with KRAS dependence
# - BMI1 activation
# - beta-catenin deactivation
# - GLI signaling deactivation
# - EGFR signaling deactivation
# - KRAS deactivation
# but not every clone response the same!

# heterogeneous expression patterns associated with
# - PTEN signaling
# - MTRO signaling

# "Loss of Mel-18 enhances breast cancer stem cell activity and tumorigenicity through activating Notch signaling mediated by the Wnt/TCF pathway"
# MEL18 (BMI1) -| Wnt signaling -> stemness

# "β-Catenin Is Required for the Tumorigenic Behavior of Triple-Negative Breast Cancer Cells"
# Beta-catenin knockdown -| stem-like population
# Beta-catenin knockdown -> Doxorubicin sensitivity

# Olaparib + liposomal doxorubicin clinical trial in ovarian cancer:
# NCT03161132

gsets.c7 <- read_msigdb("c7.all");

es.c7 <- gsva(x, gsets.c7$data, parallel.sz=mc.cores, method="ssgsea");

means <- rowMeans(es.c7)

es.c7.sub <- es.c7[means > 0, , drop=FALSE]
Heatmap(es.c7.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));

es.c7.sub <- es.c7[order(means)[1:n.sub], , drop=FALSE]
Heatmap(es.c7.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));

#

gsets.c3 <- read_msigdb("c3.tft");

es.c3 <- gsva(x, gsets.c3$data, parallel.sz=mc.cores, method="ssgsea");

summary(es.c3)

means <- rowMeans(es.c3);

es.c3.sub <- es.c3[order(means)[1:n.sub], , drop=FALSE]
Heatmap(es.c3.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));

es.c3.sub <- es.c3[order(-means)[1:n.sub], , drop=FALSE]
Heatmap(es.c3.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));

#

gsets.c2 <- read_msigdb("c2.cgp");

es.c2 <- gsva(x, gsets.c2$data, parallel.sz=mc.cores, method="ssgsea");

summary(es.c2)

means <- rowMeans(es.c2);

es.c2.sub <- es.c2[order(means)[1:n.sub], , drop=FALSE]
Heatmap(es.c2.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));

es.c2.sub <- es.c2[order(-means)[1:n.sub], , drop=FALSE]
Heatmap(es.c2.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));

