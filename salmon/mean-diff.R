library(io)
library(dplyr)
library(reshape2)
library(limma)
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

mc.cores <- 4;

pheno <- setup_pheno(qread("../annot/sample-info_parpi-resist_stage2.tsv"));
x <- qread(in.fname);

#ha <- HeatmapAnnotation(
#	df = select(pheno, clone, treatment, batch, lane, fold_resistance),
#	col = list(
#		treatment = c("DMSO" = "grey30", "None" = "grey60", "Talazoparib" = "royalblue"),
#		clone = clone.cols
#	)
#);

na <- apply(x, 1, function(z) any(is.na(z)));
y <- x[!na, ];

gsets.h <- read_msigdb("h.all");

# gene set z transform
gsetz_transform <- function(x, gsets) {
	index <- lapply(gsets,
		function(gset) {
			rownames(x) %in% gset
		}
	);

	res <- apply(x, 2, function(s) {
		unlist(lapply(index, function(idx) {
			si <- s[idx];
			z <- mean(si) / sd(si);
			z
		}))
	});
	rownames(res) <- names(index);

	res
}

es.h <- gsetz_transform(y, gsets.h$data);

summary(es.h)

colf = circlize::colorRamp2(c(-3, 0, 3), c("blue", "white", "red"));


pdf(tag(out.fname, c("gsetz", "h"), ext="pdf"), width=10, height=10);
Heatmap(es.h, col=colf, cluster_col = FALSE)
dev.off();
#Heatmap(es.h, top_annotation = ha, cluster_col = TRUE)


gsets.c6 <- read_msigdb("c6.all");

es.c6 <- gsetz_transform(y, gsets.c6$data);

Heatmap(es.c6, col=colf, cluster_col = FALSE, row_names_gp = gpar(fontsize=4))

dim(es.c6)
hist(es.c6)
summary(es.c6)

means <- rowMeans(es.c6);
n.sub <- 80;

# upregulated in resistant clones
es.c6.sub <- es.c6[order(-means)[1:n.sub], , drop=FALSE]

ha <- NULL;

pdf(tag(out.fname, c("gsetz", "c6", "parental-up"), ext="pdf"), width=10, height=15);
Heatmap(es.c6.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));
dev.off();

# downregulated in resistant clones
es.c6.sub <- es.c6[order(means)[1:n.sub], , drop=FALSE]

pdf(tag(out.fname, c("gsetz", "c6", "parental-down"), ext="pdf"), width=10, height=15);
Heatmap(es.c6.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8))
dev.off();

genes <- gsets.h$data[[grep("GLI1_UP.V1_DN", names(gsets.c6$data))]];
plot_gene_set(x, genes);

genes <- gsets.h$data[[grep("GLI1_UP.V1_UP", names(gsets.c6$data))]];
plot_gene_set(x, genes);

apply(x[rownames(x) %in% genes, ], 2, mean, na.rm=TRUE)


gsets.c7 <- read_msigdb("c7.all");

es.c7 <- gsetz_transform(y, gsets.c7$data);

summary(es.c7)

means <- rowMeans(es.c7)

es.c7.sub <- es.c7[order(-means)[1:n.sub], , drop=FALSE]
Heatmap(es.c7.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));

es.c7.sub <- es.c7[order(means)[1:n.sub], , drop=FALSE]
Heatmap(es.c7.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));

#

gsets.c3 <- read_msigdb("c3.tft");

es.c3 <- gsetz_transform(y, gsets.c3$data);

summary(es.c3)

means <- rowMeans(es.c3);

es.c3.sub <- es.c3[order(-means)[1:n.sub], , drop=FALSE]
Heatmap(es.c3.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));

es.c3.sub <- es.c3[order(means)[1:n.sub], , drop=FALSE]
Heatmap(es.c3.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));

#

gsets.c2 <- read_msigdb("c2.cgp");

es.c2 <- gsetz_transform(y, gsets.c2$data);

summary(es.c2)

means <- rowMeans(es.c2);

es.c2.sub <- es.c2[order(-means)[1:n.sub], , drop=FALSE]
Heatmap(es.c2.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));

es.c2.sub <- es.c2[order(means)[1:n.sub], , drop=FALSE]
Heatmap(es.c2.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));

