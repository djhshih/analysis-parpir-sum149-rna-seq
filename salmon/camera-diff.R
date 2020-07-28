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
#in.fname <- as.filename("parpi-resist_deseq-stat_treatment-clone-interaction_clone_treated-vs-untreated.mtx");
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

na <- apply(x, 1, function(z) any(is.na(z)));
y <- x[!na, ];

gsets.h <- read_msigdb("h.all");

#d <- data.frame(comparison = c(colnames(x), "Parental"));
#d$reference <- 0;
#d$reference[d$comparison != "Parental"] <- 1;

#design <- model.matrix(~ reference, data=d);

#gset <- gsets.h$data$HALLMARK_KRAS_SIGNALING_UP;

camera_transform <- function(y, gsets) {
	index <- lapply(gsets,
		function(gset) {
			rownames(y) %in% gset
		}
	);

	res <- apply(y, 2, function(s) {
		#d <- cameraPR(s, index, use.ranks=TRUE);
		d <- cameraPR(s, index, sort=FALSE);
		z <- log10(d$FDR);

		# flip the sign of enriched results
		idx <- d$Direction == "Up";
		z[idx] <- -z[idx]

		z
	});
	
	rownames(res) <- names(index);

	res
}

es.h <- camera_transform(y, gsets.h$data);

summary(es.h)

colf = circlize::colorRamp2(c(-10, 0, 10), c("blue", "white", "red"));


pdf(tag(out.fname, c("camera", "h"), ext="pdf"), width=10, height=10);
Heatmap(es.h, col=colf, cluster_col = FALSE)
dev.off();
#Heatmap(es.h, top_annotation = ha, cluster_col = TRUE)

es.h[grep("TNFA", rownames(es.h)), ]

library(ggplot2)

plot_gene_set <- function(x, genes) {
	y.sel <- y[rownames(y) %in% genes, ];
	y.sel.m <- melt(y.sel, varnames=c("gene", "comparison"));
	ggplot(y.sel.m, aes(x=value, fill=comparison)) +
		geom_density(alpha=0.5) +
		facet_grid(comparison ~ .)
}

genes <- gsets.h$data[[grep("TNFA", names(gsets.h$data))]];
plot_gene_set(x, genes);

genes <- gsets.h$data[[grep("OXIDATIVE", names(gsets.h$data))]];
plot_gene_set(x, genes);

genes <- gsets.h$data[[grep("KRAS_SIGNALING_UP", names(gsets.h$data))]];
plot_gene_set(x, genes);

genes <- gsets.h$data[[grep("MESENCHYMAL", names(gsets.h$data))]];
plot_gene_set(x, genes);



gsets.c6 <- read_msigdb("c6.all");

es.c6 <- camera_transform(y, gsets.c6$data);

Heatmap(es.c6, col=colf, cluster_col = FALSE, row_names_gp = gpar(fontsize=4))

dim(es.c6)
hist(es.c6)
summary(es.c6)

means <- rowMeans(es.c6);
summary(means)
n.sub <- 80;

# upregulated in resistant clones
#es.c6.sub <- es.c6[means > 0, , drop=FALSE]
es.c6.sub <- es.c6[means > 1, , drop=FALSE]

ha <- NULL;

pdf(tag(out.fname, c("camera", "c6", "parental-up"), ext="pdf"), width=10, height=15);
Heatmap(es.c6.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));
dev.off();

# downregulated in resistant clones
#es.c6.sub <- es.c6[order(means)[1:n.sub], , drop=FALSE]
#es.c6.sub <- es.c6[means < -2, , drop=FALSE]
es.c6.sub <- es.c6[means < -1, , drop=FALSE]

pdf(tag(out.fname, c("camera", "c6", "parental-down"), ext="pdf"), width=10, height=15);
Heatmap(es.c6.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8))
dev.off();

y["BMI1", ]

# downregulate EGFR signaling
# downregulating NFkappaB signaling
# BMI downregulation
# suppress response to BMI1 downregulation
# downregulating KRAS signaling
# suppress gene expression associated with KRAS dependency

# suppress signaling in response to STK33 knockdown
# suppress IL2, IL15 signaling



gsets.c7 <- read_msigdb("c7.all");

es.c7 <- camera_transform(y, gsets.c7$data);

summary(es.c7)

means <- rowMeans(es.c7)

es.c7.sub <- es.c7[order(-means)[1:n.sub], , drop=FALSE]
Heatmap(es.c7.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));

es.c7.sub <- es.c7[order(means)[1:n.sub], , drop=FALSE]
Heatmap(es.c7.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));

#

gsets.c3 <- read_msigdb("c3.tft");

es.c3 <- camera_transform(y, gsets.c3$data);

summary(es.c3)

means <- rowMeans(es.c3);

es.c3.sub <- es.c3[order(-means)[1:n.sub], , drop=FALSE]
Heatmap(es.c3.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));

es.c3.sub <- es.c3[order(means)[1:n.sub], , drop=FALSE]
#es.c3.sub <- es.c3[means < -1, , drop=FALSE]
Heatmap(es.c3.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));

#

gsets.c2 <- read_msigdb("c2.cgp");

es.c2 <- camera_transform(y, gsets.c2$data);

summary(es.c2)

means <- rowMeans(es.c2);

#es.c2.sub <- es.c2[order(-means)[1:n.sub], , drop=FALSE]
es.c2.sub <- es.c2[means > 2, , drop=FALSE]
Heatmap(es.c2.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));

#es.c2.sub <- es.c2[order(means)[1:n.sub], , drop=FALSE]
es.c2.sub <- es.c2[means < -2, , drop=FALSE]
Heatmap(es.c2.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));

#

gsets.c2 <- read_msigdb("c2.cp");

es.c2 <- camera_transform(y, gsets.c2$data);

summary(es.c2)

means <- rowMeans(es.c2);
#maxes <- apply(es.c2, 1, function(z) max(abs(z)));

#es.c2.sub <- es.c2[order(-means)[1:n.sub], , drop=FALSE]
es.c2.sub <- es.c2[means > 2, , drop=FALSE]
Heatmap(es.c2.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));

#es.c2.sub <- es.c2[order(means)[1:n.sub], , drop=FALSE]
es.c2.sub <- es.c2[means < -2, , drop=FALSE]
Heatmap(es.c2.sub, col=colf, top_annotation = ha, cluster_col = FALSE, cluster_row = FALSE, row_names_gp = gpar(fontsize=8));

