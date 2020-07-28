library(io)
library(GSVA)
library(mmalign)
library(dplyr)
library(reshape2)
library(GSVA)
library(parallel)
library(ComplexHeatmap)
library(RColorBrewer)
library(sva)

# cell line: SUM149

in.fname <- as.filename("parpi-resist_ltpm-genes-max.rds");
out.fname <- in.fname;
out.fname$ext <- NULL;

x <- qread(in.fname);
pheno <- qread("../sample-info_parpi-resist_stage2.tsv");
source("relevel-pheno.R");

stopifnot(colnames(x) == as.character(pheno$sample_id))

mc.cores <- 4;

read_msigdb <- function(collection, version="6.2") {
	release <- gsub(".", "", version, fixed=TRUE);
	qread(sprintf("~/data/msigdb/release-%s/%s.v%s.symbols.gmt", release, collection, version))
}

clone.cols <- brewer.pal(8, "Accent");
names(clone.cols) <- levels(pheno$clone);

ha <- HeatmapAnnotation(
	df = select(pheno, clone, treatment, batch, lane, fold_resistance),
	col = list(
		treatment = c("DMSO" = "grey30", "None" = "grey60", "Talazoparib" = "royalblue"),
		clone = clone.cols
	)
);


gsets.h <- read_msigdb("h.all");

es.h <- gsva(x, gsets.h$data, parallel.sz=mc.cores);
#es.h <- gsva(x, gsets.h$data, parallel.sz=mc.cores, method="ssgsea");
#es.h <- gsva(x, gsets.h$data, parallel.sz=mc.cores, method="zscore");
# KRAS signaling status gets reversed using plage...??!
#es.h <- gsva(x, gsets.h$data, parallel.sz=mc.cores, method="plage");

pdf(tag(out.fname, c("gsva", "h"), ext="pdf"), width=10, height=10);
Heatmap(es.h, top_annotation = ha, cluster_col = FALSE)
dev.off();
#Heatmap(es.h, top_annotation = ha, cluster_col = TRUE)


gsets.c6 <- read_msigdb("c6.all");

es.c6 <- gsva(x, gsets.c6$data, parallel.sz=mc.cores);

Heatmap(es.c6, top_annotation = ha, cluster_col = FALSE, row_names_gp = gpar(fontsize=4))


compare <- function(formula, x, pheno) {
	t(apply(x, 1, function(z) {
		# get the t statistic
		coef(summary(lm(formula, data = cbind(pheno, x=z))))[,3]
	}))
}

resistance.t <- compare(x ~ resistance, es.c6, pheno)[, 2];

n.sub <- 80;

dim(es.c6)

# upregulated in parental, downregulated in resistant clones
es.c6.sub <- es.c6[order(resistance.t)[1:n.sub], ]

pdf(tag(out.fname, c("gsva", "c6", "parental-up"), ext="pdf"), width=10, height=15);
Heatmap(es.c6.sub, top_annotation = ha, cluster_col = FALSE, row_names_gp = gpar(fontsize=8));
dev.off();

# downregulated in parental, upregulated in resistant clones
es.c6.sub <- es.c6[order(-resistance.t)[1:n.sub], ]

pdf(tag(out.fname, c("gsva", "c6", "parental-down"), ext="pdf"), width=10, height=15);
Heatmap(es.c6.sub, top_annotation = ha, cluster_col = FALSE, row_names_gp = gpar(fontsize=8))
dev.off();

# resistant clones have expression patterns indictative of
# - KRAS deactivation
# - beta-catenin deactivation
# - BMI1 activation
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

