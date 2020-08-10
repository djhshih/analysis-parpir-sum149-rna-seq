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
library(ggsci)

source("../R/preamble.R")

# cell line: SUM149

#in.fname <- as.filename("parpi-resist_ltpm-genes-max.rds");
in.fname <- as.filename("parpi-resist_ltpm-genes-mean.rds");
out.fname <- in.fname;
out.fname$ext <- NULL;

mc.cores <- 4;


x <- qread(in.fname);
pheno <- setup_pheno(qread("../sample-info_parpi-resist_stage2.tsv"));

# arrange samples in order to clones
idx <- order(pheno$clone);
x <- x[, idx];
pheno <- pheno[idx, ];

stopifnot(colnames(x) == as.character(pheno$sample_id))

gsets.h <- read_msigdb("h.all");

es.h <- gsva(x, gsets.h$data, parallel.sz=mc.cores);
#es.h <- gsva(x, gsets.h$data, parallel.sz=mc.cores, method="ssgsea");
#es.h <- gsva(x, gsets.h$data, parallel.sz=mc.cores, method="zscore");
# KRAS signaling status gets reversed using plage...??!
#es.h <- gsva(x, gsets.h$data, parallel.sz=mc.cores, method="plage");

ha <- HeatmapAnnotation(
	df = select(pheno,
		clone
	),
	col = list(
		clone = clone.cols
	),
	annotation_legend_param = list(
		clone = list(nrow = 1)
	),
	gp = gpar(col = "white")
);

pdf(tag(out.fname, c("gsva", "h"), ext="pdf"), width=8, height=10);
hm <- Heatmap(
	rename_hallmarks(es.h),
	rect_gp = gpar(col = "white"),
	top_annotation = ha,
	column_gap = unit(3, "mm"),
	clustering_method_rows = "average",
	clustering_distance_rows = "pearson",
	column_split = pheno$treatment,
	cluster_columns = FALSE,
	show_column_names = FALSE,
	row_dend_width = unit(20, "mm"),
	col = rev(brewer.pal(9, "RdBu")),
	border = TRUE,
	heatmap_legend_param = list(direction = "horizontal"),
	name = "enrichment score"
);
draw(hm, merge_legend=TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off();


gsets.c6 <- read_msigdb("c6.all");

es.c6 <- gsva(x, gsets.c6$data, parallel.sz=mc.cores);

Heatmap(es.c6, top_annotation = ha, cluster_columns = FALSE, row_names_gp = gpar(fontsize=4))


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
es.c6.up <- es.c6[order(resistance.t)[1:n.sub], ]

pdf(tag(out.fname, c("gsva", "c6", "parental-up"), ext="pdf"), width=10, height=15);
Heatmap(es.c6.up, top_annotation = ha, cluster_columns = FALSE, row_names_gp = gpar(fontsize=8));
dev.off();

# downregulated in parental, upregulated in resistant clones
es.c6.down <- es.c6[order(-resistance.t)[1:n.sub], ]

pdf(tag(out.fname, c("gsva", "c6", "parental-down"), ext="pdf"), width=10, height=15);
Heatmap(es.c6.down, top_annotation = ha, cluster_columns = FALSE, row_names_gp = gpar(fontsize=8))
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

