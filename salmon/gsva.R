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
pheno <- setup_pheno(qread("../sample-info_parpi-resist_stage2.tsv"), rename.clones=TRUE);

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

