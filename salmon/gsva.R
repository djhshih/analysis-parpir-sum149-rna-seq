library(io)
library(GSVA)
library(mmalign)
library(dplyr)
library(tidyr)
library(reshape2)
library(GSVA)
library(ComplexHeatmap)
library(RColorBrewer)
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

# arrange samples in order of clones
idx <- order(pheno$clone);
x <- x[, idx];
pheno <- pheno[idx, ];

stopifnot(colnames(x) == as.character(pheno$sample_id))

####

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

d.pathways <- data.frame(
	sample_id = colnames(es.h),
	kras = es.h["HALLMARK_KRAS_SIGNALING_UP", ],
	emt  = es.h["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", ]
);
rownames(d.pathways) <- NULL;
qwrite(d.pathways, insert(out.fname, c("gsva", "kras", "emt"), ext="tsv"))

d <- d.pathways %>% left_join(pheno) %>% filter(treatment != "None");

ggplot(d, aes(x = kras, y = emt, colour = clone, shape=treatment)) + theme_clean() +
	geom_point() +
	scale_colour_manual(values=clone.cols)

d.sum <- select(d, clone, treatment, kras, emt) %>% group_by(clone, treatment) %>%
	summarize(kras = mean(kras), emt = mean(emt)) %>% ungroup();

qwrite(d.sum, insert(out.fname, c("gsva", "kras", "emt", "mean"), ext="tsv"))

g <- ggplot(d.sum, aes(x = kras, y = emt, colour = clone)) + theme_clean() +
		geom_vline(xintercept = 0, colour = "grey90") +
		geom_hline(yintercept = 0, colour = "grey90") +
		geom_point(aes(shape = treatment)) + 
		scale_colour_manual(values = clone.cols) +
		scale_shape_manual(values = c(DMSO = 19, Talazoparib = 20)) +
		coord_fixed() +
		xlab("KRAS signaling score") + ylab("EMT pathway score");

add_arrow <- function(g, data, x, xmid, xend, y, ymid, yend, alpha=0.5, tip=0.02) {
	g + geom_segment(
		aes(x = x, xend = xmid, y = y, yend = ymid),
		data = data,
		alpha = alpha,
		arrow = arrow(length = unit(tip, "npc"))
	) +
	geom_segment(
		aes(x = xmid, xend = xend, y = ymid, yend = yend),
		data = data,
		alpha = alpha
	)
}

d.treat <- 
	left_join(
		pivot_wider(d.sum, clone, names_from = treatment, values_from = kras, names_prefix = "es_"),
		pivot_wider(d.sum, clone, names_from = treatment, values_from = emt, names_prefix = "es_"),
		suffix = c("_kras", "_emt"),
		by = "clone"
	);
colnames(d.treat) <- tolower(colnames(d.treat));

d.treat %<>% mutate(
	es_mid_kras = (es_dmso_kras + es_talazoparib_kras) / 2,
	es_mid_emt  = (es_dmso_emt  + es_talazoparib_emt ) / 2
);

qdraw(
	with(d.treat,
		add_arrow(g, d.treat,
			es_dmso_kras, es_mid_kras, es_talazoparib_kras,
			es_dmso_emt,  es_mid_emt,  es_talazoparib_emt
		)
	)
	,
	width = 5, height = 5,
	insert(out.fname, c("kras", "emt", "arrow-treatment"), ext="pdf")
)

d.resist.p <- filter(d.sum, clone == "P") %>% select(-clone);
d.resist <- do.call(rbind, 
	lapply(clones[-1],
		function(cl) {
			left_join(
				d.resist.p,
				filter(d.sum, clone == cl),
				by = "treatment",
				suffix = c("_p", "_r")
			)
		}
	)
);

d.resist %<>% mutate(
	kras_mid = (kras_p + kras_r) / 2,
	emt_mid  = (emt_p  + emt_r ) / 2
);

d.resist.dmso <- filter(d.resist, treatment == "DMSO");

qdraw(
	with(d.resist.dmso,
		add_arrow(g, d.resist.dmso,
			kras_p, kras_mid, kras_r,
			emt_p,  emt_mid,  emt_r
		)
	)
	,
	width = 5, height = 5,
	insert(out.fname, c("kras", "emt", "arrow-resistance", "dmso"), ext="pdf")
)

d.resist.tal <- filter(d.resist, treatment == "Talazoparib");

qdraw(
	with(d.resist.tal,
		add_arrow(g, d.resist.tal,
			kras_p, kras_mid, kras_r,
			emt_p,  emt_mid,  emt_r
		)
	)
	,
	width = 5, height = 5,
	insert(out.fname, c("kras", "emt", "arrow-resistance", "tal"), ext="pdf")
)

d.treat.p <- filter(d.treat, clone == "P");
d.resist.rc1.tal <- filter(d.resist, clone == "RC1", treatment == "Talazoparib");
d.treat.rc1 <- filter(d.treat, clone == "RC1");

g.p <- with(d.treat.p,
	add_arrow(g, d.treat.p, 
		es_dmso_kras, es_mid_kras, es_talazoparib_kras,
		es_dmso_emt,  es_mid_emt,  es_talazoparib_emt
	)
);

g.rc1.tal <- with(d.resist.rc1.tal,
	add_arrow(g.p, d.resist.rc1.tal,
		kras_p, kras_mid, kras_r,
		emt_p,  emt_mid,  emt_r
	)
);

g.rc1.dmso <- with(d.treat.rc1,
	add_arrow(g.rc1.tal, d.treat.rc1,
		es_talazoparib_kras, es_mid_kras, es_dmso_kras,
		es_talazoparib_emt,  es_mid_emt,  es_dmso_emt
	)
);

qdraw(
	g.rc1.dmso,
	width = 5, height = 5,
	insert(out.fname, c("kras", "emt", "arrow-adapt", "rc1"), ext="pdf")
)

