library(io)
library(dplyr)
library(reshape2)
library(limma)
library(parallel)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

source("../R/preamble.R")

# cell line: SUM149
# triple negative, inflammatory breast cancer
# disease progressed through chemotherapy
# mutation: BRCA1 p.P724fs*12
# microamplification in PTEN exon 2 leads to truncated PTEN transcript
# beginning in exon 2; no RNA expression of PTEN for all subsequent exons

# TODO verify PTEN loss of function in RNAseq data

in.fname <- as.filename("parpi-resist_deseq-stat_treatment-clone-interaction_clones-vs-parental.mtx");
#in.fname <- as.filename("parpi-resist_deseq-stat_treatment-clone-interaction_clone_treated-vs-untreated.mtx");


out.fname <- in.fname;
out.fname$ext <- NULL;

mc.cores <- 4;

pheno <- setup_pheno(qread("../sample-info_parpi-resist_stage2.tsv"));
x <- qread(in.fname);

overall <- qread("parpi-resist_resistance.rnk");

#ha <- HeatmapAnnotation(
#	df = select(pheno, clone, treatment, batch, lane, fold_resistance),
#	col = list(
#		treatment = c("DMSO" = "grey30", "None" = "grey60", "Talazoparib" = "royalblue"),
#		clone = clone.cols
#	)
#);

na <- apply(x, 1, function(z) any(is.na(z)));
y <- x[!na, ];

# rename columns to clone names
colnames(y) <- gsub("clone_(C\\d+)_vs_Parental", "\\1", colnames(y));
# order columns
y <- y[, clones[-1]];


gsets.h <- read_msigdb("h.all");

#d <- data.frame(comparison = c(colnames(x), "Parental"));
#d$reference <- 0;
#d$reference[d$comparison != "Parental"] <- 1;

#design <- model.matrix(~ reference, data=d);

#gset <- gsets.h$data$HALLMARK_KRAS_SIGNALING_UP;

camera_single <- function(s, gsets=NULL, index=NULL) {
	if (is.null(index)) {
		index <- lapply(gsets,
			function(gset) {
				names(s) %in% gset
			}
		);
	}

	d <- cameraPR(s, index, sort=FALSE);
	d$delta <- unlist(lapply(index,
		function(i) {
			mean(s[i]) - mean(s[!i])
		}
	))
	d$z1 <- unlist(lapply(index,
		function(i) {
			mean(s[i])
		}
	));

	d
}

camera_batch <- function(y, gsets) {
	index <- lapply(gsets,
		function(gset) {
			rownames(y) %in% gset
		}
	);

	res <- lapply(1:ncol(y), function(j) {
		camera_single(y[,j], index=index)
	});
	names(res) <- colnames(y);

	res
}

camera_transform <- function(y, gsets) {
	index <- lapply(gsets,
		function(gset) {
			rownames(y) %in% gset
		}
	);

	res <- apply(y, 2, function(s) {
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


theme_clean <- function() {
	theme_bw() +
		theme(
			strip.background = element_blank(),
			panel.grid = element_blank()
		)
}

cam_volcano_plot <- function(cam, fdr.cut=0.01, z.cut=0.5, label.size=4) {
	if (is.null(cam$gset)) {
		cam$gset = rownames(cam);
	}

	cam <- mutate(cam,
		keep = FDR < fdr.cut & abs(z1) > z.cut,
		group = factor(
			ifelse(keep, ifelse(delta < 0, "Down", "Up"), "NS"),
			c("Down", "Up", "NS")
		),
		gset = rename_hallmarks(gset)
	);

	ggplot(cam, aes(x=z1, y=FDR, alpha=keep, colour=group, label=ifelse(keep, gset, NA))) +
		theme_clean() +
		geom_vline(xintercept=0) +
		geom_vline(xintercept=c(z.cut, -z.cut), linetype=3, colour="grey60") +
		geom_hline(yintercept=fdr.cut, linetype=3, colour="grey60") +
		geom_point(show.legend=FALSE) +
		geom_label_repel(show.legend=FALSE, nudge_y=0.1, nudge_x=0.1, size=label.size) +
		scale_y_continuous(trans=revlog_trans(10), sec.axis = dup_axis(name=NULL)) +
		scale_colour_manual(values=diff.group.cols) +
		xlab("mean expression difference") + ylab("false discovery rate") +
		#annotate("text", label=levels(cam$group)[1:2], x=c(-1.5, 1.5), y=0.2) +
		xlim(-2, 2)
}


clones.cams.h <- camera_batch(y, gsets.h$data);
overall.cam.h <- camera_single(overall, gsets.h$data);

for (i in 1:length(clones.cams.h)) {
	name <- names(clones.cams.h)[i];	
	qdraw(
		cam_volcano_plot(clones.cams.h[[i]]) + ggtitle(paste0(name, " vs. Parental"))
		,
		width = 8, height = 5,
		file = insert(out.fname, c("camera", "volcano", tolower(name)), ext="pdf")
	)
}

h.remove <- c(
	"HALLMARK_ANGIOGENESIS",
	"HALLMARK_INFLAMMATORY_RESONSE",
	"HALLMARK_ALLOGRAFT_REJECTION",
	"HALLMARK_COMPLEMENT",
	"HALLMARK_UV_RESPONSE_UP"
)

c2.cam.h <- clones.cams.h$C2;
c2.cam.h$gset <- rownames(c2.cam.h);
c2.cam.h$gset[c2.cam.h$gset %in% h.remove] <- NA;

qdraw(
	cam_volcano_plot(c2.cam.h) + ggtitle("C2 vs. Parental")
	,
	width = 8, height = 5,
	file = insert(out.fname, c("camera", "volcano", "c2", "sel"), ext="pdf")
)


qdraw(
	cam_volcano_plot(overall.cam.h) + ggtitle("All resistant clones vs. Parental") +
		coord_cartesian(ylim=c(max(overall.cam.h$FDR), 1e-9))
		#annotation_logticks(side="lr")
	,
	width = 8, height = 5,
	file = insert(out.fname, c("camera", "volcano", "all"), ext="pdf")
)

cams.df <- do.call(rbind,
	mapply(function(d, name) data.frame(comparison=name, gset=rownames(d), d),
		c(list(overall.cam.h), clones.cams.h),
		paste0(c("All resistant clones", names(clones.cams.h)), " vs. Parental"),
		SIMPLIFY=FALSE
	)
);
rownames(cams.df) <- NULL;

qdraw(
	cam_volcano_plot(cams.df, label.size=2.5) + facet_wrap(~ comparison, ncol=2)
	,
	width = 8, height = 11,
	file = insert(out.fname, c("camera", "volcano", "each"), ext="pdf")
)

####

es.h <- camera_transform(y, gsets.h$data);

summary(es.h)

colf <- circlize::colorRamp2(c(-10, 0, 10), c("blue", "white", "red"));


pdf(tag(out.fname, c("camera", "h"), ext="pdf"), width=10, height=10);
Heatmap(es.h, col=colf, cluster_columns = FALSE)
dev.off();
#Heatmap(es.h, top_annotation = ha, cluster_columns = TRUE)

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

Heatmap(es.c6, col=colf, cluster_columns = FALSE, row_names_gp = gpar(fontsize=4))

dim(es.c6)
#hist(es.c6)
summary(es.c6)

means <- rowMeans(es.c6);
summary(means)
n.sub <- 80;

# upregulated in resistant clones
#es.c6.up <- es.c6[means > 0, , drop=FALSE]
es.c6.up <- es.c6[means > 1, , drop=FALSE]

ha <- NULL;

pdf(tag(out.fname, c("camera", "c6", "parental-up"), ext="pdf"), width=10, height=15);
Heatmap(es.c6.up, col=colf, top_annotation = ha, cluster_columns = FALSE, cluster_rows = FALSE, row_names_gp = gpar(fontsize=8));
dev.off();

# downregulated in resistant clones
#es.c6.sub <- es.c6[order(means)[1:n.sub], , drop=FALSE]
es.c6.sub <- es.c6[means < -2, , drop=FALSE]
#es.c6.down <- es.c6[means < -1, , drop=FALSE]

pdf(tag(out.fname, c("camera", "c6", "parental-down"), ext="pdf"), width=10, height=15);
Heatmap(es.c6.down, col=colf, top_annotation = ha, cluster_columns = FALSE, cluster_rows = FALSE, row_names_gp = gpar(fontsize=8))
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
Heatmap(es.c7.sub, col=colf, top_annotation = ha, cluster_columns = FALSE, cluster_rows = FALSE, row_names_gp = gpar(fontsize=8));

es.c7.sub <- es.c7[order(means)[1:n.sub], , drop=FALSE]
Heatmap(es.c7.sub, col=colf, top_annotation = ha, cluster_columns = FALSE, cluster_rows = FALSE, row_names_gp = gpar(fontsize=8));

#

gsets.c3 <- read_msigdb("c3.tft");

es.c3 <- camera_transform(y, gsets.c3$data);

summary(es.c3)

means <- rowMeans(es.c3);

es.c3.sub <- es.c3[order(-means)[1:n.sub], , drop=FALSE]
Heatmap(es.c3.sub, col=colf, top_annotation = ha, cluster_columns = FALSE, cluster_rows = FALSE, row_names_gp = gpar(fontsize=8));

es.c3.sub <- es.c3[order(means)[1:n.sub], , drop=FALSE]
#es.c3.sub <- es.c3[means < -1, , drop=FALSE]
Heatmap(es.c3.sub, col=colf, top_annotation = ha, cluster_columns = FALSE, cluster_rows = FALSE, row_names_gp = gpar(fontsize=8));

#

gsets.c2 <- read_msigdb("c2.cgp");

es.c2 <- camera_transform(y, gsets.c2$data);

summary(es.c2)

means <- rowMeans(es.c2);

#es.c2.sub <- es.c2[order(-means)[1:n.sub], , drop=FALSE]
es.c2.sub <- es.c2[means > 2, , drop=FALSE]
Heatmap(es.c2.sub, col=colf, top_annotation = ha, cluster_columns = FALSE, cluster_rows = FALSE, row_names_gp = gpar(fontsize=8));

#es.c2.sub <- es.c2[order(means)[1:n.sub], , drop=FALSE]
es.c2.sub <- es.c2[means < -2, , drop=FALSE]
Heatmap(es.c2.sub, col=colf, top_annotation = ha, cluster_columns = FALSE, cluster_rows = FALSE, row_names_gp = gpar(fontsize=8));

#

gsets.c2 <- read_msigdb("c2.cp");

es.c2 <- camera_transform(y, gsets.c2$data);

summary(es.c2)

means <- rowMeans(es.c2);
#maxes <- apply(es.c2, 1, function(z) max(abs(z)));

#es.c2.sub <- es.c2[order(-means)[1:n.sub], , drop=FALSE]
es.c2.sub <- es.c2[means > 2, , drop=FALSE]
Heatmap(es.c2.sub, col=colf, top_annotation = ha, cluster_columns = FALSE, cluster_rows = FALSE, row_names_gp = gpar(fontsize=8));

#es.c2.sub <- es.c2[order(means)[1:n.sub], , drop=FALSE]
es.c2.sub <- es.c2[means < -2, , drop=FALSE]
Heatmap(es.c2.sub, col=colf, top_annotation = ha, cluster_columns = FALSE, cluster_rows = FALSE, row_names_gp = gpar(fontsize=8));

