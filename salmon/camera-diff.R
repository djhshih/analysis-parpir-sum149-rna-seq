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

# remove genes with NA statistic
na <- apply(x, 1, function(z) any(is.na(z)));
y <- x[!na, ];

# rename columns to clone names
colnames(y) <- rename_clones(gsub("clone_(C\\d+)_vs_Parental", "\\1", colnames(y)));
# order columns
y <- y[, clones[-1]];


####

gsets.h <- read_msigdb("h.all");

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

cam_volcano_plot <- function(cam, fdr.cut=0.01, z.cut=0.5, label.size=4, rename_gset=identity) {
	if (is.null(cam$gset)) {
		cam$gset = rownames(cam);
	}

	cam <- mutate(cam,
		keep = FDR < fdr.cut & abs(z1) > z.cut,
		group = factor(
			ifelse(keep, ifelse(delta < 0, "Down", "Up"), "NS"),
			c("Down", "Up", "NS")
		),
		gset = rename_gset(gset)
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

overall.cam.h <- camera_single(overall, gsets.h$data);

qdraw(
	cam_volcano_plot(overall.cam.h, rename_gset=rename_hallmarks) + ggtitle("RCs vs. P") +
		coord_cartesian(ylim=c(max(overall.cam.h$FDR), 1e-9))
		#annotation_logticks(side="lr")
	,
	width = 8, height = 5,
	file = insert(out.fname, c("camera", "volcano", "all"), ext="pdf")
)

clones.cams.h <- camera_batch(y, gsets.h$data);

for (i in 1:length(clones.cams.h)) {
	name <- names(clones.cams.h)[i];	
	qdraw(
		cam_volcano_plot(clones.cams.h[[i]], rename_gset=rename_hallmarks) + ggtitle(paste0(name, " vs. P"))
		,
		width = 8, height = 5,
		file = insert(out.fname, c("camera", "volcano", tolower(name)), ext="pdf")
	)
}

h.omit <- c(
	"HALLMARK_ANGIOGENESIS",
	"HALLMARK_INFLAMMATORY_RESONSE",
	"HALLMARK_ALLOGRAFT_REJECTION",
	"HALLMARK_COMPLEMENT",
	"HALLMARK_UV_RESPONSE_UP"
)

rc1.cam.h <- clones.cams.h$RC1;
rc1.cam.h$gset <- rownames(rc1.cam.h);
rc1.cam.h$gset[rc1.cam.h %in% h.omit] <- NA;

qdraw(
	cam_volcano_plot(rc1.cam.h, rename_gset=rename_hallmarks) + ggtitle("RC1 vs. P")
	,
	width = 8, height = 5,
	file = insert(out.fname, c("camera", "volcano", "rc1", "sel"), ext="pdf")
)

cams.df <- do.call(rbind,
	mapply(function(d, name) data.frame(comparison=name, gset=rownames(d), d),
		c(list(overall.cam.h), clones.cams.h),
		paste0(c("RCs", names(clones.cams.h)), " vs. P"),
		SIMPLIFY=FALSE
	)
);
rownames(cams.df) <- NULL;

qdraw(
	cam_volcano_plot(cams.df, rename_gset=rename_hallmarks, label.size=2.5) + 
		facet_wrap(~ comparison, ncol=2)
	,
	width = 8, height = 11,
	file = insert(out.fname, c("camera", "volcano", "each"), ext="pdf")
)

####

plot_gene_set_density <- function(y, genes=NULL) {
	if (is.null(genes)) {
		y.sel <- y;
	} else {
		y.sel <- y[rownames(y) %in% genes, ];
	}
	y.sel.m <- melt(y.sel, varnames=c("gene", "comparison"));

	ggplot(y.sel.m, aes(x=value, fill=comparison)) + theme_clean() +
		geom_vline(xintercept=0, colour="grey60", linetype=2) +
		geom_density(alpha=0.5) +
		scale_fill_manual(values=clone.cols) +
		facet_wrap(~ comparison, ncol=1) +
		guides(fill=FALSE) +
		xlab("standardized difference vs. parental") +
		xlim(-10, 10)
}

plot.opts <- getOption("plot");
options(plot = within(plot.opts, {width <- 3; height <- 6;}));

qdraw(
	plot_gene_set_density(y) + ggtitle("All genes")
	,
	file = insert(out.fname, c("gene-set-density", "all"), ext="pdf")
)

qdraw(
	plot_gene_set_density(y, gsets.h$data$HALLMARK_OXIDATIVE_PHOSPHORYLATION) +
		ggtitle("Oxidative phosphorylation")
	,
	file = insert(out.fname, c("gene-set-density", "oxphos"), ext="pdf")
)

qdraw(
	plot_gene_set_density(y, gsets.h$data$HALLMARK_KRAS_SIGNALING_UP) +
		ggtitle("KRAS signaling")
	,
	file = insert(out.fname, c("gene-set-density", "nfkb"), ext="pdf")
)

qdraw(
	plot_gene_set_density(y, gsets.h$data$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION) +
		ggtitle("EMT pathway")
	,
	file = insert(out.fname, c("gene-set-density", "emt"), ext="pdf")
)

qdraw(
	plot_gene_set_density(y, gsets.h$data$HALLMARK_TNFA_SIGNALING_VIA_NFK) +
		ggtitle("NF-κB signaling")
	,
	file = insert(out.fname, c("gene-set-density", "nfkb"), ext="pdf")
)

qdraw(
	plot_gene_set_density(y, gsets.h$data$HALLMARK_DNA_REPAIR) +
		ggtitle("DNA repair")
	,
	file = insert(out.fname, c("gene-set-density", "dna-repair"), ext="pdf")
)

# restore default plot options
options(plot = plot.opts);

####

gsets.c6 <- read_msigdb("c6.all");

overall.cam.c6 <- camera_single(overall, gsets.c6$data);

qdraw(
	cam_volcano_plot(overall.cam.c6, label.size=2.5) + ggtitle("RCs vs. P")
	,
	width = 8, height = 5,
	file = insert(out.fname, c("camera", "c6", "volcano", "all"), ext="pdf")
)

clones.cams.c6 <- camera_batch(y, gsets.c6$data);

for (i in 1:length(clones.cams.c6)) {
	name <- names(clones.cams.c6)[i];	
	qdraw(
		cam_volcano_plot(clones.cams.c6[[i]], label.size=2.5) + ggtitle(paste0(name, " vs. P"))
		,
		width = 8, height = 5,
		file = insert(out.fname, c("camera", "c6", "volcano", tolower(name)), ext="pdf")
	)
}

rc1.cam.c6 <- clones.cams.c6$RC1;
rc1.cam.c6$gset <- rownames(rc1.cam.c6);
rc1.cam.c6$gset[-grep("KRAS", rownames(rc1.cam.c6))] <- NA;

qdraw(
	cam_volcano_plot(rc1.cam.c6, label.size=2.5) + ggtitle("RC1 vs. P")
	,
	width = 8, height = 5,
	file = insert(out.fname, c("camera", "c6", "volcano", "rc1", "sel"), ext="pdf")
)

# Results summary

# resistance clones:
# - downregulate EGFR signaling
# - downregulating NFkappaB signaling
# - suppress response to BMI1 downregulation
# - downregulating KRAS signaling
# - suppress gene expression associated with KRAS dependency
# - suppress signaling in response to STK33 knockdown
# - suppress IL2, IL15 signaling

####

gsets.c7 <- read_msigdb("c7.all");

overall.cam.c7 <- camera_single(overall, gsets.c7$data);

qdraw(
	cam_volcano_plot(overall.cam.c7, label.size=2.5) + ggtitle("RCs vs. P") +
		coord_cartesian(ylim=c(max(overall.cam.c7$FDR), 1e-5))
	,
	width = 8, height = 5,
	file = insert(out.fname, c("camera", "c7", "volcano", "all"), ext="pdf")
)

####

gsets.c3 <- read_msigdb("c3.tft");

overall.cam.c3 <- camera_single(overall, gsets.c3$data);

# Enriched KRCTCNNNNMANAGC and TTTNNANAGCYR motif may be bound by ILF3 (DRBP76)
# https://www.ncbi.nlm.nih.gov/pubmed/20668518
# ILF3 is involved in interleukin signaling

qdraw(
	cam_volcano_plot(overall.cam.c3, label.size=2.5) + ggtitle("RCs vs. P")
	,
	width = 8, height = 5,
	file = insert(out.fname, c("camera", "c3", "volcano", "all"), ext="pdf")
)

####

gsets.c2.cgp <- read_msigdb("c2.cgp");

overall.cam.c2.cgp <- camera_single(overall, gsets.c2.cgp$data);

qdraw(
	cam_volcano_plot(overall.cam.c2.cgp, label.size=2.5) + ggtitle("RCs vs. P") +
		coord_cartesian(ylim=c(max(overall.cam.c2.cgp$FDR), 1e-9))
	,
	width = 8, height = 10,
	file = insert(out.fname, c("camera", "c2-cgp", "volcano", "all"), ext="pdf")
)

####

gsets.c2.cp <- read_msigdb("c2.cp");

overall.cam.c2.cp <- camera_single(overall, gsets.c2.cp$data);

qdraw(
	cam_volcano_plot(overall.cam.c2.cp, label.size=2.5) + ggtitle("RCs vs. P") +
		coord_cartesian(ylim=c(max(overall.cam.c2.cp$FDR), 1e-9))
	,
	width = 8, height = 10,
	file = insert(out.fname, c("camera", "c2-cp", "volcano", "all"), ext="pdf")
)

graphics.off()

