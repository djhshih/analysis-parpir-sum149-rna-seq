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
source("../R/camera.R")

# cell line: SUM149
# triple negative, inflammatory breast cancer
# disease progressed through chemotherapy
# mutation: BRCA1 p.P724fs*12
# microamplification in PTEN exon 2 leads to truncated PTEN transcript
# beginning in exon 2; no RNA expression of PTEN for all subsequent exons

# TODO verify PTEN loss of function in RNAseq data

overall <- qread("parpi-resist_resistance.rnk");
in.fname <- as.filename("parpi-resist_deseq-stat_treatment-clone-interaction_clones-vs-parental.mtx");

out.fname <- in.fname;
out.fname$ext <- NULL;

pheno <- setup_pheno(qread("../sample-info_parpi-resist_stage2.tsv"), rename.clones=TRUE);
x <- qread(in.fname);


# remove genes with NA statistic
na <- apply(x, 1, function(z) any(is.na(z)));
y <- x[!na, ];

# rename columns to clone names
colnames(y) <- rename_clones(gsub(".*clone_?(C\\d+).*", "\\1", colnames(y)));
# order columns
y <- y[, clones[-1]];


####

gsets.h <- read_msigdb("h.all");

overall.cam.h <- camera_single(overall, gsets.h$data);

qdraw(
	cam_volcano_plot(overall.cam.h, rename_gset=rename_hallmarks) + ggtitle("RCs vs. P") +
		coord_cartesian(ylim=c(max(overall.cam.h$FDR), 1e-9)) + xlim(-2, 2)
		#annotation_logticks(side="lr")
	,
	width = 8, height = 5,
	file = insert(out.fname, c("camera", "volcano", "all"), ext="pdf")
)

clones.cams.h <- camera_batch(y, gsets.h$data);

for (i in 1:length(clones.cams.h)) {
	name <- names(clones.cams.h)[i];	
	qdraw(
		cam_volcano_plot(clones.cams.h[[i]], rename_gset=rename_hallmarks) + 
			ggtitle(paste0(name, " vs. P")) + xlim(-2, 2)
		,
		width = 8, height = 5,
		file = insert(out.fname, c("camera", "volcano", tolower(name)), ext="pdf")
	)
}

h.omit <- c(
	"HALLMARK_ANGIOGENESIS",
	"HALLMARK_INFLAMMATORY_RESPONSE",
	"HALLMARK_ALLOGRAFT_REJECTION",
	"HALLMARK_COMPLEMENT",
	"HALLMARK_UV_RESPONSE_UP",
	"HALLMARK_INTERFERON_ALPHA_RESPONSE",
	"HALLMARK_INTERFERON_GAMMA_RESPONSE",
	#"HALLMARK_TNFA_SIGNALING_VIA_NFKB",
	"HALLMARK_IL2_STAT5_SIGNALING",
	"HALLMARK_IL6_JAK_STAT3_SIGNALING",
	"HALLMARK_COAGULATION"
);

rc1.cam.h <- clones.cams.h$RC1;
rc1.cam.h$gset <- rownames(rc1.cam.h);
rc1.cam.h$gset[rc1.cam.h$gset %in% h.omit] <- NA;

qdraw(
	cam_volcano_plot(rc1.cam.h, rename_gset=rename_hallmarks) + 
		ggtitle("RC1 vs. P") + xlim(-2, 2)
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
		facet_wrap(~ comparison, ncol=2) + xlim(-2, 2)
	,
	width = 8, height = 11,
	file = insert(out.fname, c("camera", "volcano", "each"), ext="pdf")
)

####

plot.opts <- getOption("plot");
options(plot = within(plot.opts, {width <- 2; height <- 8;}));

qdraw(
	plot_gene_set_density(y) + ggtitle("All genes") + xlim(-6, 6)
	,
	file = insert(out.fname, c("gene-set-density", "all"), ext="pdf")
)

qdraw(
	plot_gene_set_density(y, gsets.h$data$HALLMARK_OXIDATIVE_PHOSPHORYLATION) +
		ggtitle("Oxphos pathway") + xlim(-6, 6)
	,
	file = insert(out.fname, c("gene-set-density", "oxphos"), ext="pdf")
)

qdraw(
	plot_gene_set_density(y, gsets.h$data$HALLMARK_KRAS_SIGNALING_UP) +
		ggtitle("KRAS signaling") + xlim(-6, 6)
	,
	file = insert(out.fname, c("gene-set-density", "kras"), ext="pdf")
)

qdraw(
	plot_gene_set_density(y, gsets.h$data$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION) +
		ggtitle("EMT pathway") + xlim(-6, 6)
	,
	file = insert(out.fname, c("gene-set-density", "emt"), ext="pdf")
)

qdraw(
	plot_gene_set_density(y, gsets.h$data$HALLMARK_TNFA_SIGNALING_VIA_NFK) +
		ggtitle("NF-κB signaling") + xlim(-6, 6)
	,
	file = insert(out.fname, c("gene-set-density", "nfkb"), ext="pdf")
)

qdraw(
	plot_gene_set_density(y, gsets.h$data$HALLMARK_DNA_REPAIR) +
		ggtitle("DNA repair") + xlim(-6, 6)
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

