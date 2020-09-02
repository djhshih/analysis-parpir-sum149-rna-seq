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

overall <- qread("parpi-resist_treatment-clone_treatment.rnk");
parental <- qread("parpi-resist_treatment-clone-interaction_treatment.rnk");
in.fname <- as.filename("parpi-resist_deseq-stat_treatment-clone-interaction_clone_treated-vs-untreated.mtx");

out.fname <- in.fname;
out.fname$ext <- NULL;

mc.cores <- 4;
fdr.cut <- 0.05;

pheno <- setup_pheno(qread("../sample-info_parpi-resist_stage2.tsv"), rename.clones=TRUE);
x <- qread(in.fname);


# remove genes with NA statistic
na <- apply(x, 1, function(z) any(is.na(z)));
y <- x[!na, ];

# rename columns to clone names
colnames(y) <- rename_clones(gsub(".*clone_?(C\\d+).*", "\\1", colnames(y)));
# order columns
y <- y[, clones[-1]];

parental <- parental[match(rownames(y), names(parental))];
y <- cbind(P = parental, y);

####

gsets.h <- read_msigdb("h.all");

overall.cam.h <- camera_single(overall, gsets.h$data);

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

overall.cam.h$gset <- rownames(overall.cam.h);
overall.cam.h$gset[overall.cam.h$gset %in% h.omit] <- NA;

qdraw(
	cam_volcano_plot(overall.cam.h, rename_gset=rename_hallmarks, fdr.cut=fdr.cut) + ggtitle("Talazoparib vs. DMSO") +
		xlim(-2, 2)
	,
	width = 8, height = 5,
	file = insert(out.fname, c("camera", "volcano", "all"), ext="pdf")
)

overall.cam.h$gset <- NULL;

clones.cams.h <- camera_batch(y, gsets.h$data);

for (i in 1:length(clones.cams.h)) {
	name <- names(clones.cams.h)[i];	
	qdraw(
		cam_volcano_plot(clones.cams.h[[i]], rename_gset=rename_hallmarks, fdr.cut=fdr.cut) + 
			ggtitle(paste0("Talazoparib vs. DMSO in ", name)) + xlim(-2, 2)
		,
		width = 8, height = 5,
		file = insert(out.fname, c("camera", "volcano", tolower(name)), ext="pdf")
	)
}

cams.df <- do.call(rbind,
	mapply(
		function(d, name) data.frame(comparison=name, gset=rownames(d), d),
		clones.cams.h,
		paste0("Talazoparib vs. DMSO in ", names(clones.cams.h)),
		SIMPLIFY=FALSE
	)
);
rownames(cams.df) <- NULL;

qdraw(
	cam_volcano_plot(cams.df, rename_gset=rename_hallmarks, label.size=2.5, fdr.cut=fdr.cut) + 
		facet_wrap(~ comparison, ncol=2) + xlim(-2, 2)
	,
	width = 8, height = 11,
	file = insert(out.fname, c("camera", "volcano", "each"), ext="pdf")
)

####

plot.opts <- getOption("plot");
options(plot = within(plot.opts, {width <- 2; height <- 6;}));

qdraw(
	plot_gene_set_density(y) + ggtitle("All genes") + xlim(-6, 6)
	,
	file = insert(out.fname, c("gene-set-density", "all"), ext="pdf")
)

qdraw(
	plot_gene_set_density(y, gsets.h$data$HALLMARK_E2F_TARGETS) +
		ggtitle("E2F targets") + xlim(-6, 6)
	,
	file = insert(out.fname, c("gene-set-density", "e2f"), ext="pdf")
)

qdraw(
	plot_gene_set_density(y, gsets.h$data$HALLMARK_G2M_CHECKPOINT) +
		ggtitle("G2M checkpoint") + xlim(-6, 6)
	,
	file = insert(out.fname, c("gene-set-density", "g2m-checkpoint"), ext="pdf")
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
	cam_volcano_plot(overall.cam.c6, label.size=2.5, fdr.cut=fdr.cut) + ggtitle("RCs vs. P")
	,
	width = 8, height = 5,
	file = insert(out.fname, c("camera", "c6", "volcano", "all"), ext="pdf")
)

clones.cams.c6 <- camera_batch(y, gsets.c6$data);

for (i in 1:length(clones.cams.c6)) {
	name <- names(clones.cams.c6)[i];	
	qdraw(
		cam_volcano_plot(clones.cams.c6[[i]], label.size=2.5, fdr.cut=fdr.cut) + 
			ggtitle(paste0(name, " vs. P"))
		,
		width = 8, height = 5,
		file = insert(out.fname, c("camera", "c6", "volcano", tolower(name)), ext="pdf")
	)
}

####

# too many to visualize!

#gsets.c7 <- read_msigdb("c7.all");
#overall.cam.c7 <- camera_single(overall, gsets.c7$data);

####

gsets.c3 <- read_msigdb("c3.tft");

overall.cam.c3 <- camera_single(overall, gsets.c3$data);

# Enriched KRCTCNNNNMANAGC and TTTNNANAGCYR motif may be bound by ILF3 (DRBP76)
# https://www.ncbi.nlm.nih.gov/pubmed/20668518
# ILF3 is involved in interleukin signaling

qdraw(
	cam_volcano_plot(overall.cam.c3, label.size=2.5, fdr.cut=fdr.cut) + ggtitle("RCs vs. P") +
		coord_cartesian(ylim=c(max(overall.cam.c3$FDR), 1e-4))
	,
	width = 8, height = 5,
	file = insert(out.fname, c("camera", "c3", "volcano", "all"), ext="pdf")
)

####

# too many to visualize!

#gsets.c2.cgp <- read_msigdb("c2.cgp");

#overall.cam.c2.cgp <- camera_single(overall, gsets.c2.cgp$data);

####

# too many to visualize!

#gsets.c2.cp <- read_msigdb("c2.cp");

#overall.cam.c2.cp <- camera_single(overall, gsets.c2.cp$data);

graphics.off()

