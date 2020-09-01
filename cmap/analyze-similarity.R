library(cmapR)
library(mmalign)
library(dplyr)
library(limma)
library(ComplexHeatmap)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(io)

source("../R/preamble.R")
source("../R/camera.R")

lincs.indir <- "~/data/cmap/phase1";

source(file.path(lincs.indir, "common.R"));

# use only landmark genes
gct <- file.path(lincs.indir, "lincs-sig_lm_n473647x978.gctx");
pheno <- qread(file.path(lincs.indir, "lincs-sig-info.rds"));
cell <- qread(file.path(lincs.indir, "lincs-cell-info.rds"));
fannot <- qread(file.path(lincs.indir, "lincs-gene-info_lm.rds"));
groups <- qread(file.path(lincs.indir, "lincs-sig-groups.rds"));

compound.annot <- qread("compound-annot.tsv");

#in.fname <- as.filename("../salmon/parpi-resist_deseq-stat_treatment-clone-interaction_clones-vs-parental.mtx");
in.fname <- as.filename("../salmon/parpi-resist_deseq-stat_treatment-clone-interaction_treatment_clones-vs-parental.mtx");

fname <- in.fname;
fname$path <- NULL;
fname$ext <- NULL;

ss.fname <- insert(fname, "lincs-lm-sim", ext="rds");
ss.null.fname <- insert(fname, "lincs-lm-sim-null", ext="rds");

x <- qread(in.fname);

ss <- qread(ss.fname);
ss.null <- qread(ss.null.fname);

group_name_to_clone <- function(x) {
	as.character(factor(sub(".*clone", "", x), levels=clones.from, labels=clones))
}

rownames(ss) <- group_name_to_clone(rownames(ss));
rownames(ss.null) <- group_name_to_clone(rownames(ss.null));

out.fname <- fname;


# Reformat data ##############################################################

# rearrange fannot to match data

rid <- as.integer(read.gctx.ids(gct));
cid <- read.gctx.ids(gct, dimension="column");

stopifnot( !any(is.na(rid)) )
stopifnot( all(fannot$pr_gene_id == rid) )
stopifnot( all(pheno$sig_id == cid) )

# rearrange x to match reference data

x.idx <- match(fannot$pr_gene_symbol, rownames(x));
table(is.na(x.idx))

# NB x is populated with NAs in order to match the reference data
x <- x[x.idx, ];

colnames(x) <- group_name_to_clone(colnames(x));


##############################################################################

mean_transform <- function(ss, group) {
	res <- apply(ss, 1,
		function(ss1) {
			unlist(lapply(group,
				function(g) {
					mean(ss1[g])
				}
			))
		}
	);
	rownames(res) <- names(group);

	res
}

max_transform <- function(ss, group) {
	res <- apply(ss, 1,
		function(ss1) {
			unlist(lapply(group,
				function(g) {
					max(ss1[g])
				}
			))
		}
	);
	rownames(res) <- names(group);

	res
}

min_transform <- function(ss, group) {
	res <- apply(ss, 1,
		function(ss1) {
			unlist(lapply(group,
				function(g) {
					min(ss1[g])
				}
			))
		}
	);
	rownames(res) <- names(group);

	res
}

median_percentile_transform <- function(ss, group) {
	res <- apply(ss, 1,
		function(ss1) {
			z <- rank(ss1) / length(ss1) * 100;
			unlist(lapply(group,
				function(g) {
					# percentile
					median(z[g])
				}
			))
		}
	);
	rownames(res) <- names(group);

	res
}

camera_fdr_transform <- function(ss, group) {
	res <- apply(ss, 1,
		function(ss1) {
			d <- cameraPR(ss1, group, sort=FALSE);
			z <- log10(d$FDR);

			# flip the sign of enriched results
			idx <- d$Direction == "Up";
			z[idx] <- -z[idx];

			z
		}
	);
	rownames(res) <- names(group);

	res
}

test_vs_null <- function(ss1, ss.null1, alternative = c("two.sided", "less", "greater")) {
	mu <- mean(ss.null1);
	sigma <- sd(ss.null1);

	alternative <- match.arg(alternative);

	if (alternative == "less") {
		pnorm(ss1, mean=mu, sd=sigma, lower.tail=TRUE)
	} else if (alternative == "greater") {
		pnorm(ss1, mean=mu, sd=sigma, lower.tail=FALSE)
	} else {
		pmin(1, pnorm(-abs(ss1), mean=mu, sd=sigma, lower.tail=TRUE) * 2)
	}
}

simplify_features <- function(features) {
	sub("^.+,", "", features)
}

####

qdraw(
	{
		plot(density(ss.null), las=1, xlab="Cosine similarity", main=NA, col="grey30", lwd=2, bty="n")
		lines(density(ss), las=1, xlab="Cosine similarity", main=NA, col="firebrick", lwd=2)
	},
	width = 4, height = 3,
	file = tag(out.fname, c("cos-sim", "density"), ext="pdf")
);

####

ss.p <- do.call(rbind, lapply(
	rownames(ss),
	function(i) {
		test_vs_null(ss[i, ], ss.null[i, ], alternative="greater")
	}
));
rownames(ss.p) <- rownames(ss);
colnames(ss.p) <- colnames(ss);

hist(ss.p, breaks=100)

ss.q <- t(apply(ss.p, 1, p.adjust, method="BH"));

hist(ss.q, breaks=100)

ss.q.min <- min_transform(ss.q, groups$trt_cp);

ss.lq <- log10(ss.q.min);

idx <- apply(ss.lq, 1, function(lq) any(lq < log10(0.05)));

ss.lq[idx, ]

idx <- order(apply(ss.lq, 1, mean));

ss.lq[idx[1:100], ]

####

ms.cp <- mean_transform(ss, groups$trt_cp);
es.cp <- camera_fdr_transform(ss, groups$trt_cp);

cams <- lapply(colnames(es.cp),
	function(clone) {
	data.frame(
		gset = simplify_features(rownames(es.cp)),
		z1 = ms.cp[, clone],
		FDR = 10^(-abs(es.cp[, clone])),
		stringsAsFactors=FALSE
	)
});
names(cams) <- colnames(es.cp);

for (clone in names(cams)) {
	d <- cams[[clone]];
	d.a <- left_join(d, rename(compound.annot, gset=pert_iname, group=compound_class));
	d.a$group <- factor(d.a$group, levels=c("MEK inhibitor", "ERK inhibitor", "RAF inhibitor", "PI3K inhibitor"));

	qdraw(
		cam_volcano_plot(filter(d.a, z1 > 0), fdr.cut=1e-06, z.cut = 0.05, nudge=0, down.sample=0.97, show.legend=TRUE, two.sided=FALSE, ylim=c(1, 1e-30)) + 
			guides(alpha=FALSE, colour="legend") + labs(colour="compound class") +
			xlim(0, 0.15) +
			xlab("mean cosine similarity") +
			scale_colour_npg(na.value="grey40") +
			ggtitle(clone)
		,
		width = 8, height = 5,
		file = insert(out.fname, c("cmap-cosine", "camera", "volcano", clone), ext="pdf")
	);
}

cam.all <- do.call(rbind, 
	mapply(
		function(d, clone) {
			data.frame(d, clone=clone)
		},
		cams,
		names(cams),
		SIMPLIFY=FALSE
	)
);
rownames(cam.all) <- NULL;
cam.all <- rename(cam.all, clone = clone, compound = gset, similarity = z1, fdr = FDR);

qwrite(cam.all, insert(out.fname, c("cmap-cosine", "camera"), ext="tsv"));

####

ss.m <- matrix(colMeans(ss), nrow=1);
colnames(ss.m) <- colnames(ss);

m.es.cp <- camera_fdr_transform(ss.m, groups$trt_cp);
m.ms.cp <- mean_transform(ss.m, groups$trt_cp);

m.cam <- data.frame(
	gset = simplify_features(rownames(m.es.cp)),
	z1 = m.ms.cp[, 1],
	FDR = 10^(-abs(m.es.cp[, 1])),
	stringsAsFactors=FALSE
);
rownames(m.cam) <- NULL;

m.cam$group <- m.cam$gset %in% c("vinblastine", "vincristine");

qdraw(
	cam_volcano_plot(m.cam, fdr.cut=1e-06, z.cut = 0.025, nudge=-0.05, down.sample=0.97, two.sided=FALSE) + 
		xlim(-0.2, 0) +
		xlab("mean cosine similarity") +
		scale_color_manual(values=c("grey30", "#2166AC")) +
		ggtitle("RCs")
	,
	width = 8, height = 5,
	file = insert(out.fname, c("cmap-cosine", "camera", "volcano", "all", "down"), ext="pdf")
);

m.cam.a <- left_join(m.cam, rename(compound.annot, gset=pert_iname, group=compound_class));
m.cam.a$group <- factor(m.cam.a$group, levels=c("MEK inhibitor", "ERK inhibitor", "RAF inhibitor", "PI3K inhibitor"));

qdraw(
	cam_volcano_plot(filter(m.cam.a, z1 > 0), fdr.cut=1e-06, z.cut = 0.03, nudge=0, down.sample=0.97, show.legend=TRUE, two.sided=FALSE, ylim=c(1, 1e-30)) + 
		xlim(0, 0.075) +
		guides(alpha=FALSE, colour="legend") + labs(colour="compound class") +
		xlab("mean cosine similarity") +
		scale_colour_npg(na.value="grey40") +
		ggtitle("RCs")
	,
	width = 8, height = 5,
	file = insert(out.fname, c("cmap-cosine", "camera", "volcano", "all"), ext="pdf")
);

####

plot_compound <- function(compound) {
	idx <- groups$trt_cp[[paste0("trt_cp,", compound)]];
	ss.sub <- ss[, idx];
	d.ss <- melt(ss.sub, varnames=c("clone", "sig"));

	qdraw(
		ggplot(d.ss, aes(x=clone, y=value)) + geom_boxplot() + theme_bw() +
			ylab("cosine similarity") +
			geom_hline(yintercept=0, colour="grey40", linetype=2) +
			ggtitle(compound)
		,
		width = 4, height = 5,
		file = insert(out.fname, c("cos-sim", compound), ext="pdf")
	);
}

# negative correlated
plot_compound("obatoclax");
plot_compound("prostratin");
plot_compound("vincristine");
plot_compound("vinblastine");

# budesonide: anti-inflammatory
# AS-605240: PI3K inhibitor
# PD-98059: MEK inhibitor
# selumetinib: MEK inhibitor
# trametinib: MEK inhibitor
# ERK-inhibitor-11E: ERK inhibitor
# PD-184352: MEK inhibitor
# PD-0325901: MEK inhibitor
# AZD-8330: MEK inhibitor
# AZ-628: RAF inhibitor
# U0126: MEK inhibitor
# AS-703026: MEK inhibitor

# positively correlated
plot_compound("selumetinib");
plot_compound("trametinib");
plot_compound("AS-605240");

# others

plot_compound("paclitaxel");
plot_compound("docetaxel");

groups$trt_cp[[paste0("trt_cp,", "paclitaxel")]];
