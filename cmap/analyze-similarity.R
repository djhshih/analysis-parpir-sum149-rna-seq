library(cmapR)
library(mmalign)
library(dplyr)
library(limma)
library(ComplexHeatmap)
library(ggplot2)
library(reshape2)
#library(viridis)
library(RColorBrewer)
library(io)


lincs.indir <- "~/data/cmap/phase1";

source(file.path(lincs.indir, "common.R"));

# use only landmark genes
gct <- file.path(lincs.indir, "lincs-sig_lm_n473647x978.gctx");
pheno <- qread(file.path(lincs.indir, "lincs-sig-info.rds"));
cell <- qread(file.path(lincs.indir, "lincs-cell-info.rds"));
fannot <- qread(file.path(lincs.indir, "lincs-gene-info_lm.rds"));
groups <- qread(file.path(lincs.indir, "lincs-sig-groups.rds"));

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


##############################################################################

clones <- sub(".*clone(C.*)", "\\1", rownames(ss));
rownames(ss) <- clones;

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

camera_transform <- function(ss, group) {
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

ms.cp <- mean_transform(ss, groups$trt_cp);
mx.cp <- max_transform(ss, groups$trt_cp);
mn.cp <- min_transform(ss, groups$trt_cp);
mp.cp <- median_percentile_transform(ss, groups$trt_cp);

es.cp <- camera_transform(ss, groups$trt_cp);

#idx <- apply(es.cp, 1, function(lq) any(lq > 6));
idx <- apply(es.cp, 1, function(lq) any(lq > 10));
es.cp.sub <- es.cp[idx, , drop=FALSE];
mp.cp.sub <- mp.cp[idx, , drop=FALSE];
ms.cp.sub <- ms.cp[idx, , drop=FALSE];
mx.cp.sub <- mx.cp[idx, , drop=FALSE];

remap_samples <- function(samples) {
	sample.id <- as.integer(sub("C", "", samples));
	si <- order(sample.id);
	factor(sample.id, levels=sample.id[si], label=paste0("RC", 1:length(sample.id)))
}

simplify_features <- function(features) {
	sub("^.+,", "", features)
}

simplify_names <- function(z) {
	colnames(z) <- remap_samples(colnames(z));
	rownames(z) <- simplify_features(rownames(z));
	z
}

pal <- rev(brewer.pal(7, "RdBu"));

#colf <- circlize::colorRamp2(c(0, 100), c(pal[4], pal[7]));
colf <- circlize::colorRamp2(c(0, 50, 100), c(pal[1], pal[4], pal[7]));

#pdf(file=tag(out.fname, c("cos-sim-percentile", "up", "heatmap"), ext="pdf"), width=5, height=12);
pdf(file=tag(out.fname, c("cos-sim-percentile", "up", "heatmap"), ext="pdf"), width=4, height=6);
Heatmap(simplify_names(mp.cp.sub), col=colf, heatmap_legend_param=list(title="percentile"))
dev.off();

colf <- circlize::colorRamp2(c(-0.15, 0, 0.15), c(pal[1], pal[4], pal[7]));
Heatmap(mx.cp.sub, col=colf)


qdraw(
	{
		plot(density(ss.null), las=1, xlab="Cosine similarity", main=NA, col="grey30", lwd=2, bty="n")
		lines(density(ss), las=1, xlab="Cosine similarity", main=NA, col="firebrick", lwd=2)
	},
	width = 4, height = 3,
	file = tag(out.fname, c("cos-sim", "density"), ext="pdf")
);

####

idx <- apply(es.cp, 1, function(lq) any(lq < -6));
#idx <- apply(es.cp, 1, function(lq) any(lq < -15));
es.cp.sub <- es.cp[idx, , drop=FALSE];
mp.cp.sub <- mp.cp[idx, , drop=FALSE];
ms.cp.sub <- mp.cp[idx, , drop=FALSE];
mn.cp.sub <- mn.cp[idx, , drop=FALSE];

pal <- rev(brewer.pal(7, "RdBu"));

#colf <- circlize::colorRamp2(c(0, 100), c(pal[4], pal[7]));
colf <- circlize::colorRamp2(c(0, 50, 100), c(pal[1], pal[4], pal[7]));
pdf(file=tag(out.fname, c("cos-sim-percentile", "down", "heatmap"), ext="pdf"), width=4, height=24);
Heatmap(simplify_names(mp.cp.sub), col=colf, heatmap_legend_param=list(title="percentile"))
dev.off();

colf <- circlize::colorRamp2(c(-0.15, 0, 0.15), c(pal[1], pal[4], pal[7]));
Heatmap(mn.cp.sub, col=colf)


####

idx <- apply(es.cp, 1, min) > 1;
es.cp.min.sub <- es.cp[idx, , drop=FALSE];
mp.cp.min.sub <- mp.cp[idx, , drop=FALSE];



varnames <- c("treatment", "clone");
d.es <- melt(es.cp.min.sub, varnames=varnames);
d.mp <- melt(mp.cp.min.sub, varnames=varnames);

# up FDR needs to be negated
d <- left_join(mutate(d.es, q=10^(-value)) %>% select(-value), rename(d.mp, percentile=value));

# rank treatments by percentile
treatments <- levels(d$treatment);
means <- d %>% group_by(treatment) %>% summarize(y=mean(percentile));
idx <- order( -means$y );
d$treatment <- factor(d$treatment, levels=treatments[idx]);

qdraw(
	ggplot(d, aes(x=clone, y=percentile, fill=q)) +
		geom_bar(stat="identity") + theme_bw() +
		facet_wrap(~ treatment, ncol=1, strip.position="top") +
		#facet_wrap(~ treatment, nrow=1, strip.position="top") +
		scale_fill_gradient(low="#56B1F7", high="#132B43", trans="log10") +
		ylab("median percentile of cosine similarity")
	,
	width = 3.5, height = 10,
	file = insert(out.fname, c("trt_cp", "up", "min-q"), ext="pdf")
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
