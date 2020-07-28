library(io)
library(cmapR)
library(mmalign)
library(dplyr)


lincs.indir <- "~/data/cmap/phase1";

source(file.path(lincs.indir, "common.R"));

# use only landmark genes
gct <- file.path(lincs.indir, "lincs-sig_lm_n473647x978.gctx");
pheno <- qread(file.path(lincs.indir, "lincs-sig-info.rds"));
cell <- qread(file.path(lincs.indir, "lincs-cell-info.rds"));
fannot <- qread(file.path(lincs.indir, "lincs-gene-info_lm.rds"));

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


control.idx <- pheno$pert_type %in% control.types;
ss.ctl <- ss[, control.idx];

# distributions of similarity scores are remarkably Gaussian
qdraw(
	{
		par(mfrow=c(3,1))
		xlim <- c(-0.2, 0.2);
		hist(ss, breaks=100, freq=FALSE, xlim=xlim)
		curve(dnorm(x, mean(ss), sd(ss)), col="red", add=TRUE)
		hist(ss.ctl, breaks=100, freq=FALSE, xlim=xlim)
		curve(dnorm(x, mean(ss.ctl), sd(ss.ctl)), col="red", add=TRUE)
		hist(ss.null, breaks=100, freq=FALSE, xlim=xlim)
		curve(dnorm(x, mean(ss.null), sd(ss.null)), col="red", add=TRUE)
	},
	width = 5, height = 8,
	file = insert(out.fname, c("lincs-lm-sim", "hist"), ext="pdf")
);


qdraw(
	{
		par(mfrow=c(3,1))
		qqnorm(sort(sample(c(ss), 1000)), type="l");
		qqnorm(sort(sample(c(ss.ctl), 1000)), type="l");
		qqnorm(sort(sample(c(ss.null), 1000)), type="l");
	},
	width = 5, height = 8,
	file = insert(out.fname, c("lincs-lm-sim", "qqnorm"), ext="pdf")
);


null.mean <- mean(ss.null);
null.sd <- sd(ss.null);

primary.site <- "breast";
cells <- as.character(cell$cell_id[which(cell$primary_site == primary.site)]);

pert.idx <-
	! pheno$pert_type %in% control.types & 
	pheno$cell_id %in% cells &
	pheno$pert_type %in% c("trt_cp");

table(pert.idx)

pheno.f <- pheno[pert.idx, ];
ss.f <- ss[, pert.idx];

ss.p.neq2 <- pmin(2 * pnorm(-abs(ss.f), mean=mean(ss.ctl), sd=sd(ss.ctl)), 1);
ss.p.gt2 <- pnorm(ss.f, mean=mean(ss.ctl), sd=sd(ss.ctl), lower.tail=FALSE);

hist(ss.p.neq2, breaks=100, xlim=c(0, 1))
hist(ss.p.gt2, breaks=100, xlim=c(0, 1))

ss.p.neq <- pmin(2 * pnorm(-abs(ss.f), mean=null.mean, sd=null.sd), 1);

qdraw(
	{
		hist(ss.p.neq[ss.p.neq < 1], breaks=100, xlim=c(0, 1));
	},
	file = insert(out.fname, c("lincs-lm-sim-p", "hist"), ext="pdf")
)

# NB We don't expect most correlations to be zero

source("qqunif.plot.R");

#qqunif.plot(as.numeric(ss.p.neq[ss.p.neq < 1]), aspect="fill")

qdraw(
	{
		ggd.qqplot(as.numeric(ss.p.neq[ss.p.neq < 1]))
	},
	file = insert(out.fname, c("lincs-lm-sim-p", "qq"), ext="png")
)

#qqunif.plot(as.numeric(ss.p.neq2));
#ggd.qqplot(as.numeric(ss.p.neq2))

#inflation_factor <- function(p) {
#	chisq <- qchisq(1 - p, 1);
#  median(chisq) / qchisq(0.5, 1)
#}

#inflation_factor(ss.p.neq[ss.p.neq < 1])
#inflation_factor(ss.p.neq2)


fdr <- 0.1;

hits <- lapply(
	1:nrow(ss.p.neq),
	function(i) {
		# vehicle samples (DMSO) occur in samples,
		# presumably because their expression change directions are random...
		# FIXME how to verify this?
		q <- p.adjust(ss.p.neq[i, ], method="BH");
		idx <- q < fdr;
		d <- cbind(
			similarity = ss.f[i, idx],
			p = ss.p.neq[i, idx],
			q = q[idx],
			#pheno.f[idx, c("pert_iname", "pert_type", "cell_id", "pert_idose", "pert_time", "sig_id")]
			pheno.f[idx, c("pert_iname", "cell_id", "sig_id")]
		);
		rownames(d) <- NULL;
		d[order(d$similarity), ]
	}
);


data.frame(sample=colnames(x), n=unlist(lapply(hits, nrow)))

n.tops <- 20;
tops <- lapply(
	1:nrow(ss.p.neq),
	function(i) {
		q <- p.adjust(ss.p.neq[i, ], method="BH");
		idx <- order(ss.p.neq[i, ])[1:n.tops];
		d <- cbind(
			similarity = ss.f[i, idx],
			p = ss.p.neq[i, idx],
			q = q[idx],
			#pheno.f[idx, c("pert_iname", "pert_type", "cell_id", "pert_idose", "pert_time", "sig_id")]
			pheno.f[idx, c("pert_iname", "cell_id", "sig_id")]
		);
		rownames(d) <- NULL;
		d[order(d$similarity), ]
	}
);



mek <- unlist(lapply(hits, function(h) {
	if (nrow(h) > 0) {
		h2 <- h[h$similarity > 0, ];
		"selumetinib" %in% as.character(h2$pert_iname)
	} else {
		FALSE
	}
}));


na.idx <- apply(x, 1, function(x1) any(is.na(x1)));
x.f <- x[!na.idx, ];

x.f.a <- cbind(x.f, Parental=0)

x.pca <- pca(x.f.a);

x.pheno <- data.frame(
	clone = sub(".*clone(.*)", "\\1", colnames(x.f.a)),
	mek = c(mek, FALSE)
);

pca_text_plot(x.pca, pheno=x.pheno, aes(label=clone, colour=mek))
pca_text_plot(x.pca, pheno=x.pheno, aes(label=clone, colour=mek), dims=2:3)
pca_text_plot(x.pca, pheno=x.pheno, aes(label=clone, colour=mek), dims=3:4)
pca_text_plot(x.pca, pheno=x.pheno, aes(label=clone, colour=mek), dims=4:5)
pca_text_plot(x.pca, pheno=x.pheno, aes(label=clone, colour=mek), dims=5:6)
pca_text_plot(x.pca, pheno=x.pheno, aes(label=clone, colour=mek), dims=6:7)


# look for dose-dependent effect or at least consistent effect

j <- 2;

hits.sig1 <- hits[[j]]$sig_id;

pheno.hits <- filter(pheno.f, sig_id %in% hits.sig1) %>% select(-distil_id) %>% left_join(hits[[j]])
filter(pheno.hits, cell_id == "MCF10A", pert_iname == "selumetinib")

pert.ids <- unique(pheno.f$pert_id[pheno.f$sig_id %in% hits.sig1]);

sigs.sel <- pheno.f$sig_id[pheno.f$pert_id %in% pert.ids];

y <- parse.gctx(gct, cid=sigs.sel, matrix_only=TRUE)@mat;

y.pheno <- pheno[match(sigs.sel, pheno$sig_id), ]

y.f <- y[!na.idx, ];


y.pca <- pca(y.f);
pca_text_plot(y.pca, pheno=y.pheno, aes(label=pert_iname))

x.t <- pca_transform(x.f, y.pca);
pca_text_plot(x.t, pheno=x.pheno, aes(label=clone))

x.pheno[j, ]

z <- cbind(x.f, y.f);
z.pca <- pca(z);
z.pheno <- rbind(
	data.frame(
		label = x.pheno$clone,
		group = "case"
	),
	data.frame(
		label = y.pheno$pert_iname,
		group = "cmap"
	)
);
	
pca_text_plot(z.pca, pheno=z.pheno, aes(label=label, colour=group))
pca_text_plot(z.pca, pheno=z.pheno, aes(label=label, colour=group), dims=3:4)
pca_text_plot(z.pca, pheno=z.pheno, aes(label=label, colour=group), dims=5:6)


s.xy <- cos_sim_all(x.f[, j, drop=FALSE], y.f);
hist(s.xy)

c.xy <- cor_all(x.f[, j, drop=FALSE], y.f);
hist(c.xy)

d.xy <- euclidean_distance_all(x.f[, j, drop=FALSE], y.f);
hist(d.xy)

stopifnot(y.pheno$sig_id == rownames(d.xy))

y.pheno2 <- data.frame(y.pheno, distance=as.numeric(d.xy), similarity=as.numeric(s.xy), correlation=as.numeric(c.xy));

y.pheno2.f <- filter(y.pheno2, pert_iname == "selumetinib");

ggplot(y.pheno2.f, aes(x=pert_dose, y=distance, colour=cell_id)) +
	geom_point() + theme_bw() + scale_x_log10()

ggplot(filter(y.pheno2.f, cell_id == "MCF10A"), aes(x=pert_dose, y=distance, colour=factor(pert_time))) +
	geom_point() + theme_bw() + scale_x_log10()

ggplot(filter(y.pheno2.f, cell_id == "MCF10A"), aes(x=pert_dose, y=similarity, colour=factor(pert_time))) +
	geom_point() + theme_bw() + scale_x_log10()

dd <- filter(y.pheno2.f, cell_id == "MCF10A");
select(dd, -distil_id)

with(dd, plot(distance, similarity))
with(dd, plot(correlation, similarity))

ggplot(filter(y.pheno2.f, cell_id == "MCF7"), aes(x=pert_dose, y=distance, colour=factor(pert_time))) +
	geom_point() + theme_bw() + scale_x_log10()

ggplot(filter(y.pheno2.f, cell_id == "MCF7"), aes(x=pert_dose, y=similarity, colour=factor(pert_time))) +
	geom_point() + theme_bw() + scale_x_log10()

dd <- filter(y.pheno2.f, cell_id == "MCF7");

summary( lm(distance ~ pert_dose, data=dd) )

ggplot(filter(y.pheno2.f, cell_id == "MDAMB231"), aes(x=pert_dose, y=distance, colour=factor(pert_time))) +
	geom_point() + theme_bw() + scale_x_log10()

ggplot(filter(y.pheno2.f, cell_id == "MDAMB231"), aes(x=pert_dose, y=similarity, colour=factor(pert_time))) +
	geom_point() + theme_bw() + scale_x_log10()


# NB perturbation vector can be in same direction as target vector, but
# the eucliean distance can increase with magnitude, when the perturbation
# vector is longer than the target vector

###



