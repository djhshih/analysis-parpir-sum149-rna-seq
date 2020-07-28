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

out.fname <- in.fname;
out.fname$path <- NULL;
out.fname$ext <- NULL;

x <- qread(in.fname);

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


# Analysis ###################################################################

# cov(X, Y) = E[ (X - E[X]) (Y - E[Y]) ]
# cor(X, Y) = cov(X, Y) / ( sd(X) sd(Y) )
#           = E[ ((X - E[X]) / sd(X)) ((Y - E[Y]) / sd(Y)) ]
#           = cos(\tilde{X}, \tilde{Y})
#           where \tilde{X} = X - E[X]
#                 \tilde{Y} = Y - E[Y]
# cos(X, Y) = X \cdot Y / ( || X || || Y || )

# Both cov and cor subtract the mean, which is undesirable
# because in our case, the mean contains useful signal

# cosine similarity is the uncentered pearson correlation

fnorm <- function(x) {
	sqrt(sum( x * x, na.rm=TRUE ))
}

cos_sim <- function(x, y) {
	sum(x * y, na.rm=TRUE) / ( fnorm(x) * fnorm(y) )
}

columns_by_columns <- function(x, y, f, ...) {
	apply(x, 2,
		function(x1) {
			apply(y, 2,
				function(y1) {
					f(x1, y1, ...)
				}
			)
		}
	)
}

euclidean_distance <- function(x, y) {
	d <- x - y;
	sqrt(sum(d*d))
}

# similarity matrix between each column of x and each column of r
cos_sim_all <- function(r, x) {
	columns_by_columns(r, x, cos_sim)
}

cor_all <- function(r, x) {
	columns_by_columns(r, x, cor)
}

euclidean_distance_all <- function(r, x) {
	columns_by_columns(r, x, euclidean_distance)
}

#ss <- qread(tag(out.fname, "lincs-lm-sim", ext="rds"));
ss <- apply_gctx(gct, cos_sim_all, length(cid), x=x);
ss <- do.call(cbind, ss);

qwrite(ss, insert(out.fname, "lincs-lm-sim", ext="rds"));

min.idx <- apply(ss, 1, which.min);
ss[, min.idx]
pheno[min.idx, ]

max.idx <- apply(ss, 1, which.max);
ss[, max.idx]
pheno[max.idx, ]

# make null data matrix by bootrapping both samples and features
# x is m by n matrix with m features and n samples
bootstrap <- function(x, n=ncol(x)) {
	#na.idx <- apply(x, 1, function(x1) any(is.na(x1)));
	#x <- x[!na.idx, ];
	#cidx <- sample(1:ncol(x), n, replace=TRUE);
	ridx <- sample(1:nrow(x), nrow(x), replace=TRUE);
	#x[ridx, cidx]
	x[ridx, ]
}

set.seed(1234);
xb <- bootstrap(x);

# confirm that marginal distributions are preserved
hist(x, breaks=100)
hist(xb, breaks=100)
qqplot(x, xb, type="l");
abline(a=0, b=1, col="red");

# confirm that null sample is uncorrelated with original data
plot(c(x), c(xb), pch=".")
cor(c(x), c(xb), use="complete.obs")

ss.null <- apply_gctx(gct, cos_sim_all, length(cid), x=xb);
ss.null <- do.call(cbind, ss.null);

qwrite(ss.null, insert(out.fname, "lincs-lm-sim-null", ext="rds"));

