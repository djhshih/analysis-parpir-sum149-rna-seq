library(io)
library(cmapR)

lincs.indir <- "~/data/cmap/phase1";

gct <- file.path(lincs.indir, "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx");
pheno <- qread(file.path(lincs.indir, "lincs-sig-info.rds"));
fannot <- qread(file.path(lincs.indir, "lincs-gene-info.rds"));

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


chunk.size <- 10000;

n <- length(cid);

bstarts <- seq(1, n, by=chunk.size);

ss <- lapply(bstarts,
	function(bstart) {

		message("block start ", bstart)

		istart <- bstart;
		iend <- min(bstart + chunk.size - 1, n);

		ds <- parse.gctx(gct, cid=istart:iend, matrix_only=TRUE)@mat;

		# similarity matrix between x and d
		s <- apply(ds, 2,
			function(r1) {
				apply(x, 2,
					function(x1) {
						cos_sim(x1, r1)
					}
				)
			}
		);
		
		s
	}
);

ss2 <- do.call(cbind, ss);

qwrite(ss2, insert(out.fname, "lincs-sim", ext="rds"));

min.idx <- apply(ss2, 1, which.min);
ss2[, min.idx]
pheno[min.idx, ]

max.idx <- apply(ss2, 1, which.max);
ss2[, max.idx]
pheno[max.idx, ]

