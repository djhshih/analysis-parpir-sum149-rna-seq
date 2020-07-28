library(io)

in.fname <- "~/exigens/parpir-gdsc/expr/parpir-gdsc_limma_parpi-sen_rand-eff-cancer-type.tsv";
x <- qread(in.fname);

# resistant clone vs. parental (sensitive)
y <- qread("../salmon/parpi-resist_deseq-stat_treatment-clone-interaction_clones-vs-parental.mtx");

# match x to y
idx <- match(rownames(y), x$symbol);
xm <- x[idx, ];

stopifnot(all(xm$symbol == rownames(y), na.rm=TRUE))

# negative because x is sensitive vs. resistant
nxt <- - xm$t;

cors <- apply(y, 2,
	function(y1) {
		cor(y1, nxt, use="complete.obs")
	}
);
cors


y2 <- qread("../salmon/parpi-resist_resistance.rnk");

idx <- match(rownames(y), names(y2));
y2m <- y2[idx];

y2.cor <- cor(y2m, nxt, use="complete.obs");

stopifnot(all(names(y2m) == rownames(y), na.rm=TRUE))

my_plot_cors <- function(plotf) {
	par(mfrow=c(4,2), mai=c(0.5, 0.5, 0.2, 0.2), mgp=c(2.2, 1, 0))

	ylim <- c(-5, 5);
	xlim <- c(-12, 12);
	#ylab <- "GDSC resistant vs. sensitive t";
	#xlab <- "SUM149 resistant vs. parental t";
	ylab <- "GDSC t";
	xlab <- "%s t";

	clones <- sub("clone_(.*)_vs_Parental", "\\1", colnames(y));

	for (i in 1:ncol(y)) {
		plotf(y[, i], nxt, pch='.', ylim=ylim, xlim=xlim,
			ylab=ylab,
			xlab=sprintf(xlab, paste0("SUM149-", clones[i])),
			las=1,
		)
		legend("topleft", legend=sprintf("r = %.3f", cors[i]), bty="n", inset=0.01)
	}
	plotf(y2m, nxt, pch='.', ylab=ylab, xlab=sprintf(xlab, "SUM149 overall"), ylim=ylim, xlim=xlim, las=1)
	legend("topleft", legend=sprintf("r = %.3f", y2.cor), bty="n", inset=0.01)
}

qdraw(
	{
		my_plot_cors(smoothScatter);
	},
	width = 5, height = 9,
	file = "gdsc-vs-sum149_resistant-vs-sensitive_smooth-scatter.pdf"
);

