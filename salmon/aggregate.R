library(io)
library(GSVA)
library(mmalign)
library(dplyr)
library(reshape2)

# cell line: SUM149

out.fname <- filename("parpi-resist");

x <- qread("parpi-resist_tpm.rds");
fannot <- qread("parpi-resist_fannot.rds");

rownames(fannot) <- fannot$ensembl_transcript;

hist(x, breaks=1000)

aggregate_genes <- function(x, fannot, stat=mean) {
	# aggregate genes
	require(data.table);
	y <- data.table(x);
	y <- y[, lapply(.SD, stat), by = fannot$gene_name];

	# convert back into matrix
	genes <- y[[1]];
	y <- as.matrix(y[,-1]);
	rownames(y) <- genes;

	y
}


y.mean <- aggregate_genes(log(x + 1), fannot);

hist(y.mean, breaks=1000)
hist(y.mean[y.mean > 0], breaks=1000)

qwrite(y.mean, insert(out.fname, "ltpm-genes-mean", ext="rds"));
qwrite(y.mean, insert(out.fname, "ltpm-genes-mean", ext="mtx"));

y.max <- aggregate_genes(log(x + 1), fannot, max);

hist(y.max, breaks=1000)
hist(y.max[y.max > 0], breaks=1000)

qwrite(y.max, insert(out.fname, "ltpm-genes-max", ext="rds"));
qwrite(y.max, insert(out.fname, "ltpm-genes-max", ext="mtx"));

