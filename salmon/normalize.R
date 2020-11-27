library(io)
library(DESeq2)
library(mmalign)
library(sva)
library(vsn)
library(dplyr)
library(reshape2)

source("../R/preamble.R")

# cell line: SUM149

out.fname <- filename("parpi-resist");

x <- qread("parpi-resist_counts.rds");
fannot <- qread("parpi-resist_fannot.rds");
pheno <- setup_pheno(qread("../annot/sample-info_parpi-resist_stage2.tsv"));

stopifnot(colnames(x) == as.character(pheno$sample_id))

counts <- round(x);
mode(counts) <- "integer";

# compare all levels against DMSO
pheno$treatment <- relevel(pheno$treatment, "DMSO");

# TODO change s.t. the most extreme transcript is chosen!
deseqresults_to_rnk <- function(res, fannot) {
	rnk <- data.frame(gene_name = fannot$gene_name, stat = res$stat);
	rnk <- rnk[!duplicated(rnk$gene_name) & !is.na(rnk$stat), ];
	rnk[order(rnk$stat, decreasing=TRUE), ]
}

#----

dds.treatment <- DESeqDataSetFromMatrix(
	countData = counts, colData = pheno,
	design = ~ flowcell + treatment
);
dds.treatment <- DESeq(dds.treatment);
qwrite(dds.treatment, insert(out.fname, c("deseq", "treatment"), ext="rds"));

resultsNames(dds.treatment)

res <- results(dds.treatment, name="treatment_Talazoparib_vs_DMSO");
qwrite(res, insert(out.fname, c("deseq-res", "treatment"), ext="rds"));

rnk <- deseqresults_to_rnk(res, fannot)
qwrite(rnk, insert(out.fname, c("treatment"), ext="rnk"));

res.s <- res[order(res$padj), ];
res.sf <- res.s[which(res.s$padj < 0.05), ];
head(res.sf)

txs.up <- rownames(res.sf)[res.sf$log2FoldChange > 1];
as.data.frame(cbind(fannot[txs.up, "gene_name"], res.sf[txs.up, ]))

txs.down <- rownames(res.sf)[res.sf$log2FoldChange < -1];
as.data.frame(cbind(fannot[txs.down, "gene_name"], res.sf[txs.down, ]))

####

dds.resistance <- DESeqDataSetFromMatrix(
	countData = counts, colData = pheno,
	design = ~ flowcell + resistance
);
dds.resistance <- DESeq(dds.resistance);
qwrite(dds.resistance, insert(out.fname, c("deseq", "resistance"), ext="rds"));

resultsNames(dds.resistance)

res <- results(dds.resistance, name="resistance_Resistant_vs_Sensitive");
qwrite(res, insert(out.fname, c("deseq-res", "resistance"), ext="rds"));

rnk <- deseqresults_to_rnk(res, fannot)
qwrite(rnk, insert(out.fname, c("resistance"), ext="rnk"));

res.s <- res[order(res$padj), ];
res.sf <- res.s[which(res.s$padj < 0.05), ];
head(res.sf)

txs.up <- rownames(res.sf)[res.sf$log2FoldChange > 1];
as.data.frame(cbind(fannot[txs.up, "gene_name"], res.sf[txs.up, ]))

txs.down <- rownames(res.sf)[res.sf$log2FoldChange < -1];
as.data.frame(cbind(fannot[txs.down, "gene_name"], res.sf[txs.down, ]))

####

dds.treatment.b.clone <- DESeqDataSetFromMatrix(
	countData = counts, colData = pheno,
	design = ~ treatment + clone
);
dds.treatment.b.clone <- DESeq(dds.treatment.b.clone);
qwrite(dds.treatment.b.clone, insert(out.fname, c("deseq", "treatment-clone"), ext="rds"));

resultsNames(dds.treatment.b.clone)

# these contrasts based on the DMSO treatment
contrasts.clones <- grep("clone_", resultsNames(dds.treatment.b.clone), value=TRUE);
names(contrasts.clones) <- contrasts.clones;

rnks <- lapply(contrasts.clones,
	function(contrast) {
		res <- results(dds.treatment.b.clone, name=contrast);
		tag <- gsub("_", "-", tolower(contrast));
		qwrite(res, insert(out.fname, c("deseq-res", tag), ext="rds"));
		rnk <- deseqresults_to_rnk(res, fannot=fannot);
		qwrite(rnk, insert(out.fname, tag, ext="rnk"));

		rnk
	}
);

genes <- levels(rnks[[1]]$gene_name);

rnk.mtx <- do.call(cbind, lapply(rnks, function(x) x[[2]][match(genes, x[[1]])] ));
rownames(rnk.mtx) <- genes;

hist(rnk.mtx, breaks=1000)

qwrite(rnk.mtx, insert(out.fname, c("deseq-stat", "clones-vs-parental"), ext="mtx"));

#

# overall treatment effect across parental and all resistant clones
res <- results(dds.treatment.b.clone, name="treatment_Talazoparib_vs_DMSO");
qwrite(res, insert(out.fname, c("deseq-res", "treatment-clone", "treatment"), ext="rds"));

rnk <- deseqresults_to_rnk(res, fannot)
qwrite(rnk, insert(out.fname, c("treatment-clone", "treatment"), ext="rnk"));

#


#----

# remove treatment = None in order to make full rank design matrix
idx <- pheno$treatment != "None"
pheno2 <- pheno[idx, ];
pheno2$treatment <- droplevels(pheno2$treatment);
counts2 <- counts[, idx];

dds.treatment.c.clone <- DESeqDataSetFromMatrix(
	countData = counts2, colData = pheno2,
	design = ~ treatment * clone
);
dds.treatment.c.clone <- DESeq(dds.treatment.c.clone);
qwrite(dds.treatment.c.clone, insert(out.fname, c("deseq", "treatment-clone-interaction"), ext="rds"));

resultsNames(dds.treatment.c.clone)

# these contrasts based on the DMSO treatment
contrasts.clones <- grep("clone_", resultsNames(dds.treatment.c.clone), value=TRUE);
names(contrasts.clones) <- contrasts.clones;

rnks <- lapply(contrasts.clones,
	function(contrast) {
		res <- results(dds.treatment.c.clone, name=contrast);
		tag <- gsub("_", "-", tolower(contrast));
		qwrite(res, insert(out.fname, c("deseq-res", "treatment-clone-interaction", tag), ext="rds"));
		rnk <- deseqresults_to_rnk(res, fannot=fannot);
		qwrite(rnk, insert(out.fname, c("treatment-clone-interaction", tag), ext="rnk"));

		rnk
	}
);

# collect all the main clone effects together,
# ensuring that the genes are in the correct order

rnk2.mtx <- do.call(cbind, lapply(rnks, function(x) x[[2]][match(genes, x[[1]])] ));
rownames(rnk2.mtx) <- genes;

hist(rnk2.mtx, breaks=1000)

cor(as.numeric(rnk.mtx), as.numeric(rnk2.mtx), use="complete.obs")
plot(as.numeric(rnk.mtx), as.numeric(rnk2.mtx), pch=".")

qwrite(rnk2.mtx, insert(out.fname, c("deseq-stat", "treatment-clone-interaction", "clones-vs-parental"), ext="mtx"));

#

# treatment effect in parental
res <- results(dds.treatment.c.clone, name="treatment_Talazoparib_vs_DMSO");
qwrite(res, insert(out.fname, c("deseq-res", "treatment-clone-interaction", "treatment"), ext="rds"));

rnk <- deseqresults_to_rnk(res, fannot)
qwrite(rnk, insert(out.fname, c("treatment-clone-interaction", "treatment"), ext="rnk"));

#

contrasts.inters <- grep("treatmentTalazoparib.clone", resultsNames(dds.treatment.c.clone), value=TRUE);
names(contrasts.inters) <- contrasts.inters;

rnks <- lapply(contrasts.inters,
	function(contrast) {
		res <- results(dds.treatment.c.clone, name=contrast);
		tag <- gsub(".", "-", tolower(contrast), fixed=TRUE);
		qwrite(res, insert(out.fname, c("deseq-res", "treatment-clone-interaction", tag), ext="rds"));
		rnk <- deseqresults_to_rnk(res, fannot=fannot);
		qwrite(rnk, insert(out.fname, c("treatment-clone-interaction", tag), ext="rnk"));

		rnk
	}
);

rnk3.mtx <- do.call(cbind, lapply(rnks, function(x) x[[2]][match(genes, x[[1]])] ));
rownames(rnk3.mtx) <- genes;

hist(rnk3.mtx, breaks=1000)

qwrite(rnk3.mtx, insert(out.fname, c("deseq-stat", "treatment-clone-interaction", "interactions"), ext="mtx"));

