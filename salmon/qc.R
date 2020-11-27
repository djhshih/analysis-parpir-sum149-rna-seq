library(io)
library(mmalign)
library(sva)

out.fname <- filename("parpi-resist");

tpm <- qread("parpi-resist_tpm.rds");
fannot <- qread("parpi-resist_fannot.rds");
pheno <- qread("../annot/sample-info_parpi-resist_stage2.tsv");

barcodes <- strsplit(as.character(pheno$barcode), "+", fixed=TRUE);
pheno$barcode1 <- unlist(lapply(barcodes, function(x) x[1]));
pheno$barcode2 <- unlist(lapply(barcodes, function(x) x[2]));

stopifnot(colnames(tpm) == as.character(pheno$sample_id))

ltpm <- log(tpm + 1);

p <- pca(ltpm);

qdraw(
	mmalign::pca_plot(p, pheno, aes(colour=batch, shape=factor(lane)))
	,
	width = 6,
	file = insert(out.fname, c("pca", "1-2", "batch", "lane"), ext="pdf")
);

mmalign::pca_plot(p, pheno, aes(colour=batch, shape=factor(lane)), dims=3:4)

shape.scale <- scale_shape_manual(values=c(20, 21));
mmalign::pca_plot(p, pheno, aes(shape=batch, colour=factor(lane))) + shape.scale
mmalign::pca_plot(p, pheno, aes(shape=batch, colour=factor(lane)), dims=3:4) + shape.scale

####

tpm.prot <- tpm[fannot$gene_type == "protein_coding", ];
ltpm.prot <- log(tpm.prot + 1);
rsums <- rowSums(ltpm.prot);

hist(ltpm.prot, breaks=100)
hist(rsums, breaks=100)

hist(ltpm.prot[ltpm.prot > 0], breaks=500)
hist(ltpm.prot[ltpm.prot > 0], breaks=500, xlim=c(0, 0.5))

summary(ltpm.prot)

prop.expressed.cut <- 0.10;
prop.expressed <- rowMeans(ltpm.prot > 0.2);
hist(prop.expressed, breaks=50)

ltpm.prot.c <- ltpm.prot[prop.expressed > prop.expressed.cut, ];

dim(ltpm.prot.c)

p <- pca(ltpm.prot);

mmalign::pca_plot(p, pheno, aes(colour=batch, shape=factor(lane)))

mmalign::pca_plot(p, pheno, aes(colour=batch, shape=factor(lane)), dims=3:4)

shape.scale <- scale_shape_manual(values=c(20, 21));
mmalign::pca_plot(p, pheno, aes(shape=batch, colour=factor(lane))) + shape.scale
mmalign::pca_plot(p, pheno, aes(shape=batch, colour=factor(lane)), dims=3:4) + shape.scale


shape.scale <- scale_shape_manual(values=c(20, 24, 4));

qdraw(
	mmalign::pca_plot(p, pheno, aes(shape=treatment, colour=resistance)) +
		shape.scale
	,
	width = 6,
	file = insert(out.fname, c("pca", "1-2"), ext="pdf")
);

qdraw(
	mmalign::pca_plot(p, pheno, aes(shape=treatment, colour=resistance), dims=3:4) +
		shape.scale
	,
	width = 6,
	file = insert(out.fname, c("pca", "3-4"), ext="pdf")
);

qdraw(
	mmalign::pca_plot(p, pheno, aes(shape=treatment, colour=resistance), dims=5:6) +
		shape.scale
	,
	width = 6,
	file = insert(out.fname, c("pca", "5-6"), ext="pdf")
);

shape.scale <- scale_shape_manual(values=c(1:8));
mmalign::pca_plot(p, pheno, aes(shape=clone, colour=resistance)) + shape.scale
mmalign::pca_plot(p, pheno, aes(shape=clone, colour=resistance), dims=3:4) + shape.scale

shape.scale <- scale_shape_manual(values=c(1:8));
mmalign::pca_plot(p, pheno, aes(shape=clone, colour=batch)) + shape.scale
mmalign::pca_plot(p, pheno, aes(shape=clone, colour=batch), dims=3:4) + shape.scale

shape.scale <- scale_shape_manual(values=c(1:8));

qdraw(
	mmalign::pca_plot(p, pheno, aes(shape=clone, colour=treatment)) + shape.scale
	,
	width = 6,
	file = insert(out.fname, c("pca", "1-2", "clone", "treatment"), ext="pdf")
);

qdraw(
	mmalign::pca_plot(p, pheno, aes(shape=clone, colour=treatment), dims=3:4) + shape.scale
	,
	width = 6,
	file = insert(out.fname, c("pca", "3-4", "clone", "treatment"), ext="pdf")
);

qdraw(
	mmalign::pca_plot(p, pheno, aes(shape=clone, colour=treatment), dims=4:5) + shape.scale
	,
	width = 6,
	file = insert(out.fname, c("pca", "4-5", "clone", "treatment"), ext="pdf")
);

shape.scale <- scale_shape_manual(values=c(1:8));
mmalign::pca_plot(p, pheno, aes(shape=clone, colour=barcode1)) + shape.scale
mmalign::pca_plot(p, pheno, aes(shape=clone, colour=barcode1), dims=3:4) + shape.scale

shape.scale <- scale_shape_manual(values=c(1:8));
mmalign::pca_plot(p, pheno, aes(shape=clone, colour=barcode2)) + shape.scale
mmalign::pca_plot(p, pheno, aes(shape=clone, colour=barcode2), dims=3:4) + shape.scale

mmalign::pca_plot(p, pheno, aes(colour=log(fold_resistance)))
mmalign::pca_plot(p, pheno, aes(colour=log(fold_resistance)), dims=3:4)
mmalign::pca_plot(p, pheno, aes(colour=log(fold_resistance)), dims=5:6)
mmalign::pca_plot(p, pheno, aes(colour=log(fold_resistance)), dims=7:8)
mmalign::pca_plot(p, pheno, aes(colour=log(fold_resistance)), dims=9:10)

mmalign::pca_plot(p, pheno, aes(colour=log(n_pf_clusters)))
mmalign::pca_plot(p, pheno, aes(colour=log(n_pf_clusters)), dims=3:4)
mmalign::pca_plot(p, pheno, aes(colour=log(n_pf_clusters)), dims=5:6)


d <- as.dist(1 - cor(ltpm.prot));
hc <- hclust(d, method="average");
plot(hc)

qdraw(
	{
		plot(as.dendrogram(hc, hang=0.01), horiz=TRUE, xlab="1 - correlation")
	},
	height = 10,
	file = insert(out.fname, "hclust", ext="pdf")
);

#----

hist(ltpm.prot, breaks=1000)
hist(ltpm.prot.c, breaks=1000)

boxplot(ltpm, outline=FALSE, las=2)
boxplot(ltpm.prot, outline=FALSE, las=2)
boxplot(ltpm.prot.c, outline=FALSE, las=2)

#---

library(vsn)
library(DESeq2)

qdraw(
	meanSdPlot(ltpm),
	file = insert(out.fname, c("ltpm", "meansd"), ext="pdf")
);
#meanSdPlot(ltpm.prot)
#meanSdPlot(ltpm.prot.c)

#vtpm.prot <- vsn2(tpm.prot);
# miserable!
#qdraw(
#	meanSdPlot(vtpm.prot),
#	file = insert(out.fname, c("vtpm", "meansd"), ext="pdf")
#);

mod <- model.matrix(~ treatment + clone, data=pheno)
#mod <- model.matrix(~ clone, data=pheno)
#mod <- model.matrix(~ 1, data=pheno)
#ltpm.prot.bc <- ComBat(ltpm.prot.c, pheno$batch, mod, par.prior=FALSE);
ltpm.prot.bc <- ComBat(ltpm.prot.c, pheno$batch, mod);

hist(ltpm.prot.bc, breaks=100)

p <- pca(ltpm.prot.bc);

qdraw(
	mmalign::pca_plot(p, pheno, aes(colour=batch, shape=factor(lane)))
	,
	width = 6,
	file = insert(out.fname, c("cb", "pca", "1-2", "batch", "lane"), ext="pdf")
);

