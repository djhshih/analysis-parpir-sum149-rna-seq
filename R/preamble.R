library(ggsci)

shift <- function(x, n=1) {
	if (n == 0) {
		x
	} else {
		c(tail(x, -n), head(x, n))
	}
}

clones <- c("Parental", "C2", "C17", "C19", "C20", "C26", "C29", "C30");
clone.cols <- shift(pal_d3()(length(clones)), -1);
names(clone.cols) <- clones;

setup_pheno <- function(pheno) {
	# DMSO is the reference
	pheno$treatment = factor(pheno$treatment, levels=c("None", "DMSO", "Talazoparib"));

	# the parental clone is sensitive
	pheno$resistance <- relevel(pheno$resistance, "Sensitive");

	# the parent clone is the reference
	pheno$clone <- factor(pheno$clone, clones);

	pheno$flowcell <- factor(gsub("_.*", "", as.character(pheno$batch)));

	pheno
}

read_msigdb <- function(collection, version="6.2") {
	release <- gsub(".", "", version, fixed=TRUE);
	qread(sprintf("~/data/msigdb/release-%s/%s.v%s.symbols.gmt", release, collection, version))
}

