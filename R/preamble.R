library(ggsci)
library(RColorBrewer)

shift <- function(x, n=1) {
	if (n == 0) {
		x
	} else {
		c(tail(x, -n), head(x, n))
	}
}

diff.cols <- rev(brewer.pal(9, "RdBu"));
diff.group.cols <- c(NS = "#333333", Down = diff.cols[1], Up = diff.cols[9]);

clones.from <- c("Parental", "C2", "C17", "C19", "C20", "C26", "C29", "C30");
clones <- c("P", "RC1", "RC2", "RC3", "RC4", "RC5", "RC6", "RC7");

rename_clones <- function(x) {
	as.character(factor(x, levels=clones.from, labels=clones))
}

clone.cols <- shift(pal_d3()(length(clones)), -1);
names(clone.cols) <- clones;

setup_pheno <- function(pheno, rename.clones=FALSE) {
	# DMSO is the reference
	pheno$treatment = factor(pheno$treatment, levels=c("None", "DMSO", "Talazoparib"));

	# the parental clone is sensitive
	pheno$resistance <- relevel(pheno$resistance, "Sensitive");

	# the parent clone is the reference
	if (rename.clones) {
		pheno$clone <- factor(pheno$clone, levels=clones.from, labels=clones);
	} else {
		pheno$clone <- factor(pheno$clone, levels=clones.from);
	}

	pheno$flowcell <- factor(gsub("_.*", "", as.character(pheno$batch)));

	pheno
}

read_msigdb <- function(collection, version="6.2") {
	release <- gsub(".", "", version, fixed=TRUE);
	qread(sprintf("~/data/msigdb/release-%s/%s.v%s.symbols.gmt", release, collection, version))
}

capitalize <- function(s) {
	ifelse(is.na(s),
		NA, 
		paste0(toupper(substring(s, 1, 1)), substring(s, 2))	
	)
}

rename_hallmarks <- function(x, ...) {
	UseMethod("rename_hallmarks", x)
}

rename_hallmarks.default <- function(x, ...) {
	x	<- tolower(gsub("_", " ", sub("HALLMARK_", "", x)));

	library(magrittr)
	x %<>% gsub("kras ", "KRAS ", .) %>%
		gsub("nfkb", "NFKB", .) %>%
		gsub("e2f ", "E2F ", .) %>%
		gsub("tnfa ", "TNF-alpha ", .) %>%
		gsub("myc ", "MYC ", .) %>%
		gsub("il(\\d) ", "IL\\1 ", .) %>%
		gsub("stat(\\d) ", "STAT\\1 ", .) %>%
		gsub("mtorc", "mTORC", .) %>%
		gsub("mtor", "mTOR", .) %>%
		gsub("pi3k ", "PI3K ", .) %>%
		gsub("beta catenin", "beta-catenin", .) %>%
		gsub("tgf beta", "TGF-beta", .) %>%
		gsub("akt ", "Akt ", .) %>%
		gsub("jak ", "Jak ", .) %>%
		gsub("wnt ", "Wnt ", .) %>%
		gsub("notch ", "Notch ", .) %>%
		gsub("uv ", "UV ", .) %>%
		gsub("g2m ", "G2M ", .) %>%
		gsub("dna ", "DNA ", .);
	x <- capitalize(x) %>%
		gsub("P53", "p53", .) %>%
		gsub("MTORC", "mTORC", .);

	x
}

rename_hallmarks.matrix <-
rename_hallmarks.data.frame <- function(x, ...) {
	rownames(x) <- rename_hallmarks(rownames(x));
	x
}

revlog_trans <- function(base = exp(1)) {
	library(scales)
	trans <- function(x) -log(x, base)
	inv <- function(x) base^(-x)
	trans_new(paste0("reverselog-", format(base)), trans, inv, 
						log_breaks(base = base), 
						domain = c(1e-100, Inf))
}

theme_clean <- function() {
	library(ggplot2)

	theme_bw() +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank()
	)
}

