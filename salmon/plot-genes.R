library(io)
library(ggplot2)
library(reshape2)
library(dplyr)

source("../R/preamble.R");

# cell line: SUM149

out.fname <- filename("parpi-resist");
pdf.fname <- insert(out.fname, ext="pdf");

x <- qread("parpi-resist_counts.rds");
fannot <- qread("parpi-resist_fannot.rds");
pheno <- setup_pheno(qread("../annot/sample-info_parpi-resist_stage2.tsv"));

stopifnot(colnames(x) == as.character(pheno$sample_id))

counts <- round(x);
mode(counts) <- "integer";


plot_gene_counts <- function(gene, transcripts=NULL) {
	if (is.null(transcripts)) {
		transcripts <- fannot$ensembl_transcript[fannot$gene_name == gene];
	}

	counts.s.m <- melt(counts[transcripts, , drop=FALSE] + 1);
	colnames(counts.s.m) <- c("transcript", "sample_id", "count");

	d <- left_join(pheno, counts.s.m, by="sample_id");
	ggplot(d, aes(x=clone, y=count, colour=treatment)) + 
		geom_jitter(width=0.1, alpha=0.35) + theme_bw() + 
		facet_wrap(~ transcript) + 
		#scale_y_log10() + ylab("count + 1") +
		ylab("count") +
		theme(axis.text.x = element_text(angle=45, hjust=1)) 
}


qdraw(
	plot_gene_counts("BRCA1") +
		theme(strip.background = element_blank()) +
		facet_wrap(~ transcript, ncol=5)
	,
	width = 10, height = 10,
	file = insert(pdf.fname, c("brca1", "all-transcripts"))
);

qdraw(
	plot_gene_counts("BRCA1", c("ENST00000352993.7", "ENST00000357654.8")) +
		theme(strip.background = element_blank()) +
		facet_wrap(~ transcript, ncol=5)
	,
	file = insert(pdf.fname, c("brca1", "main-transcripts"))
);

qdraw(
	plot_gene_counts("BRCA2") +
		theme(strip.background = element_blank()) +
		facet_wrap(~ transcript, ncol=5)
	,
	width = 10, height = 10,
	file = insert(pdf.fname, c("brca2", "all-transcripts"))
);

qdraw(
	plot_gene_counts("TP53BP1") +
		theme(strip.background = element_blank()) +
		facet_wrap(~ transcript, ncol=5)
	,
	width = 10, height = 10,
	file = insert(pdf.fname, c("tp53bp1", "all-transcripts"))
);

plot_gene_counts("CHEK1")
plot_gene_counts("TP53")

plot_gene_counts("CDKN1A")
plot_gene_counts("CDKN1A") + scale_y_continuous()

plot_gene_counts("CHEK1") + scale_y_continuous()
plot_gene_counts("CHEK2") + scale_y_continuous()

plot_gene_counts("MMP14")
plot_gene_counts("GREM1")
plot_gene_counts("MSX1")
plot_gene_counts("TGM2")

plot_gene_counts("PARP1")

plot_gene_counts("BRCA2")

plot_gene_counts("GABBR2")

plot_gene_counts("HIST1H1E")
plot_gene_counts("HIST1H2BH")

plot_gene_counts("KRT1")
plot_gene_counts("KRT16")

plot_gene_counts("ATR")
plot_gene_counts("ATM")

