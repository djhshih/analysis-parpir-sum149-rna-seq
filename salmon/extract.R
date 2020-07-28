library(io)

samples <- qread("../samples.vtr");
indir <- "cutadapt2"

out.fname <- filename("parpi-resist");

# These fields are based on the Gencode transcript fasta
fields <- c(
	"ensembl_transcript", "ensembl_gene",
	"havana_gene", "havana_transcript",
	"transcript_name", "gene_name", "transcript_length", "gene_type"
);

infiles <- unlist(lapply(samples, function(s) file.path(indir, s, "quant.sf")));
names(infiles) <- samples;

x <- lapply(infiles, qread, type="tsv");

annot <- x[[1]][, c("Name", "Length", "EffectiveLength")];

meta <- do.call(rbind, strsplit(as.character(annot$Name), "|", fixed=TRUE));
meta[meta == "-"] <- NA;
meta <- as.data.frame(meta, stringsAsFactors=FALSE);
names(meta) <- fields;

# ensure all row and column dimensions are equal across samples
stopifnot(unlist(lapply(x, nrow)) == nrow(x[[1]]))
stopifnot(unlist(lapply(x, ncol)) == ncol(x[[1]]))

tpm <- as.matrix(as.data.frame(lapply(x, function(xx) xx$TPM)));
colnames(tpm) <- samples;
rownames(tpm) <- meta[, 1];

counts <- as.matrix(as.data.frame(lapply(x, function(xx) xx$NumReads)));
colnames(counts) <- samples;
rownames(counts) <- meta[, 1];

stopifnot(meta$transcript_length == annot$Length)

meta$effective_length <- annot$EffectiveLength;

qwrite(meta, insert(out.fname, "fannot", ext="rds"));
qwrite(meta, insert(out.fname, "fannot", ext="tsv"));

qwrite(tpm, insert(out.fname, "tpm", ext="rds"));
qwrite(tpm, insert(out.fname, "tpm", ext="mtx"));

qwrite(counts, insert(out.fname, "counts", ext="rds"));
qwrite(counts, insert(out.fname, "counts", ext="mtx"));

