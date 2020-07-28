library(io)

before <- qread("no-cutadapt/10-C20-V-1/quant.sf", type="tsv");
after <- qread("cutadapt2/10-C20-V-1/quant.sf", type="tsv");

dim(before)
dim(after)

stopifnot(as.character(before$Name) == as.character(after$Name))

x <- log2(before$TPM + 1);
y <- log2(after$TPM + 1);

print(cor(x, y))

smoothScatter(x, y)

summary(x)
summary(y)

features <- as.character(before$Name);
hugos <- unlist(lapply(strsplit(features, "|", fixed=TRUE), function(x) x[6]));

d <- x - y;
names(d) <- hugos;

summary(abs(d))

sort(d, decreasing=TRUE)[1:50]

