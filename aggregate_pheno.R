library(io)
library(dplyr)

pheno.fname <- "sample-info_parpi-resist.tsv";

pheno <- qread(pheno.fname);
batch <- qread("batch-info_parpi-resist.tsv");

pheno <- full_join(pheno, batch, by="sample_id")

qwrite(pheno, tag(pheno.fname, "stage2"));
