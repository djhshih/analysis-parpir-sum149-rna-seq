library(io)
library(sva)

source("../R/preamble.R")

in.fname <- as.filename("parpi-resist_ltpm-genes-max.rds");
out.fname <- insert(in.fname, "cb");

x <- qread(in.fname);
fannot <- qread("parpi-resist_fannot.rds");
pheno <- setup_pheno(qread("../sample-info_parpi-resist_stage2.tsv"));

stopifnot(colnames(x) == as.character(pheno$sample_id))

#mod <- model.matrix(~ treatment + clone, data=pheno)
#x.cb <- ComBat(x, pheno$batch, mod, par.prior=FALSE);   # generates all NA!

#mod <- model.matrix(~ clone, data=pheno)

mod <- model.matrix(~ 1, data=pheno)
x.cb <- ComBat(x, pheno$batch);
# fails!

qwrite(x.cb, out.fname);

