# DMSO is the reference
pheno$treatment <- relevel(pheno$treatment, "DMSO");

# the parental clone is sensitive
pheno$resistance <- relevel(pheno$resistance, "Sensitive");

# the parent clone is the reference
pheno$clone <- relevel(pheno$clone, "Parental");

rownames(fannot) <- fannot$ensembl_transcript;
pheno$flowcell <- factor(gsub("_.*", "", as.character(pheno$batch)));
