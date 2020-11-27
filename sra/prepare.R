library(io)
library(dplyr)

pheno <- qread("../annot/sample-info_parpi-resist_stage2.tsv");

r <- with(pheno,
	paste0(clone, "_r", gsub(".+-", "", sample_id))
);

biosample.d <- transmute(pheno,
	sample_name = sample_id,
	sample_title = paste("SUM149PT", clone, "clone,", treatment),
	organism = "Homo sapiens",
	isolate = r,
	age = 40,
	biomaterial_provider = "BioIVT",
	sex = "female",
	tissue = "mammary gland",
	cell_line = "SUM149PT",
	treatment = treatment
);

qwrite(biosample.d, "biosample-attributes.tsv")

sra.d <- transmute(pheno,
	sample_name = sample_id,
	library_ID = paste0("L", lane, "_", barcode),
	title = paste0("RNA-Seq of Homo sapiens: SUM149PT"),
	library_strategy = "RNA-Seq",
	library_source = "TRANSCRIPTOMIC",
	library_selection = "Oligo-dT",
	library_layout = "paired",
	platform = "ILLUMINA",
	instrument_model = "Illumina HiSeq 4000",
	design_description = "Paired-end mRNA-Seq",
	filetype = "fastq",
);

# pair fastq files together

files <- list_files("../cutadapt2/fastq", pattern=".+fastq.gz");
files <- sub("cut2.", "", files);

files.d <- data.frame(
	sample_id = sub("_.*", "", files),
	filename = files,
	read = sub(".*(R\\d).*", "\\1", files)
);

r1.d <- filter(files.d, read == "R1");
r2.d <- filter(files.d, read == "R2");

paired.d <- left_join(
	select(r1.d, -read),
	select(r2.d, -read),
	suffix = c("", "2"),
	by = "sample_id"
);

sra.d <- left_join(sra.d, paired.d, by = c(sample_name="sample_id"));

qwrite(sra.d, "sra-metadata.tsv");

