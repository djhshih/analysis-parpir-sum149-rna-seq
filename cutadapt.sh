#!/bin/bash

set -oe pipefail


outdir=cutadapt/fastq
logdir=log/cutadapt

cores=4

fastq_r1=$1
fname_r1=${fastq_r1##*/}
fstem_r1=${fname_r1%.fastq.gz}

fastq_r2=$2
fname_r2=${fastq_r2##*/}
fstem_r2=${fname_r2%.fastq.gz}

echo $fstem_r1
echo $fstem_r2

mkdir -p ${outdir} ${logdir}

# Trim 3' adapter sequences:
#  Illumina Multiplexing Adapter 1 in R1 reads
#  Reverse complement of Illumina Multiplexing Adapter 2 without
#  the last nucleotide in the adapter in R2 reads

# Illumina Multiplexing Adapter 2
# ACACTCTTTCCCTACACGACGCTCTTCCGATCT
# Reverse complement
#	AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

# A prefix of 'A' is added since the adapter often occurs
# with it (but only ~50% with our particular library prep...)
# We will follow-up with another round of filter without the prefix.

cutadapt \
	-j $cores \
	-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	-o ${outdir}/${fstem_r1}.cut.fastq.gz \
	-p ${outdir}/${fstem_r2}.cut.fastq.gz \
	${fastq_r1} ${fastq_r2} \
	| tee ${logdir}/${fstem_r1}.log

#	-a GATCGGAAGAGCACACGTCT \
#	-A GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \

