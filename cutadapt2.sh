#!/bin/bash
# Second round of adapter sequence removal
# In this round, reads with remaining adapter sequences are removed.
# Minimum length is also imposed.

set -oe pipefail


outdir=cutadapt2/fastq
logdir=log/cutadapt2

cores=4

fastq_r1=$1
fname_r1=${fastq_r1##*/}
fstem_r1=${fname_r1%.cut.fastq.gz}

fastq_r2=$2
fname_r2=${fastq_r2##*/}
fstem_r2=${fname_r2%.cut.fastq.gz}

echo $fstem_r1
echo $fstem_r2

mkdir -p ${outdir} ${logdir}

# Look remaining 3' adapter in R1
# Illumina Multiplexing Adapter 1
# This occurs due to tandem duplication of the adapter
# We don't use the prefix A before the adapter sequence
# here, because it does not appear 50% of the time with
# the particular library prep used (KAPA Stranded mRNA seq Kit)
# The adapter also occurs with 5' truncation.
# So, we look also look for the sequence as a 5' unanchored adapter
# using the -b (both) flag.

# Additionally, we need to remove reads containing artificial sequence
# AATGATACGGCGACCACCGAGATCTACACG
# which partially matches the Illumina RNA PCR Primer

cutadapt \
	-j $cores \
	-b GATCGGAAGAGCACACGTCTGAACTCCAGTCA \
	-B CGTGTAGATCTCGGTGGTCGCCGTATCATT \
	--overlap 20 \
	--minimum-length 20 \
	-o ${outdir}/${fstem_r1}.cut2.fastq.gz \
	-p ${outdir}/${fstem_r2}.cut2.fastq.gz \
	${fastq_r1} ${fastq_r2} \
	| tee ${logdir}/${fstem_r1}.log

#	-a GATCGGAAGAGCACACGTCT \

