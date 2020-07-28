#!/bin/bash

indir=$1
sample=$2

outdir=salmon/cutadapt2/${sample}

index=~/data/gencode/release-28/salmon_quasi

fastq_r1=${indir}/${sample}*_R1*.fastq.gz
fastq_r2=${indir}/${sample}*_R2*.fastq.gz


mkdir -p $outdir

# library is pair-end (_i_nward oriented), _s_tranded, and first read is _r_eversed
# therefore, library type is ISR

salmon quant \
	-i $index \
	-l ISR -1 $fastq_r1 -2 $fastq_r2 \
	-o $outdir
