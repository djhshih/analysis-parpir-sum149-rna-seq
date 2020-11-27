#!/bin/bash

list=annot/samples.vtr

script=bin/salmon.sh
indir=cutadapt2/fastq
outdir=salmon/cutadapt2

for sample in $(cat $list); do
	echo $sample;
	$script $sample $indir $outdir
done
