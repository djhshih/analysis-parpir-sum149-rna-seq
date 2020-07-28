#!/bin/bash

list=samples.vtr

script=./salmon.sh
indir=cutadapt2/fastq
outdir=salmon/cutadapt2

for sample in $(cat $list); do
	echo $sample;
	$script $sample $indir $outdir
done
