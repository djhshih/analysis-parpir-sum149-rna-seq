#!/bin/bash

list=samples.vtr

script=bin/cutadapt.sh
indir=fastq

for x in $(cat $list); do
	echo $x;
	echo "$script ${indir}/$x*"
	$script ${indir}/$x*
done
