#!/bin/bash

#list=samples.vtr
list=samples_remnant.vtr

script=./cutadapt.sh
indir=fastq

for x in $(cat $list); do
	echo $x;
	echo "$script ${indir}/$x*"
	$script ${indir}/$x*
done
