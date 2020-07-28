#!/bin/bash

set -oe pipefail

mkdir -p cutadapt2/fastqc
fastqc cutadapt2/fastq/* -o cutadapt2/fastqc

