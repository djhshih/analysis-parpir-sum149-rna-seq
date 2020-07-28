#!/bin/bash

set -oe pipefail

mkdir -p cutadapt/fastqc
fastqc cutadapt/fastq/* -o cutadapt/fastqc

