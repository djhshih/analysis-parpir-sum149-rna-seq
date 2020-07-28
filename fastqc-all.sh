#!/bin/bash

set -oe pipefail

mkdir -p fastqc
fastqc fastq/* -o fastqc

