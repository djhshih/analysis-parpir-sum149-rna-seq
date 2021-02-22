#!/bin/bash

threads=4
genome_dir=~/data/gencode/grch38/star/sjdbOverhang75

indir=cutadapt2/fastq

#fastq_r1=$indir/1-Parental-1_S64_L008_R1_001.cut2.fastq.gz
#fastq_r2=$indir/1-Parental-1_S64_L008_R2_001.cut2.fastq.gz

#fastq_r1=$indir/R450-26-C2-V-2_S10_L001_R1_001.cut2.fastq.gz
#fastq_r2=$indir/R450-26-C2-V-2_S10_L001_R2_001.cut2.fastq.gz

#fastq_r1=$indir/R450-14-C17-V-2_S1_L001_R1_001.cut2.fastq.gz
#fastq_r2=$indir/R450-14-C17-V-2_S1_L001_R2_001.cut2.fastq.gz

#fastq_r1=$indir/R450-46-C19-P-1_S27_L003_R1_001.cut2.fastq.gz
#fastq_r2=$indir/R450-46-C19-P-1_S27_L003_R2_001.cut2.fastq.gz

#fastq_r1=$indir/R450-34-C20-P-1_S20_L002_R1_001.cut2.fastq.gz
#fastq_r2=$indir/R450-34-C20-P-1_S20_L002_R2_001.cut2.fastq.gz

#fastq_r1=$indir/R450-19-C26-V-1_S11_L001_R1_001.cut2.fastq.gz
#fastq_r2=$indir/R450-19-C26-V-1_S11_L001_R2_001.cut2.fastq.gz

#fastq_r1=$indir/R450-40-C29-P-1_S29_L003_R1_001.cut2.fastq.gz
#fastq_r2=$indir/R450-40-C29-P-1_S29_L003_R2_001.cut2.fastq.gz

fastq_r1=$indir/7-C30-V-1_S62_L008_R1_001.cut2.fastq.gz
fastq_r2=$indir/7-C30-V-1_S62_L008_R2_001.cut2.fastq.gz



fname=${fastq_r1##*/}
sample=${fname%%_*}

out_path=bam/$sample

mkdir -p $out_path

# use ENCODE standard options for (long) RNA-seq pipeline for insert size ~200 bp
STAR \
	--runThreadN $threads \
	--genomeDir $genome_dir \
	--readFilesCommand gunzip -c \
	--readFilesIn $fastq_r1 $fastq_r2 \
	--outFilterType BySJout \
	--outFilterMultimapNmax 20 \
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--outFilterMismatchNmax 999 \
	--outFilterMismatchNoverReadLmax 0.04 \
	--alignIntronMin 20 \
	--alignIntronMax 1000000 \
	--alignMatesGapMax 1000000 \
	--outFileNamePrefix $out_path/ \
	--outSAMtype BAM SortedByCoordinate \
  --twopassMode Basic \
  --chimSegmentMin 12 \
  --chimJunctionOverhangMin 12 \
  --chimSegmentReadGapMax 3 \
  --alignSJstitchMismatchNmax 5 -1 5 5 \
  --outSAMstrandField intronMotif \
  --chimOutJunctionFormat 1 \

# conflicting recommendation:
#  --alignSJDBoverhangMin 10
#  --alignIntronMax 100000
#  --alignMatesGapMax 100000

samtools index $out_path/*.bam


