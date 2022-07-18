#!/bin/bash


# Input parameters: Directory with the projetc;
INDIR=$1;
OUTDIR=$2;
TRIMOMATIC=$3;


# run trimmomatic and PEAR
for file in $(ls $INDIR/*_R1_*.fastq.gz); do
  #raw read 1 file
  R1=$file;
  #raw read 2 file
  R2=${file/_R1_/_R2_};
  #quality controled read 1 file containing paired reads
  R1Paired=$OUTDIR$(echo $R1 | cut -d"/" -f3 | cut -d"." -f1)"_paired.fastq.gz"
  #quality controled read 2 file containing paired reads
  R2Paired=$OUTDIR$(echo $R2 | cut -d"/" -f3 | cut -d"." -f1)"_paired.fastq.gz"
  #quality controled read 1 file containing unpaired reads
  R1Unpaired=$OUTDIR$(echo $R1 | cut -d"/" -f3 | cut -d"." -f1)"_unpaired.fastq.gz"
  #quality controled read 2 file containing unpaired reads
  R2Unpaired=$OUTDIR$(echo $R2 | cut -d"/" -f3 | cut -d"." -f1)"_unpaired.fastq.gz"
  #prefix to save PEAR outputs
  outPEAR=$OUTDIR$(echo $R1 | cut -d"/" -f3 | cut -d"_" -f1);
  #trimmomatic code
  java -jar $TRIMOMATIC/trimmomatic-0.39.jar PE -threads 10 -phred33 $R1 $R2 $R1Paired $R1Unpaired $R2Paired $R2Unpaired ILLUMINACLIP:$TRIMOMATIC/adapters/SmarterThruPlex.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50;
  #PEAR code
  ~/programs/pear-0.9.11-linux-x86_64/bin/pear -f $R1Paired -r $R2Paired -o $outPEAR -j 10 -y 60G -n 50;
  #Converting fastq to fasta
  for fq in $(ls $outPEAR*.fastq); do
	  echo $fq;
	  fa=${fq/.fastq/.fasta};
	  sed -n '1~4s/^@/>/p;2~4p' $fq > $fa;
  done
done
