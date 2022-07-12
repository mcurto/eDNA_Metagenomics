#!/bin/bash


# Input parameters: Directory with the projetc;
INDIR=$1;
OUTDIR=$2;


# run trimmomatic and PEAR
for file in $(ls $INDIR/*_R1_*.fastq.gz); do
  R1=$file;
  R2=${file/_R1_/_R2_};
  R1Paired=$OUTDIR$(echo $R1 | cut -d"/" -f3 | cut -d"." -f1)"_paired.fastq.gz"
  R2Paired=$OUTDIR$(echo $R2 | cut -d"/" -f3 | cut -d"." -f1)"_paired.fastq.gz"
  R1Unpaired=$OUTDIR$(echo $R1 | cut -d"/" -f3 | cut -d"." -f1)"_unpaired.fastq.gz"
  R2Unpaired=$OUTDIR$(echo $R2 | cut -d"/" -f3 | cut -d"." -f1)"_unpaired.fastq.gz"
  outPEAR=$OUTDIR$(echo $R1 | cut -d"/" -f3 | cut -d"_" -f1);
  java -jar /home/dsi/programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 10 -phred33 $R1 $R2 $R1Paired $R1Unpaired $R2Paired $R2Unpaired ILLUMINACLIP:/home/dsi/programs/Trimmomatic-0.39/adapters/SmarterThruPlex.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50;
  /home/dsi/programs/pear-0.9.11-linux-x86_64/bin/pear -f $R1Paired -r $R2Paired -o $outPEAR -j 10 -y 60G -n 50;
done
