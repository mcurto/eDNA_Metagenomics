#!/bin/bash

SRA_file=$1
Out_dir=$2

#Please adapt the paths to the programs
TRIMOMATIC="/home/envmetagen/programs/Trimmomatic-0.39"
FASTQ_DUMP="~/programs/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump"
PEAR="~/programs/pear-0.9.11-linux-x86_64/bin/pear"


for acc in $(cat $SRA_file); do
  echo "Download fastq files for "$acc;

  #Define output files
  fq1=$acc"_1.fastq.gz";
  fq2=$acc"_2.fastq.gz";
  fq1_paired=$Out_dir"/"$acc"_1.paired.fastq.gz";
  fq2_paired=$Out_dir"/"$acc"_2.paired.fastq.gz";
  fq1_unpaired=$Out_dir"/"$acc"_1.unpaired.fastq.gz";
  fq2_unpaired=$Out_dir"/"$acc"_2.unpaired.fastq.gz";
  fastq1_1000=$Out_dir"/"$acc"_1.paired.1000.fastq";
  fastq2_1000=$Out_dir"/"$acc"_2.paired.1000.fastq";
  assembled=$Out_dir"/"$acc".paired.1000.assembled.fastq";
  un_forward=$Out_dir"/"$acc".paired.1000.unassembled.forward.fastq";
  un_reverse=$Out_dir"/"$acc".paired.1000.unassembled.reverse.fastq";
  fasta_ass=$Out_dir"/"$acc".paired.1000.assembled.fasta";
  fasta_for=$Out_dir"/"$acc".paired.1000.unassembled.forward.fasta";
  fasta_rev=$Out_dir"/"$acc".paired.1000.unassembled.reverse.fasta";

  #Download fastq
  $FASTQ_DUMP -X 1000000 --gzip  --split-3 $acc;
  #Clear SRA temporary files
  rm /home/envmetagen/data/SRA_data/sra/*;

  #Quality control fastq
  java -jar $TRIMOMATIC/trimmomatic-0.39.jar PE -threads 10 -phred33 $fq1 $fq2 $fq1_paired $fq1_unpaired $fq2_paired $fq2_unpaired ILLUMINACLIP:$TRIMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50;
  #Exclude unpaired files
  rm $Out_dir"/"*.unpaired.fastq.gz;

  #Subsample 1000 reads
  seqtk sample -s100 $fq1_paired 1000 > $fastq1_1000;
  seqtk sample -s100 $fq2_paired 1000 > $fastq2_1000;

  #Merge reads 1 and reads 2
  $PEAR -f $fastq1_1000  -r $fastq2_1000  -o $Out_dir"/"$acc".paired.1000";

  #Convert fastq to fasta
  sed -n '1~4s/^@/>/p;2~4p' $assembled > $fasta_ass;
  sed -n '1~4s/^@/>/p;2~4p' $un_forward > $fasta_for;
  sed -n '1~4s/^@/>/p;2~4p' $un_reverse > $fasta_rev;

  #Clean non necessary files
  rm $fq1 $fq2 $fq1_paired $fq2_paired $fastq1_1000 $fastq2_1000 $assembled $un_forward $un_reverse $Out_dir/*.fastq;
done
