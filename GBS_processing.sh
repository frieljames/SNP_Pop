#!/bin/bash

dir_raw=$1
dir_pro=$2
dir_map=$3
adapter=$4
dir_ref=$5

if [ -z $dir_ref ]
   then

   echo "";
   echo "GBS_processing.sh script version 0.1 by James Friel";
   echo "USAGE: GBS_processing.sh <dir_raw_reads> <dir_processed_reads> <dir_mapped_reads> <adapter.fasta> <reference.fasta>"
   echo "";
   echo "Note: Remember to set the environmental variable BAMADDRG_BIN with the bamaddrg executable";
   echo "      example: export BAMADDRG_BIN=/my_path_to/bamaddrg";
   echo "";

   exit;
fi 

seed_cmd_bamaddrg=$BAMADDRG_BIN;

basename=$(ls $dir_raw| grep fq| sed -r 's/.R[12].fq.gz//'| uniq );

for base in $basename
        do
        echo "1- Running fastq-stats on $base"
        #fastq-stats $dir_raw"/"$base".R1.fq.gz" > $dir_raw"/"$base".R1.stats.txt";
        #fastq-stats $dir_raw"/"$base".R2.fq.gz" > $dir_raw"/"$base".R2.stats.txt";
        
        echo "2- Trimming raw file $base"
        #fastq-mcf -q30 -l50 -o $dir_pro"/"$base"_q30l50.R1.fq.gz" -o $dir_pro"/"$base"_q30l50.R2.fq.gz" $adapter $dir_raw"/"$base".R1.fq.gz" $dir_raw"/"$base".R2.fq.gz" 
        
        echo "3- Running fastq-stats on processed file $base"
        #fastq-stats $dir_pro/"$base"_q30l50.R1.fq.gz > $dir_pro/"$base"_q30l50.R1.stats.txt;
        #fastq-stats $dir_pro/"$base"_q30l50.R2.fq.gz > $dir_pro/"$base"_q30l50.R2.stats.txt;

        echo "4- Running the mapping on $base"
        #bwa mem -t60 $dir_ref $dir_pro"/"$base"_q30l50.R1.fq.gz" $dir_pro"/"$base"_q30l50.R2.fq.gz" | samtools view -F4 -bS -o $dir_map"/"$base".bwa.bam" - ;
        #samtools sort -o $dir_map"/"$base".bwa.bam" $dir_map"/"$base".bwa.bam"

        seed_cmd_bamaddrg=$seed_cmd_bamaddrg" -s "$base" -b "$dir_map"/"$base".bwa.bam"
done

bamaddrg_cmd=$seed_cmd_bamaddrg" > merged.bam";

echo "5- Running bamaddrg";

$bamaddrg_cmd
