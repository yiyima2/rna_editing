#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -j y

dir=$1
input=$2
ID=$3
PL=$4
PU=$5

ref="${dir}/GRCh38.primary_assembly.genome.fa"

picard="${dir}/picard.jar"
gatk35="${dir}/gatk-3.5/GenomeAnalysisTK.jar"

echo $ref

echo "STEP 1: add read groups, sort, mark duplicates, and create index"
java -jar $picard AddOrReplaceReadGroups \
I=${dir}/${input}.bam \
O=${input}_rgAddedSorted.bam \
SO=coordinate \
RGID=${ID} \
RGLB=${ID} \
RGPL=${PL} \
RGPU=${PU} \
RGSM=${ID} 

echo "STEP 2: markDuplicates"
java -jar $picard MarkDuplicates \
I=${input}_rgAddedSorted.bam \
O=${input}_dedupped.bam \
CREATE_INDEX=TRUE \
VALIDATION_STRINGENCY=SILENT \
M=output.metrics

echo "STEP 3: split'N'Trim and reassign mapping qualities"
cmd="java -jar $gatk35 -T SplitNCigarReads -R $ref -I ${input}_dedupped.bam -o ${input}_splitted.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 --unsafe ALLOW_N_CIGAR_READS"
echo $cmd
eval $cmd


rm -rf ${input}_rgAddedSorted.bam
rm -rf ${input}_rgAddedSorted.bai
rm -rf ${input}_dedupped.bam
rm -rf ${input}_dedupped.bai









