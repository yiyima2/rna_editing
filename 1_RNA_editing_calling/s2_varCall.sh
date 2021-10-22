#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -j y
source /mnt/mfs/cluster/bin/HgrcPathSetup.sh

input=$1
dir=$2

dbsnp="${dir}/dbsnp_138.hg38.vcf.gz"
ref="${dir}/GRCh38_Gencode24/GRCh38.primary_assembly.genome.fa"

picard="${dir}/picard.jar"
gatk36="${dir}/gatk-3.6/GenomeAnalysisTK.jar"

echo "variant calling"
cmd2="java -Djava.io.tmpdir=${dir} -jar $gatk36 -T HaplotypeCaller -R $ref -I ${input}_splitted.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 --dbsnp $dbsnp -o ${input}.vcf"
echo $cmd2
eval $cmd2









