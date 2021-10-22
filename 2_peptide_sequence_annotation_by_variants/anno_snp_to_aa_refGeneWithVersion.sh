#!/bin/bash -l
#
#$ -cwd
#$ -j y
#$ -S /bin/bas
scriptDir=$1
fileName=$2
build=$3

annovar=${scriptDir}/annovar/

wd=`pwd`
echo $wd 

#step 1: run the ANNOVAR #

perl ${annovar}/annotate_variation.pl --downdb refGeneWithVer ${scriptDir}/annovar/humandb -build hg${build}

perl ${annovar}/annotate_variation.pl -out ${fileName}_hg${build}_refGeneVer -build hg${build} $fileName.bed ${annovar}/humandb/ -dbtype refGeneWithVer
perl ${annovar}/coding_change.pl -out ${fileName}_hg${build}_refGeneVer_aa.fa --newevf ${fileName}_hg${build}_refGeneVer_updated.exonic_variant_function --includesnp --onlyAltering --alltranscript ${fileName}_hg${build}_refGeneVer.exonic_variant_function ${annovar}/humandb/hg${build}_refGeneWithVer.txt ${annovar}/humandb/hg${build}_refGeneWithVerMrna.fa

#step 2.1: get annotation files#
awk -F "\t" '($2=="nonsynonymous SNV")||($2=="stoploss"){print $0}' ${fileName}_hg${build}_refGeneVer_updated.exonic_variant_function > ${fileName}_hg${build}_refGeneVer_updated.exonic_variant_function.nonsynStoploss
r1="${scriptDir}/lib/data_process.R"
cmd="R --slave --vanilla --args ${fileName}_hg${build}_refGeneVer_updated.exonic_variant_function.nonsynStoploss ${fileName}_hg${build}_refGeneVer_anno.txt ${fileName}_hg${build}_refGeneVer_anno_uniqNM.txt < $r1"
echo $cmd
eval $cmd

#step 2.2: get NM index file for the amino acid file created by the ANNOVA#
grep -n ">line" ${fileName}_hg${build}_refGeneVer_aa.fa > t1.txt
cut -d ":" -f 1 t1.txt > start.txt
cut -d ">" -f 2 t1.txt | cut -d " " -f 2,3 > nm.txt
paste nm.txt start.txt > tmp.txt
n=`wc -l ${fileName}_hg${build}_refGeneVer_aa.fa | awk '{print $1}'`
r2="${scriptDir}/lib/fasta_index.R"
cmd="R --slave --vanilla --args tmp.txt $n tmp2.txt < $r2"
echo $cmd
eval $cmd
grep "WILDTYPE" tmp2.txt > ${fileName}_hg${build}_refGeneVer_aa_fasta_NM_index.txt
rm -rf t1.txt
rm -rf start.txt
rm -rf nm.txt
rm -rf tmp.txt
rm -rf tmp2.txt

#step 3: get all possible peptides#
r="${scriptDir}/lib/protein_fasta_edit_2n.R"
n=`wc -l ${fileName}_hg${build}_refGeneVer_anno_uniqNM.txt | awk '{print $1}'`
for((i=2;i<=$n;i++))
do
j=$((i-1))
NM=`awk -F "\t" -v l=$i '(NR==l){print $1}' ${fileName}_hg${build}_refGeneVer_anno_uniqNM.txt`
gene=`awk -F "\t" -v NM=$NM '($16==NM){print $10}' ${fileName}_hg${build}_refGeneVer_anno.txt | awk -F "\t" '(NR==1){print $1}'`
s=`awk -F " " -v NM=$NM '($1==NM){print $3}' ${fileName}_hg${build}_refGeneVer_aa_fasta_NM_index.txt | awk -F " " '(NR==1){print $1}'`
e=`awk -F " " -v NM=$NM '($1==NM){print $4}' ${fileName}_hg${build}_refGeneVer_aa_fasta_NM_index.txt | awk -F " " '(NR==1){print $1}'`
awk -F "\t" -v NM=$NM '(NR==1)||($16==NM){print $0}' ${fileName}_hg${build}_refGeneVer_anno.txt > anno_${j}.txt
if (($s>=1))
then
echo -e ">${NM}_${gene}_WT |${gene} OS=Homo sapiens OX=9690 GN=${gene}" > ref_${j}.fasta
awk -v s=$s -v e=$e '(NR>s)&&(NR<=e){print $0}' ${fileName}_hg${build}_refGeneVer_aa.fa >> ref_${j}.fasta
cmd="R --slave --vanilla --args anno_${j}.txt ref_${j}.fasta alt_${j}.fasta < $r"
echo $cmd
eval $cmd
cat alt_${j}.fasta >> result_${fileName}_hg${build}_refGeneVer.fasta
fi
rm -rf anno_${j}.txt
rm -rf ref_${j}.fasta
rm -rf alt_${j}.fasta
done




