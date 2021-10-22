rm(list=ls())
t.args<-commandArgs();i.args<-which(t.args=="--args")
inFile<-t.args[i.args+1]
outA<-t.args[i.args+2]
outB<-t.args[i.args+3]
cat("\n\n")
print(paste(t.args,collapse=" "))
cat("\n\n")

library(stringr)
library(readr)
a<-read.table(inFile,sep="\t",header=F)
dim(a)
colnames(a)<-c("lineNumber","SNV_type","anno","chr","start","end","ref","alt")
head(a)
a$anno_number<-str_count(a$anno,",")

a1<-a[a$anno_number==1,]
a2<-a[a$anno_number>1,]
a3<-NULL
for(i in c(1:nrow(a2))){
    tmp<-a2$anno[i]
    tmp2<-as.data.frame(str_split(tmp,","))
    n<-a2$anno_number[i]
    lineNumber.i<-a2$lineNumber[i]
    SNV_type.i<-a2$SNV_type[i]
    chr.i<-a2$chr[i]
    start.i<-a2$start[i]
    end.i<-a2$end[i]
    ref.i<-a2$ref[i]
    alt.i<-a2$alt[i]
    anno_number.i<-a2$anno_number[i]
    tmp3<-NULL
    for(j in c(1:n)){
        tmp3.j<-NULL
        anno.j<-tmp2[j,1]
        tmp3.j<-cbind.data.frame(lineNumber.i,SNV_type.i,anno.j,chr.i,start.i,end.i,ref.i,alt.i,anno_number.i)
        tmp3<-rbind.data.frame(tmp3,tmp3.j)
    }
    a3.i<-tmp3
    a3<-rbind.data.frame(a3,a3.i)
}
dim(a3)
colnames(a3)<-c("lineNumber","SNV_type","anno","chr","start","end","ref","alt","anno_number")
final<-rbind.data.frame(a1,a3)
dim(final)

tmp4<-as.data.frame(t(as.data.frame(strsplit(as.character(final$anno),split="[,,:,.]"))))
tmp5<-tmp4[c(1:4,6,8)]
dim(tmp5)
head(tmp5)
colnames(tmp5)<-c("gene","NM","NM_version","exon_number","dna_substitute","aa_substitute")
tmp5$NM_wVer<-paste(tmp5$NM,tmp5$NM_version,sep=".")
head(tmp5)
row.names(tmp5)<-c(1:nrow(tmp5))
head(tmp5)
final2<-cbind.data.frame(final,tmp5)
dim(final2)
head(final2)

write.table(final2,"final2.txt",col.names=T,row.names=F,quote=F,sep="\t")
final2<-read.table("final2.txt",header=T,sep="\t",stringsAsFactors = F)

tmp<-final2
final2$dna_ref<-str_sub(final2$dna_substitute,1,1)
final2$dna_alt<-str_sub(final2$dna_substitute,-1,-1)
final2$aa_ref<-str_sub(final2$aa_substitute,1,1)
final2$aa_alt<-str_sub(final2$aa_substitute,-1,-1)
final2$dna_substitute_pos<-parse_number(as.character(final2$dna_substitute))
final2$aa_substitute_pos<-parse_number(as.character(final2$aa_substitute))
dim(final2)
head(final2)
final2<-unique(final2)
dim(final2)
#remove those sites with * for either aa_ref or aa_alt#
final3<-final2[! final2$aa_ref %in% c("*"),]
final4<-final3[! final3$aa_alt %in% c("*"),]
#remove those I->L or L->I mutation#
final4$aa_change<-paste(final4$aa_ref,final4$aa_alt,sep="")
final5<-final4[! final4$aa_change %in% c("IL","LI"),]
n<-ncol(final3)
final6<-final5[c(1:n)]
dim(final6)

write.table(final6,outA,col.names=T,row.names=F,quote=F,sep="\t")

nm_list<-as.data.frame(unique(final6$NM_wVer))
colnames(nm_list)<-c("NM_withVer")
dim(nm_list)
head(nm_list)
write.table(nm_list,outB,col.names=T,row.names=F,quote=F,sep="\t")



