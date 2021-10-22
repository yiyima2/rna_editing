rm(list=ls())
t.args<-commandArgs();i.args<-which(t.args=="--args")
annoFile<-t.args[i.args+1]
refFasta<-t.args[i.args+2]
outFasta<-t.args[i.args+3]
cat("\n\n")
print(paste(t.args,collapse=" "))
cat("\n\n")

library(seqinr)
anno<-read.table(annoFile,sep="\t",quote="\n",comment.char="",header=T, stringsAsFactors = F)
ref<-read.fasta(refFasta,seqtype = c("AA"))
anno$aa_ref[anno$aa_ref %in% c("TRUE")]<-"T"
anno$aa_ref[anno$aa_ref %in% c("FALSE")]<-"F"
anno$aa_alt[anno$aa_alt %in% c("TRUE")]<-"T"
anno$aa_alt[anno$aa_alt %in% c("FALSE")]<-"F"

np<-anno$NM_wVer[1]
gene<-anno$gene[1]
name1<-paste(np,gene,sep="_")
n<-nrow(anno)
name2<-strsplit(attr(ref[[1]],"Annot"),"\\|")[[1]][2]

if(n>1){
a<-anno$aa_substitute_pos
b<-list()
for(i in c(1:length(a))){
    iter=t(combn(a,i))
    b[[i]] <- iter
}

e<-anno$aa_alt
f<-list()
for(i in c(1:length(e))){
    iter=t(combn(e,i))
    f[[i]] <- iter
}

g<-anno$aa_substitute
h<-list()
for(i in c(1:length(g))){
    iter=t(combn(g,i))
    h[[i]] <- iter
}

final<-ref
for(i in c(1:length(b))){
    for(j in c(1:nrow(b[[i]]))){
        altName<-paste(name1,".",i,"var.tot",n,"var",sep="")
        alt<-ref
        for(k in c(1:i)){
            altName<-paste(altName,h[[i]][j,k],sep="_")
            altAnno<-paste(altName,name2,sep=" |")
            names(alt)<-altName
            attr(alt[[1]],"name")<-altName
            attr(alt[[1]],"Annot")<-altAnno
            alt[[1]][b[[i]][j,k]]<-as.character(f[[i]][j,k])
        }
        final<-c(final,alt)
    }
}
write.fasta(sequences=final,names=paste(names(final),name2,sep=" |"),file.out=outFasta)
}

if(n==1){
    final<-ref
    alt<-ref
    b<-anno$aa_substitute_pos[1]
    substitute<-anno$aa_substitute[1]
    alt_aa<-anno$aa_alt[1]
    altName<-paste(name1,".1var.tot1var_",substitute,sep="")
    altAnno<-paste(altName,name2,sep=" |")
    names(alt)<-altName
    attr(alt[[1]],"name")<-altName
    attr(alt[[1]],"Annot")<-altAnno
    alt[[1]][b]<-as.character(alt_aa)
    final<-c(final,alt)
    write.fasta(sequences=final,names=paste(names(final),name2,sep=" |"),file.out=outFasta)
}

