rm(list=ls())
t.args<-commandArgs();i.args<-which(t.args=="--args")
inFile<-t.args[i.args+1]
N<-t.args[i.args+2]
outFile<-t.args[i.args+3]
cat("\n\n")
print(paste(t.args,collapse=" "))
cat("\n\n")

a<-read.table(inFile,sep="\t",header=F)
b<-a$V2-1
a$end<-c(b[2:length(b)],N)
colnames(a)<-c("NM_withVer","start_line","end_line")
head(a)
write.table(a,outFile,col.names=T,row.names=F,quote=F,sep="\t")
