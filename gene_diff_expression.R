#3.1
#Creation of an element that contains the data
GED_norm<-read.delim("GeneExpressionDataset_normalized.tsv",row.names=1)
#Creation of a matrix where the logged differential expressions will be stored
logmat<-matrix(NA,length(GED_norm[,1]),5)
row.names(logmat)<-row.names(GED_norm)
colnames(logmat)<-c("log2FC(TG)","log2FC(TherA)","log2FC(TherB)","log2FC(TherC)","log2FC(TherD)")
#Application of a while-loop and a for-loop to compute the logged differential expressions
k=11
l=20
m=1
while (l<=60){
  for(i in 1:length(logmat[,1])){
    logmat[i,m]<-mean(as.numeric(GED_norm[i,k:l]))-mean(as.numeric(GED_norm[i,1:10]))
  }
  k<-k+10
  l<-l+10
  m<-m+1}
print(head(logmat))

#3.2
library(stats)
#Creation of a matrix where the p-values will be stored
pmat<-matrix(NA,length(GED_norm[,1]),5)
row.names(pmat)<-row.names(GED_norm)
colnames(pmat)<-c("TG","TherA","TherB","TherC","TherD")
#Application of a while-loop and a for-loop in order to perform the t-test and compute the p-values
o=11
p=20
r=1
while (p<=60){
  for(i in 1:length(pmat[,1])){
    pmat[i,r]<-t.test(as.numeric(GED_norm[i,o:p]),y=(as.numeric(GED_norm[i,1:10])))$p.value
  }
  o<-o+10
  p<-p+10
  r<-r+1
}
print(head(pmat))

#Creation of a matrix where the adjusted p-values (q-values) will be stored
qmat<-matrix(NA,length(pmat[,1]),5)
row.names(qmat)<-row.names(pmat)
colnames(qmat)<-colnames(pmat)
#Application of a for loop in order to compute the q-values
for (j in 1:5){
   qmat[,j]<-p.adjust(pmat[,j],method="BH")
}
print(head(qmat))

#3.3
#Creation of a matrix where the genes' qualitative expression status will be stored
genemat<-matrix(NA,length(GED_norm[,1]),5)
row.names(genemat)<-row.names(qmat)
colnames(genemat)<-colnames(qmat)
#Application of a double for-loop and an if-loop in order to compute the genes' statuses
for (j in 1:5){
  for (i in 1:length(genemat[,1])){
    if ((abs(logmat[i,j])>1)&(qmat[i,j])<=0.05){
      genemat[i,j]<-T
    }
    else{
      genemat[i,j]<-F
    }
  }
}
print(head(genemat))

#3.4
#Comparison of the four medicines' intervention to the genes' expression
lessmat<-matrix( c(min(table(genemat[,2])),min(table(genemat[,3])),min(table(genemat[,4])),min(table(genemat[,5]))),1,4)
colnames(lessmat)<-c("Therapy A","Therapy B","Therapy C","Therapy D")
print(lessmat)
barplot(lessmat,main="Therapy Intervention Levels", ylab="Differentially Expressed Genes(#)")