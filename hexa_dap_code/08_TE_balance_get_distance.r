library(rdist)

file<-commandArgs(T)
nor<-read.table(file,row.names=1)
ref<-read.table("ref.dist",row.names=1,header=T)
dis<-data.frame(cdist(X=nor,Y=ref))
nor<-cbind(nor,dis$X7)
colnames(nor)[4]="dist"
write.table(nor,paste0(file,".dist"),quote=F,row.names=T,col.names=F)



#################
##### ref.dist ##
#################
type A B D
A-dom 1 0 0
B-dom 0 1 0
D-dom 0 0 1
A-sup 0 0.5 0.5
B-sup 0.5 0 0.5
D-sup 0.5 0.5 0
balanced 0.33 0.33 0.33



