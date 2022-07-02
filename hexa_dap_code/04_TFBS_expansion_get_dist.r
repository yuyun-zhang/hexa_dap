file<-commandArgs(T)

max<-as.matrix(read.table(file,sep = "\t"))
rownames(max)<-paste0("V",1:nrow(max))
table<-na.omit(reshape2::melt(max))
write.table(table,paste0(file,".table"),quote=F,row.names=F,col.names=F)


