file<-commandArgs(T)
data<-read.table(file,row.names=1,header=T)
data<-abs(data)
ss<-data.frame(apply(data,1,sum))
names(ss)<-"divergence"
write.table(ss,paste0(file,".abssum"),row.names=T,col.names=T,quote=F)
