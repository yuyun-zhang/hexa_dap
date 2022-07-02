args <- commandArgs(T)
chr  <- args[1]
library(dplyr)
library(data.table)
dat <- fread(paste(chr,".merge.bed",sep=""))
id <- fread(paste(chr,".merge.id",sep=""),head=F)
dat <- left_join(dat,id,by=c("V4"="V1"))
dat <- dat[,c(1:3,7,5,6)]
write.table(dat,paste(chr,".merge.1.bed",sep=""),quote=F,sep="\t",col.names=F,row.names=F)
