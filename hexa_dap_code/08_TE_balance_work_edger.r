library(edgeR)
info<-commandArgs(T)
exprSet<- read.table(file = info[1], header = TRUE, row.names = 1, check.names=F)
group_list <- factor(c(rep("Control",info[2]),rep("Treat",info[3])))
exprSet <- DGEList(counts = exprSet, group = group_list)
bcv = 0.2  
et <- exactTest(exprSet, dispersion=bcv^2)
d <- calcNormFactors(exprSet)
cps <- cpm(d, normalized.lib.sizes=TRUE)
write.table(cbind(et$table,cps), paste0("edgeR.",sub(".count.txt","",info[1]),".result"), quote = FALSE, sep="\t") 

