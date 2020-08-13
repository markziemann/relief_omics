# The purpose here is todetermine whether parts of the genome have a difference distribution 
# in gene expression.
# The purpose of this was to analyse whether genes on the Y chromosome have a lower variation

ex <- t0_crp

ex <- ex[,7:ncol(ex)]

mysd <- data.frame( apply(ex,1,function(x) {sd(x)/mean(x) } ) ) 

mychr <- read.table("../refgenome/chromosomes2.tsv",sep="\t",row.names=2)

mychr2 <- merge(mysd,mychr,by=0)

mychr2[,1]=NULL

agg <- aggregate(. ~ V1, mychr2, mean)

agg <- agg[order(agg[,2]),]

rownames(agg) <- agg[,1]

agg[,1]=NULL

colnames(agg) = "mean variance"

tmp <- sapply(unique(mychr2[,2]), function(x) { tmp <- mychr2[which(mychr2[,2] %in% x),1] } )

mymeans <- sapply(tmp,mean)

new <- tmp[order(mymeans)]

boxplot(new,horizontal=TRUE,las=2,ylim=c(0,.2),main="CV of gene expression by chromosome")
# chrY genes have an overall higher CV in contrast to what the heatma shows
