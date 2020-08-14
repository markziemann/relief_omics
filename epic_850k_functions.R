
# Here is a function to make a volcano plot
make_volcano <- function(dm,name,mx) {
    sig <- subset(dm,adj.P.Val<0.05)
    N_SIG=nrow(sig)
    N_UP=nrow(subset(sig,logFC>0))
    N_DN=nrow(subset(sig,logFC<0))
    HEADER=paste(N_SIG,"@5%FDR,", N_UP, "up", N_DN, "dn")
    plot(dm$logFC,-log10(dm$P.Val),cex=0.5,pch=19,col="darkgray",
        main=name, xlab="log FC", ylab="-log10 pval")
    mtext(HEADER)
    grid()
    points(sig$logFC,-log10(sig$P.Val),cex=0.5,pch=19,col="red")
}
# Here is a function to make heatmaps 
make_heatmap <- function(dm,name,mx,n, groups) {
  topgenes <-  rownames(head(dm[order(dm$P.Value),],n))
  ss <- mx[which(rownames(mx) %in% topgenes),]
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
  colCols <- as.character(groups)
  colCols <- gsub("1","orange",colCols)
  colCols <- gsub("0","yellow",colCols)
  heatmap.2(ss,scale="row",margin=c(10, 10),cexRow=0.4,trace="none",cexCol=0.4,
  ColSideColors=colCols ,  col=my_palette, main=name)
}
# make beeswarm charts
# dm = a limma differential meth object
# name = character name of the limma dm object
# mx = matrix of normalised data
# groups = a vector of factors corresponding to the cols in mx
# n = the number of top significant genes to plot (default = 15) 
make_beeswarms <- function(dm,name,beta,groups,n=15) {
    par(mar=c(3,3,1,1))
    NCOLS=5
    NROWS=floor(n/NCOLS)
    if (n %% NCOLS > 0) { NROWS <- NROWS + 1 }
    par(mfrow=c(NROWS, NCOLS))
    topgenes <-  rownames(head(dm[order(dm$P.Value),],n))
    ss <- beta[which(rownames(beta) %in% topgenes),]
    n <- 1:n
    g1name=levels(groups)[1]
    g2name=levels(groups)[2]
    g1dat <- ss[n,which(groups == g1name)]
    g2dat <- ss[n,which(groups == g2name)]
    g1l <-lapply(split(g1dat, row.names(g1dat)), unlist)
    g2l <-lapply(split(g2dat, row.names(g2dat)), unlist)
    for (i in n) {
      mydat <- list(g1l[[i]],g2l[[i]])
        beeswarm(mydat,ylim=c(0,1),cex=0.2, pch=19,
        las=2, cex.lab=0.6, main=names( g1l )[i] , 
        ylab="",labels = c(g1name,g2name))
      grid()
    }
}
# make beeswarm charts for best confects
# dm = a limma differential meth object
# name = character name of the limma dm object
# mx = matrix of normalised data
# groups = a vector of factors corresponding to the cols in mx
# n = the number of top significant genes to plot (default = 15) 
make_beeswarms_confects <- function(confects,name,beta,groups,n=15) {
    par(mar=c(3,3,1,1))
    NCOLS=5
    NROWS=floor(n/NCOLS)
    if (n %% NCOLS > 0) { NROWS <- NROWS + 1 }
    par(mfrow=c(NROWS, NCOLS))
    topgenes <-  head(confects$table,n)$name
    ss <- beta[which(rownames(beta) %in% topgenes),]
    n <- 1:n
    g1name=levels(groups)[1]
    g2name=levels(groups)[2]
    g1dat <- ss[n,which(groups == g1name)]
    g2dat <- ss[n,which(groups == g2name)]
    g1l <-lapply(split(g1dat, row.names(g1dat)), unlist)
    g2l <-lapply(split(g2dat, row.names(g2dat)), unlist)
    for (i in n) {
      mydat <- list(g1l[[i]],g2l[[i]])
        beeswarm(mydat,ylim=c(0,1),cex=0.2, pch=19,
        las=2, cex.lab=0.6, main=names( g1l )[i] , 
        ylab="",labels = c(g1name,g2name))
      grid()
    }
}
# this is a wrapper which creates three charts
# We will be adding more
make_dm_plots <- function(dm,name,mx,groups,confects,beta) {
    make_volcano(dm=dm,name=name,mx=mx)
    make_beeswarms(dm=dm ,name=name , beta=beta , groups=groups , n= 15)
    make_heatmap(dm=dm , name=name , mx=mx ,n = 50, groups=groups)
    make_beeswarms_confects(confects=confects, name=name, beta=beta, groups=groups, n=15)
}  
