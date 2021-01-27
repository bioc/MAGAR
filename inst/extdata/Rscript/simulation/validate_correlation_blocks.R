library(methQTL)
library(rmutil)
library(foreach)
library(doParallel)

anno.genome <- rnb.get.annotation("probesEPIC","hg19")

#' Plot cluster.cor.threshold against the number of clusters
params <- seq(0,1,by=0.05)
cluster.sizes <- list()
qtlSetOption(standard.deviation.gauss = 3000)
parallel.setup(10)
cluster.sizes <- foreach(param=1:length(params),.combine="cbind")%dopar%{
  param <- params[param]
  qtlSetOption(cluster.cor.threshold = param)
  cluster.sizes.param <- unlist(lapply(1:100,function(foo){
    sel.chr <- sample(1:23,1)
    anno.chr <- anno.genome[[sel.chr]]
    sel.cpgs <- sample(1:(length(anno.chr)-1000),1)
    sel.cpgs <- anno.chr[(sel.cpgs:(sel.cpgs+999))]
    anno.cpg <- data.frame(Chromosome=seqnames(sel.cpgs),Start=start(sel.cpgs),End=end(sel.cpgs))

    meth.values.cpgs <- rbetabinom(1000,100,0.4,0.1)
    num.clusters <- sample(200:500,1)
    is.changed <- rep(FALSE,length(100))
    for(i in 1:num.clusters){
      if(all(is.changed)) break
      clust.size <- sample(0:9,1)
      clust.location <- sample(which(!is.changed),1)
      meth.values.cpgs[clust.location:(clust.location+clust.size)] <- meth.values.cpgs[clust.location]
      is.changed[clust.location:(clust.location+clust.size)] <- TRUE
    }
    meth.matrix <- matrix(rep(meth.values.cpgs,100),nrow = 1000,ncol = 100)/100
    meth.matrix <- apply(meth.matrix,c(1,2),function(cpg){
      cpg <- cpg+rnorm(1,sd=0.05)
      ifelse(cpg>1,1,ifelse(cpg<0,0,cpg))
    })
    cor.blocks <- computeCorrelationBlocks(meth.matrix,anno.cpg)
    num.clusters-length(cor.blocks)
  }))
  cluster.sizes.param
}
mean.sizes <- apply(cluster.sizes,2,mean)
se.sizes <- apply(cluster.sizes,2,function(x){
  sd(x)/sqrt(length(x))
})
to.plot <- data.frame(params,mean.sizes,se.sizes)
plot <- ggplot(to.plot,aes(x=params,y=mean.sizes,ymin=mean.sizes-2*se.sizes,ymax=mean.sizes+2*se.sizes))+
  geom_point()+geom_line()+geom_errorbar()+theme_bw()+geom_vline(xintercept = 0.25)+
  ylab("Number of expected clusters - number of clusters")+xlab("Correlation threshold")+theme(text=element_text(size=15))


#' Parameter standard deviation Gauss
params <- seq(4000,5000,by=100)
qtl.setOption(cluster.cor.threshold = 0.2)
parallel.setup(10)
cluster.sizes <- foreach(param=params,.combine="cbind")%dopar%{
  qtl.setOption(standard.deviation.gauss = param)
  cluster.sizes.param <- unlist(lapply(1:100,function(foo){
    sel.chr <- sample(1:23,1)
    anno.chr <- anno.genome[[sel.chr]]
    sel.cpgs <- sample(1:(length(anno.chr)-1000),1)
    sel.cpgs <- anno.chr[(sel.cpgs:(sel.cpgs+999))]
    anno.cpg <- data.frame(Chromosome=seqnames(sel.cpgs),Start=start(sel.cpgs),End=end(sel.cpgs))

    meth.values.cpgs <- rbetabinom(1000,100,0.4,0.1)
    num.clusters <- sample(200:500,1)
    is.changed <- rep(FALSE,length(100))
    for(i in 1:num.clusters){
      if(all(is.changed)) break
      clust.size <- sample(0:9,1)
      clust.location <- sample(which(!is.changed),1)
      meth.values.cpgs[clust.location:(clust.location+clust.size)] <- meth.values.cpgs[clust.location]
      is.changed[clust.location:(clust.location+clust.size)] <- TRUE
    }
    meth.matrix <- matrix(rep(meth.values.cpgs,100),nrow = 1000,ncol = 100)/100
    meth.matrix <- apply(meth.matrix,c(1,2),function(cpg){
      cpg <- cpg+rnorm(1,sd=0.05)
      ifelse(cpg>1,1,ifelse(cpg<0,0,cpg))
    })
    cor.blocks <- compute.correlation.blocks(meth.matrix,anno.cpg)
    num.clusters-length(cor.blocks)
  }))
  cluster.sizes.param
}
mean.sizes <- apply(cluster.sizes,2,mean)
se.sizes <- apply(cluster.sizes,2,function(x){
  sd(x)/sqrt(length(x))
})
to.plot <- data.frame(params,mean.sizes,se.sizes)
plot <- ggplot(to.plot,aes(x=params,y=mean.sizes,ymin=mean.sizes-2*se.sizes,ymax=mean.sizes+2*se.sizes))+
  geom_point()+geom_line()+geom_errorbar()+theme_bw()+
  ylab("Number of expected clusters - number of clusters")+xlab("Standard deviation Gauss")+theme(text=element_text(size=15))

#' Parameter absolute distance cutoff
params <- seq(100000,1000000,by=100000)
cluster.sizes <- list()
qtl.setOption(cluster.cor.threshold = 0.2)
qtl.setOption(standard.deviation.gauss = 3000)
cluster.sizes <- foreach(param=params,.combine="cbind")%dopar%{
  param <- params[param]
  qtl.setOption(absolute.distance.cutoff = param)
  cluster.sizes.param <- unlist(lapply(1:1000,function(foo){
    sel.chr <- sample(1:23,1)
    anno.chr <- anno.genome[[sel.chr]]
    sel.cpgs <- sample(1:(length(anno.chr)-1000),1)
    sel.cpgs <- anno.chr[(sel.cpgs:(sel.cpgs+999))]
    anno.cpg <- data.frame(Chromosome=seqnames(sel.cpgs),Start=start(sel.cpgs),End=end(sel.cpgs))

    meth.values.cpgs <- rbetabinom(1000,100,0.4,0.1)
    num.clusters <- sample(200:500,1)
    is.changed <- rep(FALSE,length(100))
    for(i in 1:num.clusters){
      if(all(is.changed)) break
      clust.size <- sample(0:9,1)
      clust.location <- sample(which(!is.changed),1)
      meth.values.cpgs[clust.location:(clust.location+clust.size)] <- meth.values.cpgs[clust.location]
      is.changed[clust.location:(clust.location+clust.size)] <- TRUE
    }
    meth.matrix <- matrix(rep(meth.values.cpgs,100),nrow = 1000,ncol = 100)/100
    meth.matrix <- apply(meth.matrix,c(1,2),function(cpg){
      cpg <- cpg+rnorm(1,sd=0.05)
      ifelse(cpg>1,1,ifelse(cpg<0,0,cpg))
    })
    cor.blocks <- compute.correlation.blocks(meth.matrix,anno.cpg)
    num.clusters-length(cor.blocks)
  }))
  cluster.sizes.param
}
mean.sizes <- apply(cluster.sizes,2,mean)
se.sizes <- apply(cluster.sizes,2,function(x){
  sd(x)/sqrt(length(x))
})
to.plot <- data.frame(params,mean.sizes,se.sizes)
plot <- ggplot(to.plot,aes(x=params,y=mean.sizes,ymin=mean.sizes-2*se.sizes,ymax=mean.sizes+2*se.sizes))+
  geom_point()+geom_line()+geom_errorbar()+theme_bw()+
  ylab("Number of expected clusters - number of clusters")+xlab("Distance cutoff")+theme(text=element_text(size=15))

