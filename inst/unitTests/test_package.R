test.meth.qtl.functions <- function(){
  logger.start("Testing methQTL package functions")
  set.seed(Sys.time())
  sel.chr <- paste0("chr",sample(1:23,1))

  # Generate methylation matrix
  sel.cpgs <- sort(sample(1:1e6,1000))
  anno.cpg <- data.frame(Chromosome=rep(sel.chr,1000),Start=sel.cpgs,End=sel.cpgs+1)

  if(requireNamespace("rmutil")){
    meth.values.cpgs <- rmutil::rbetabinom(1000,100,0.4,0.1)
  }
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

  # Generate genotype matrix
  sel.snps <- sort(sample(1:1e6,1000))
  anno.snps <- data.frame(Chromosome=rep(sel.chr,1000),Start=sel.snps,End=sel.snps)
  snp.matrix <- matrix(0,nrow = nrow(anno.snps),ncol = 100)
  for(i in 1:length(sel.snps)){
    maf <- rnbinom(1,5,0.4)/20
    if(maf>0){
      if(maf>1) maf <- 1
      num.homo <- rbinom(1,100,maf*maf)
      if(num.homo>100) num.homo <- 100
      snp.matrix[i,sample(1:100,num.homo)] <- 2
      num.hetero <- round(maf*100)
      if(num.hetero>0){
        snp.matrix[i,sample(1:100,num.hetero)] <- 1
      }
    }
  }

  num.methQTL <- 100
  is.methQTL.CpG <- rep(FALSE,1000)
  map.cpg.snp <- list()
  is.methQTL.SNP <- rep(FALSE,nrow(snp.matrix))
  mafs <- c()
  for(i in 1:ncol(snp.matrix)){
    snp <- snp.matrix[i,]
    count.snp <- count(as.numeric(snp))
    count.snp <- count.snp[count.snp$x>0,]
    if(1 %in% count.snp$x){
      homo <- (count.snp[count.snp$x==1,"freq"]*2)
    }else{
      homo <- 0
    }
    if(2 %in% count.snp$x){
      hetero <- (count.snp[count.snp$x==2,"freq"])
    }else{
      hetero <- 0
    }
    maf <- (homo+
              hetero)/
      (length(snp)*2)
    mafs <- c(mafs,maf)
  }

  # Generate methQTL interactions
  for(qtl in 1:num.methQTL){
    sel.snp <- sample(which(!is.methQTL.SNP[mafs>0.2]),1)
    snp.location <- anno.snps[sel.snp,]
    potential.cpgs <- abs(snp.location$Start-anno.cpg$Start)<500000
    sel.cpg <- sample(which(potential.cpgs),1)
    map.cpg.snp[[sel.cpg]] <- sel.snp
    is.methQTL.CpG[sel.cpg] <- TRUE
    is.methQTL.SNP[sel.snp] <- TRUE
    vals.snp <- snp.matrix[sel.snp,]
    # include cutoff for minimum number of heterozyotes required
    vals.cpg <- meth.matrix[sel.cpg,]
    effect.size <- rnorm(1,mean=0.2,sd=0.05)
    if(mean(vals.cpg)>0.95){
      effect.size <- -effect.size
    }else  if(sample(c(0,1),1)>0 & mean(vals.cpg)>0.5){
      effect.size <- -effect.size
    }
    new.values <- vals.cpg + vals.snp*effect.size + rnorm(1,sd=0.05)
    new.values[new.values>1] <- 1
    new.values[new.values<0] <- 0
    meth.matrix[sel.cpg,] <- new.values
  }

  row.names(meth.matrix) <- paste0("cg",1:nrow(meth.matrix))
  row.names(snp.matrix) <- paste0("rs",1:nrow(snp.matrix))
  row.names(anno.cpg) <- paste0("cg",1:nrow(meth.matrix))
  row.names(anno.snps) <- paste0("rs",1:nrow(snp.matrix))
  qtl.setOption(meth.qtl.type="allVSall")
  meth.qtl <- new("methQTLInput",
                  meth.data=meth.matrix,
                  geno.data=snp.matrix,
                  pheno.data=data.frame(SampleID=paste0("Sample",1:100)),
                  anno.meth=anno.cpg,
                  anno.geno=anno.snps,
                  samples=paste0("Sample",1:100),
                  platform="probesEPIC")
  passed <- inherits(meth.qtl,"methQTLInput")
  meth.qtl.res <- doMethQTLChromosome(meth.qtl,chrom = sel.chr,sel.covariates = NULL,p.val.cutoff = 1e-2)
  true.meth.qtl.cpgs <- row.names(meth.matrix)[is.methQTL.CpG]
  true.meth.qtl.snps <- row.names(snp.matrix)[is.methQTL.SNP]
  cor.blocks <- getCorrelationBlocks(meth.qtl.res)
  res <- getResult(meth.qtl.res,cor.blocks)
  pred.methQTL.cpgs <- unlist(res$CorrelationBlock)
  pred.methQTL.snps <- res$SNP
  tp.cpg <- sum(pred.methQTL.cpgs %in% true.meth.qtl.cpgs)
  tp.snp <- sum(pred.methQTL.snps %in% true.meth.qtl.snps)

  passed <- passed & inherits(meth.qtl.res,"methQTLResult") & (tp.cpg>0) & (tp.snp>0)
  checkTrue(passed)
  logger.completed()
}
