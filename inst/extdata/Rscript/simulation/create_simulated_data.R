library(methQTL)
library(rmutil)

#'Please download the common SNPs common_all_20180423_hg19.vcf.gz database from dbSNP/UCSC at ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/
snp.all <- read.table("common_all_20180423_hg19.vcf.gz")

qtl.setOption(compute.cor.blocks = FALSE)

anno.genome <- rnb.get.annotation("probesEPIC","hg19")

parallel.setup(3)
res.all <- foreach(i=1:100,.combine = "rbind") %dopar%{
  set.seed(Sys.time())
  sel.chr <- sample(1:23,1)
  anno.chr <- anno.genome[[sel.chr]]
  sel.cpgs <- sample(1:(length(anno.chr)-1000),1)
  sel.cpgs <- anno.chr[(sel.cpgs:(sel.cpgs+999))]
  anno.cpg <- data.frame(Chromosome=seqnames(sel.cpgs),Start=start(sel.cpgs),End=end(sel.cpgs))

  snp.annot <- snp.all
  names.snps <- as.character(snp.annot$V3)
  snp.annot <- makeGRangesFromDataFrame(snp.annot,seqnames.field="V1",start.field="V2",end.field="V2")
  names(snp.annot) <- names.snps
  sel.chr <- ifelse(sel.chr==23,"X",ifelse(sel.chr==24,"Y",sel.chr))
  snp.annot <- snp.annot[seqnames(snp.annot)%in%sel.chr]

  seqlevelsStyle(snp.annot) <- "UCSC"
  keep.snps <- rep(FALSE,length(snp.annot))
  for(i in 1:length(sel.cpgs)){
   distances <- distance(sel.cpgs[i],snp.annot)
   keep.snps[distances<500000] <- TRUE
  }
  sel.snps <- snp.annot[keep.snps]
  sel.snps <- sel.snps[sample(1:length(sel.snps),2000)]

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

  snp.matrix <- matrix(0,nrow = length(sel.snps),ncol = 100)
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
    count.snp <- plyr::count(as.numeric(snp))
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
  #hist(mafs)
  for(qtl in 1:num.methQTL){
    sel.snp <- sample(which(!is.methQTL.SNP[mafs>0.2]),1)
    snp.location <- sel.snps[sel.snp]
    potential.cpgs <- distance(snp.location,sel.cpgs)<500000
    sel.cpg <- sample(which(potential.cpgs),1)
    map.cpg.snp[[sel.cpg]] <- sel.snp
    is.methQTL.CpG[sel.cpg] <- TRUE
    is.methQTL.SNP[sel.snp] <- TRUE
    vals.snp <- snp.matrix[sel.snp,]
    #' include cutoff for minimum number of heterozyotes required
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

  sel.cpgs.frame <- data.frame(Chromosome=seqnames(sel.cpgs),Start=start(sel.cpgs),End=end(sel.cpgs))
  row.names(sel.cpgs.frame) <- names(sel.cpgs)
  row.names(meth.matrix) <- names(sel.cpgs)
  sel.snps.frame <- data.frame(Chromosome=seqnames(sel.snps),Start=start(sel.snps),End=end(sel.snps))
  row.names(sel.snps.frame) <- names(sel.snps)
  row.names(snp.matrix) <- names(sel.snps)
  meth.qtl <- new("MethQTLInput",
                  meth.data=meth.matrix,
                  geno.data=snp.matrix,
                  pheno.data=data.frame(SampleID=paste0("Sample",1:100)),
                  anno.meth=sel.cpgs.frame,
                  anno.geno=sel.snps.frame,
                  samples=paste0("Sample",1:100),
                  platform="probesEPIC")
  meth.qtl.res <- do.methQTL(meth.qtl)
  cor.blocks <- getCorrelationBlocks(meth.qtl.res)[[1]]
  meth.qtl.res <- filter.pval(meth.qtl.res)
  res <- getResult(meth.qtl.res)
  true.meth.qtl.cpgs <- row.names(meth.matrix)[is.methQTL.CpG]
  true.meth.qtl.snps <- row.names(snp.matrix)[is.methQTL.SNP]
  false.meth.qtl.cpgs <- row.names(meth.matrix)[!is.methQTL.CpG]
  false.meth.qtl.snps <- row.names(snp.matrix)[!is.methQTL.SNP]
  pred.methQTL.cpgs <- unlist(lapply(res$CpG,function(cpg){
    unlist(lapply(cor.blocks,function(block){
      if(cpg %in% row.names(meth.matrix)[block]) row.names(meth.matrix)[block]
    }))
  }))
  pred.methQTL.snps <- res$SNP
  tp.cpg <- sum(pred.methQTL.cpgs %in% true.meth.qtl.cpgs)
  fp.cpg <- sum(!(pred.methQTL.cpgs %in% true.meth.qtl.cpgs))
  fn.cpg <- sum(!(true.meth.qtl.cpgs %in% pred.methQTL.cpgs))
  tn.cpg <- sum(!(false.meth.qtl.cpgs %in% pred.methQTL.cpgs))
  sens.cpg <- tp.cpg/(tp.cpg+fn.cpg)
  spec.cpg <- tn.cpg/(fp.cpg+tn.cpg)

  tp.snp <- sum(pred.methQTL.snps %in% true.meth.qtl.snps)
  fp.snp <- sum(!(pred.methQTL.snps %in% true.meth.qtl.snps))
  fn.snp <- sum(!(true.meth.qtl.snps %in% pred.methQTL.snps))
  tn.snp <- sum(!(false.meth.qtl.snps %in% pred.methQTL.snps))
  sens.snp <- tp.snp/(tp.snp+fn.snp)
  spec.snp <- tn.snp/(fp.snp+tn.snp)
  return(c(sens.cpg,spec.cpg,sens.snp,spec.snp))
}

colnames(res.all) <- c("Sensitivity.CpGs","Specificity.CpGs","Sensitivity.SNPs","Specificity.SNPs")
stats.all <- apply(res.all,2,function(x){c(mean(x),sd(x)/sqrt(length(x)))})
stats.all <- as.data.frame(stats.all)
row.names(stats.all) <- c("Mean","StandardError")
to.plot <- as.data.frame(t(stats.all))
to.plot$Stat <- rep(c("Sensitivity","Specificity"),2)
to.plot$Type <- c(rep("CpG",2),rep("SNP",2))
plot <- ggplot(to.plot,aes(x=Type,y=Stat,fill=Mean))+geom_tile()+
  geom_text(aes(label=paste0(round(Mean,3),"\u00B1",round(StandardError,3))))+
  theme_bw()+theme(text=element_text(size=15,color="black"))
ggsave("sensiticity_specificity.pdf")
# to.plot <- melt(snp.matrix)
# hist(to.plot$value)

# mafs <- c()
# for(i in 1:ncol(snp.matrix)){
#   snp <- snp.matrix[i,]
#   count.snp <- plyr::count(as.numeric(snp))
#   count.snp <- count.snp[count.snp$x>0,]
#   if(1 %in% count.snp$x){
#     homo <- (count.snp[count.snp$x==1,"freq"]*2)
#   }else{
#     homo <- 0
#   }
#   if(2 %in% count.snp$x){
#     hetero <- (count.snp[count.snp$x==2,"freq"])
#   }else{
#     hetero <- 0
#   }
#   maf <- (homo+
#             hetero)/
#     (length(snp)*2)
#   mafs <- c(mafs,maf)
# }
# hist(mafs)
