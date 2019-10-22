##########################################################################################
# compute_methQTL.R
# created: 2019-08-29
# creator: Michael Scherer
# ---------------------------------------------------------------------------------------
# Methods to call methQTL from DNA methylation and genotyping data.
##########################################################################################

#' do.methQTL
#'
#' Function to compute methQTL given DNA methylation and genotyping data.
#'
#' @param meth.qtl An object of type \code{\link{methQTLInput-class}} on which methQTL computation is to be performed
#' @param sel.covariates Covariates as column names of the sample annotation sheet stored in \code{meth.qtl} to be
#'          used for covariate adjustment.
#' @param p.val.cutoff The p-value cutoff used for methQTL calling
#' @param ncores The number of cores used.
#' @param cluster.submit Flag indicating if jobs are to be distributed among a SGE compute cluster
#' @param out.dir Output directory
#' @return An object of type \code{\link{methQTLResult-class}} containing the called methQTL interactions.
#' @details The process is split into 4 steps:
#'          \describe{
#'            \item{1}{First the two matrices are split according to the chromosomes.}
#'            \item{2}{We then compute correlations among the CpGs and compute CpG correlation blocks.}
#'            \item{3}{In each of the CpG correlation blocks, linear models according to the \code{linear.model.type}
#'              \code{\link{qtl.setOption}} with the CpG methylation state of the reference CpG specified by
#'              \code{representative.cpg.computation} as output and the SNP genotype state and all possible covariates
#'              as input are computed.}
#'            \item{4}{For each of the CpG correlation blocks, we report the p-value of the representative CpG.}
#'          }
#' @author Michael Scherer
#' @export
do.methQTL <- function(meth.qtl,sel.covariates=NULL,p.val.cutoff=1e-5,ncores=1,cluster.submit=F,out.dir=getwd()){
  if(!inherits(meth.qtl,"methQTLInput")){
    stop("Invalid value for meth.qtl, needs to be of type methQTLInput")
  }
  if(!meth.qtl@imputed){
    meth.qtl <- impute.meth(meth.qtl)
  }
  all.chroms <- unique(getAnno(meth.qtl)$Chromosome)
  if(ncores>1){
    parallel.setup(ncores)
  }
  logger.start("Computing methQTLs")
  # res.all <- foreach(chrom=all.chroms,.combine="cbind",.export=c()) %dopar%{
  #   res.chrom <- do.methQTL.chromosome(meth.qtl,chrom,sel.covariates,p.val.cutoff)
  # }
  res.all <- list()
  if(!cluster.submit){
    for(chrom in all.chroms){
      res.chrom <- do.methQTL.chromosome(meth.qtl,chrom,sel.covariates,p.val.cutoff)
      res.all[[chrom]] <- res.chrom
    }
    res.all <- join.methQTLResult(res.all)
  }else{
    res.all <- submit.cluster.jobs(meth.qtl,sel.covariates,p.val.cutoff,out.dir)
  }
  logger.completed()
  return(res.all)
}

#' do.methQTL.chromosome
#'
#' This functions computes the methQTL interactions for a single chromosome
#'
#' @param meth.qtl An Object of type \code{\link{methQTLInput-class}}.
#' @param chrom Character vector represeting the chromosome to be investigated.
#' @param sel.covariates Covariates as column names of the sample annotation sheet stored in \code{meth.qtl} to be
#'          used for covariate adjustment.
#' @param p.val.cutoff The p-value used for methQTL calling
#' @return A data frame with seven columns:
#'        \describe{
#'          \item{CpGs}{The CpG ID chosen to be the representative CpG in the methQTL}
#'          \item{SNP}{The SNP ID (as rsNNNNNN) involved in the methQTL}
#'          \item{Beta}{The coefficient estimate of the linear model}
#'          \item{P.value}{The p-value associated with the coefficient estimate}
#'          \item{Chromosome}{The chromosome name}
#'          \item{Position.CpG}{The genomic position of the CpG}
#'          \item{Position.SNP}{The genomic position of the SNP}
#'          \item{Distance}{The distance between the CpG and the SNP}
#'        }
#' @author Michael Scherer
#' @export
do.methQTL.chromosome <- function(meth.qtl,chrom,sel.covariates,p.val.cutoff){
  logger.start(paste("Computing methQTL for chromosome",chrom))
  anno <- getAnno(meth.qtl,"meth")
  sel.meth <- which(anno$Chromosome %in% chrom)
  sel.anno <- anno[sel.meth,]
  sel.meth <- getMethData(meth.qtl)[sel.meth,]
  cor.blocks <- compute.correlation.blocks(sel.meth,sel.anno)
  cor.blocks <- lapply(cor.blocks,as.numeric)
  anno.geno <- getAnno(meth.qtl,"geno")
  sel.geno <- which(anno.geno$Chromosome %in% chrom)
  sel.anno.geno <- anno.geno[sel.geno,]
  sel.geno <- getGeno(meth.qtl)[sel.geno,]
  ph.dat <- getPheno(meth.qtl)
  if(is.null(sel.covariates)){
    ph.dat <- NULL
  }else{
    if(!all(sel.covariates %in% colnames(ph.dat))){
      logger.warning("Not all the covariates are present in the sample annotation sheet.")
      sel.covariates <- sel.covariates[sel.covariates %in% colnames(ph.dat)]
      if(is.null(sel.covariates)){
        stop("No valid covariate present")
      }
    }
    ph.dat <- ph.dat[,sel.covariates,drop=FALSE]
  }
  logger.start("Compute methQTL per correlation block")
  res.chr.p.val <- t(sapply(cor.blocks,call.methQTL.block,sel.meth,sel.geno,ph.dat,sel.anno,sel.anno.geno))
  for(j in 1:ncol(res.chr.p.val)){
    if(j %in% c(1,2,5)){
      res.chr.p.val[,j] <- unlist(lapply(res.chr.p.val[,j],as.character))
    }else{
      res.chr.p.val[,j] <- unlist(lapply(res.chr.p.val[,j],as.numeric))
    }
  }
  res.chr.p.val <- as.data.frame(apply(res.chr.p.val,2,unlist))
  logger.completed()
  res.chr.p.val <- res.chr.p.val[as.numeric(as.character(res.chr.p.val$P.value))<p.val.cutoff,]
  if(nrow(res.chr.p.val)==0){
    logger.info(paste("No methQTL found for chromosome",chrom))
    chrom.frame <- data.frame()
  }else{
    tests.performed <- length(cor.blocks)*nrow(sel.geno)
    chrom.frame <- data.frame(CpG=as.character(res.chr.p.val$CpG),
                            SNP=as.character(res.chr.p.val$SNP),
                            Beta=as.numeric(as.character(res.chr.p.val$Beta)),
                            P.value=as.numeric(as.character(res.chr.p.val$P.value)),
                            Chromosome=as.character(res.chr.p.val$Chromosome),
                            Position.CpG=as.numeric(as.character(res.chr.p.val$Position_CpG)),
                            Position.SNP=as.numeric(as.character(res.chr.p.val$Position_SNP)))
    chrom.frame$Distance <- chrom.frame$Position.CpG - chrom.frame$Position.SNP
    chrom.frame$p.val.adj.fdr <- p.adjust(chrom.frame$P.value,method="fdr",n=tests.performed)
  }
  methQTL.result <- new("methQTLResult",
                        result.frame=chrom.frame,
                        anno.meth=sel.anno,
                        anno.geno=sel.anno.geno,
                        method=qtl.getOption("linear.model.type"),
                        rep.type=qtl.getOption("representative.cpg.computation"))
  logger.completed()
  return(methQTL.result)
}

#' compute.correlation.blocks
#'
#' This function computes CpG correlation blocks from correlations of CpGs across samples by Louvian
#' clustering.
#'
#' @param meth.data A \code{data.frame} containing the methylation data with CpGs in the rows and samples in the columns.
#' @param annotation The genomic annotations of the CpG positions.
#' @param cor.threshold The correlation threshold used to discard edges from the correlation-based network.
#' @param sd.gauss Standard deviation of the Gauss distribution used to weight the distance
#' @param absolute.cutoff Absolute distance cutoff after which no methQTL interaction is to be considered.
#' @param max.cpgs Maximum number of CpGs used in the computation (used to save memory). 40,000 is a reasonable
#'             default for machines with ~128GB of main memory. Should be smaller for smaller machines and larger
#'             for larger ones.
#' @return A list representing the clustering of CpGs into correlation blocks. Each element is a cluster, which contains
#'      row indices of the DNA methylation matrix that correspond to this cluster.
#' @details This method performs clustering of the correlation matrix obtaind from the DNA methylation matrix. Correlations
#'      are computed for each pair of CpGs across all the samples. We then compute a similarity matrix from this correlation
#'      matrix and set correlations lower than the given threshold to 0. In the next step, we weight the correlations
#'      by the distance between the CpGs: smaller distances get higher weights according to Gaussian distribution with
#'      mean 0 and standard deviation as specified above. Furthermore, similarities of CpGs that are further than
#'      \code{absolute.distance.cutoff} away from one another are discarded.
#'
#'      We then compute the associated weighted, undirected graph from the similarity matrix and execute Louvain clustering
#'      on the graph. The resulting clusters of CpGs are returned.
#'
#' @author Michael Scherer
#' @export
compute.correlation.blocks <- function(meth.data,annotation,cor.threshold=qtl.getOption("cluster.cor.threshold"),
                                       sd.gauss=qtl.getOption("standard.deviation.gauss"),
                                       absolute.cutoff=qtl.getOption("absolute.distance.cutoff"),
                                       max.cpgs=qtl.getOption("max.cpgs")){
  logger.start("Compute correlation blocks")
  if(nrow(annotation)>max.cpgs){
    logger.info(paste("Split workload, since facing",nrow(annotation),"CpGs (Maximum is",max.cpgs,")"))
    bin.split <- round(nrow(annotation)/2)
    return(c(compute.correlation.blocks(meth.data=meth.data[1:bin.split,],
                                        annotation=annotation[1:bin.split,],
                                        cor.threshold = cor.threshold,
                                        sd.gauss = sd.gauss,
                                        absolute.cutoff = absolute.cutoff,
                                        max.cpgs = max.cpgs),
             compute.correlation.blocks(meth.data=meth.data[(bin.split+1):nrow(annotation),],
                                        annotation=annotation[(bin.split+1):nrow(annotation),],
                                        cor.threshold = cor.threshold,
                                        sd.gauss = sd.gauss,
                                        absolute.cutoff = absolute.cutoff,
                                        max.cpgs = max.cpgs)
             ))
  }
  logger.start("Compute correlation matrix")
  cor.all <- big_cor(as_FBM(t(as.matrix(meth.data)),type="double"))
  rm(meth.data)
  logger.completed()
  cor.all <- cor.all[,,drop=F]
  if(qtl.getOption("hdf5dump")){
    cor.all <- writeHDF5Array(cor.all)
  }
  logger.start("Compute correlation distance")
  dist.all <- cor2dist(cor.all)
  if(qtl.getOption("hdf5dump")){
    dist.all <- writeHDF5Array(dist.all)
  }
  logger.completed()
  rep.vals <- cor.all<cor.threshold
  if(qtl.getOption("hdf5dump")){
    rep.vals <- writeHDF5Array(rep.vals)
  }
  dist.all[rep.vals] <- 0
  genomic.positions <- annotation$Start
  logger.start("Compute pairwise distances")
  gc()
  pairwise.distance <- abs(as.data.frame(lapply(genomic.positions,function(x)x-genomic.positions)))
  logger.completed()
  rep.vals <- pairwise.distance>absolute.cutoff
  if(qtl.getOption("hdf5dump")){
    rep.vals <- writeHDF5Array(rep.vals)
  }
  dist.all[rep.vals] <- 0
  gc()
  logger.start("Weight distances")
  if(qtl.getOption("hdf5dump")){
    weighted.distances <- matrix(nrow=nrow(dist.all),ncol=ncol(dist.all))
    weighted.distances <- writeHDF5Array(weighted.distances)
    chunk.size <- 10000
    i <- 1
    while(i < nrow(dist.all)){
      if((i + chunk.size)>nrow(dist.all)){
        do.work <- i:nrow(dist.all)
        weighted.distances[do.work,] <- as.matrix(dist.all[do.work,])*dnorm(as.matrix(pairwise.distance[do.work,]),0,sd.gauss)
        break
      }
      do.work <- i:(i+chunk.size)
      weighted.distances[do.work,] <- as.matrix(dist.all[do.work,])*dnorm(as.matrix(pairwise.distance[do.work,]),0,sd.gauss)
      i <- i+chunk.size+1
    }
  }else{
    weighted.distances <- dist.all*dnorm(as.matrix(pairwise.distance),0,sd.gauss)
  }
  logger.completed()
  colnames(weighted.distances) <- as.character(1:ncol(weighted.distances))
  rownames(weighted.distances) <- as.character(1:nrow(weighted.distances))
  rm(dist.all)
  rm(rep.vals)
  rm(cor.all)
  gc()
  logger.start("Compute graph")
  graph.ad <- graph.adjacency(as.matrix(weighted.distances),weighted = T,mode = "undirected",diag=F)
  logger.completed()
  logger.start("Compute clustering")
  clust <- cluster_louvain(graph.ad)
  logger.completed()
  logger.completed()
  return(groups(clust))
}

#' call.methQTL.block
#'
#' NOTE by TL: Choosing the best CpG might lead to detecting only outliers. It might thus be better to choose a
#' representative by correlation block.
#'
#' This function computes a methQTL per correlation block.
#' @param cor.block A single correlation block as determined by \code{\link{compute.correlation.blocks}}
#' @param meth.data The methylation data matrix
#' @param geno.data The genotype data matrix
#' @param covs A \code{data.frame} specifying the covariates used for the analysis
#' @param anno.meth A \code{data.frame} containing genomic annotations of the CpGs
#' @param anno.geno A \code{data.frame} containing genomic annotations for the SNPs
#' @param model.type Type of the linear model to be used. Supported options are:
#'    \describe{
#'      \item{categorical.anova}{Performs linear regression after transforming the genotype states (0=homozygous
#'      reference allele, 1=heteroygous, 2=homoygous alternate allele) into categorical variables and then perform
#'      ANOVA}
#'      \item{classical.linear}{Performs linear regression with keeping the genotype calls as numeric values. This
#'      induces an ordering between heterozygous and homozygous reference. We then use the p-value from the t-statistic
#'      and the associated beta value.}
#'    }
#' @param repr.type Argument specifying how reference CpGs per correlation block are to be computed. Available
#'            options are \code{"row.medians"} for the site that is the row median across the samples within the
#'            correlation block (for ties a random selection is performed), \code{"mean.center"} for an artifical site in the geometric center of the block with
#'            the average methylation level or \code{"best.all"} for the CpG with the best p-value across all of the
#'            CpGs in the correlation block.
#' @return A vector of three elements:
#'    \describe{
#'      \item{1}{The p-value of the best methQTL}
#'      \item{2}{The estimated coefficient of the linear model}
#'      \item{3}{The position of the SNP}
#'      \item{4}{The position of the CpG}
#'    }
#' @details methQTL are called using a linear modeling strategy: The DNA methylation state of the CpG is considered
#'     the output, while the SNP genotype and the potential covariates are the input. We compute a linear model
#'     for each CpG in the correlation block and each SNP. We then select the best methQTL according to the lowest
#'     p-value and return it.
#' @author Michael Scherer
#' @export
call.methQTL.block <- function(cor.block,meth.data,geno.data,covs,anno.meth,anno.geno,
                               model.type=qtl.getOption("linear.model.type"),
                               repr.type=qtl.getOption("representative.cpg.computation")){
  sel.meth <- meth.data[cor.block,,drop=F]
  anno.meth <- anno.meth[cor.block,,drop=F]
  if(repr.type == "row.medians"){
    reps <- apply(as.matrix(sel.meth),1,median,na.rm=T)
    order.reps <- order(reps,decreasing = T)
    n.cpgs <- length(order.reps)
    if(n.cpgs%%2 ==0){
      sel.site <- sample(c(n.cpgs/2,n.cpgs/2 + 1),1)
    }else{
      sel.site <- n.cpgs/2 + 0.5
    }
    reps <- as.matrix(sel.meth)[sel.site,]
    sel.anno <- anno.meth[sel.site,]
  }else if(repr.type == "mean.center"){
    reps <- apply(as.matrix(sel.meth),2,mean,na.rm=T)
    sel.anno <- data.frame(Chromosome=unique(anno.meth$Chromosome),Start=mean(anno.meth$Start))
    row.names(sel.anno) <- pasteo("mean_of_",nrow(sel.meth))
  }
  distances <- abs(anno.geno$Start - sel.anno$Start)
  if(all(distances>=qtl.getOption("absolute.distance.cutoff"))){
    logger.info(paste("No SNP closer than",qtl.getOption("absolute.distance.cutoff")))
    ret <- data.frame(as.character(row.names(sel.anno)),
                      NA,
                      NA,
                      NA,
                      as.character(sel.anno$Chromosome),
                      NA,
                      sel.anno$Start)
    colnames(ret) <- c("CpG","SNP","P.value","Beta","Chromosome","Position_SNP","Position_CpG")
    return(ret)
  }
  geno.data <- geno.data[distances<qtl.getOption("absolute.distance.cutoff"),]
  anno.geno <- anno.geno[distances<qtl.getOption("absolute.distance.cutoff"),]
  all.snps <- apply(geno.data,1,function(snp){
    if(is.null(covs)){
      form <- as.formula("CpG~SNP")
    }else{
      form <- as.formula(paste("CpG~SNP",colnames(covs),sep="+"))
    }
    if(model.type == "categorical.anova"){
      in.mat <- data.frame(CpG=reps,SNP=as.factor(snp))
      if(!is.null(covs)){
        in.mat <- data.frame(in.mat,covs)
      }
      lm.model <- lm(form,data=in.mat)
      an.model <- anova(lm.model)
      p.val <- an.model["SNP","Pr(>F)"]
      return(c(p.val=p.val,beta=NA))
    }else if(model.type == "classical.linear"){
      in.mat <- data.frame(CpG=reps,SNP=snp)
      if(!is.null(covs)){
        in.mat <- data.frame(in.mat,covs)
      }
      lm.model <- lm(form,data=in.mat)
      p.val <- summary(lm.model)$coefficients["SNP","Pr(>|t|)"]
      beta <- summary(lm.model)$coefficients["SNP","Estimate"]
      return(c(p.val=p.val,beta=beta))
    }
  })
  all.snps <- t(all.snps)
  is.min <- which.min(all.snps[,'p.val'])
  min.p.val <- data.frame(as.character(row.names(sel.anno)),
                          as.character(row.names(anno.geno)[is.min]),
                          all.snps[is.min,'p.val'],
                          all.snps[is.min,'beta'],
                          as.character(sel.anno$Chromosome),
                          anno.geno$Start[is.min],
                          sel.anno$Start)
  colnames(min.p.val) <- c("CpG","SNP","P.value","Beta","Chromosome","Position_SNP","Position_CpG")
  return(min.p.val)
}
