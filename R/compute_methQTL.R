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
#' @return An object of type \code{\link{methQTLResult}} containing the called methQTL interactions.
#' @details The process is split into 4 steps:
#'          \describe{
#'            \item{1}{First the two matrices are split according to the chromosomes.}
#'            \item{2}{We then compute correlations among the CpGs and compute CpG correlation blocks.}
#'            \item{3}{In each of the CpG correlation blocks, simple linear models with the CpG methylation state
#'              of all CpGs as output and the SNP genotype state and all possible covariates as input are computed.}
#'            \item{4}{For each of the CpG correlation blocks, the CpGs with the lowest p-value is selected as the
#'              representative methQTL per block.}
#'          }
#' @author Michael Scherer
#' @export
do.methQTL <- function(meth.qtl,sel.covariates=NULL,p.val.cutoff=1e-5,ncores=1){
  if(!inherits(meth.qtl,"methQTLInput")){
    stop("Invalid value for meth.qtl, needs to be of type methQTLInput")
  }
  all.chroms <- unique(getAnno(meth.qtl)$Chromosome)
  if(ncores>1){
    require(doParallel)
    parallel.setup(ncores)
  }
  res.all <- foreach(chrom=all.chroms,.combine="cbind",.export=c()) %dopar%{
    res.chrom <- do.methQTL.chromosome(meth.qtl,chrom,sel.covariates,p.val.cutoff)
  }
}

#' do.methQTL.chromosome
#'
#' This functions computes the methQTL interactions for a single chromosome
#'
#' @param meth.qtl An Object of type \code{\link{methQTLInput-class}}.
#' @param chrom Character vector represeting the chromosome to be investigated.
#' @param sel.covariates A \code{data.fr

do.methQTL.chromosome <- function(meth.qtl,chrom,sel.covariates,p.val.cutoff){
  anno <- getAnno(meth.qtl,"meth")
  sel.meth <- which(anno$Chromosome %in% chrom)
  sel.anno <- anno[sel.meth,]
  sel.meth <- getMethData(meth.qtl)[sel.meth,]
  cor.blocks <- compute.correlation.blocks(sel.meth,sel.anno)
  anno.geno <- getAnno(meth.qtl,"geno")
  sel.geno <- which(anno.geno$Chromosome %in% chrom)
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
      ph.dat <- ph.dat[,sel.covariates]
    }
  }
  res.chr.p.val <- sapply(cor.blocks,call.methQTL.block,sel.meth,sel.geno,ph.dat)
  tests.performed <- length(cor.blocks)*nrow(sel.geno)
  chrom.frame <- data.frame(CpG=row.names(anno)[res.chr.p.val$Position_CpG],
                            SNP=row.names(anno.geno)[res.chr.p.val$Position_SNP],
                            Beta=res.chr.p.val$Beta,
                            P.palue=res.chr.p.val$P.value,
                            Chromosome=anno$Chromosome[res.chr.p.val$Position_CpG],
                            Position.CpG=anno$Start[res.chr.p.val$Position_CpG],
                            Position.SNP=anno$Start[res.chr.p.val$Position_SNP])
}

compute.correlation.blocks <- function(meth.data,annotation,cor.threshold=qtl.getOption("cluster.cor.threshold"),
                                       sd.gauss=qtl.getOption("standard.deviation.gauss"),
                                       absolute.cutoff=qtl.getOption("absolute.distance.cutoff")){
  require("RnBeads")
  require("psych")
  require("igraph")
  meth.data <- rnb.execute.imputation(meth.data)
  cor.all <- cor(t(meth.data))
  dist.all <- cor2dist(cor.all)
  dist.all[cor.all<cor.threshold] <- 0
  genomic.positions <- annotation$Start
  pairwise.distance <- abs(as.data.frame(lapply(genomic.positions,function(x)x-genomic.positions)))
  dist.all[pairwise.distance>absolute.cutoff] <- 0
  weighted.distances <- dist.all*dnorm(as.matrix(pairwise.distance),0,sd.gauss)
  colnames(weighted.distances) <- 1:ncol(weighted.distances)
  rownames(weighted.distances) <- 1:nrow(weighted.distances)
  graph.ad <- graph.adjacency(weighted.distances,weighted = T,mode = "undirected",diag=T)
  clust <- cluster_louvain(graph.ad)
  return(groups(clust))
}

call.methQTL.block <- function(cor.block,meth.data,sel.geno,covs){
  sel.meth <- meth.data[cor.block,]
  min.p.val <- apply(sel.meth,1,function(cpg){
    all.snps <- apply(sel.geno,1,function(snp){
      in.mat <- data.frame(CpG=cpg,SNP=snp,covs)
      if(is.null(covs)){
        form <- "CpG~SNP"
      }else{
        form <- paste("CpG~SNP",colnames(covs),sep="+")
      }
      desgn.mat <- model.matrix(form,data = in.mat)
      lm.model <- lm(form,desgn.mat)
    })
    is.min <- which.min(all.snps$p.val)
    c(all.snps$p.val[is.min],all.snps$beta[is.min],is.min)
  })
  min.p.val <- data.frame(min.p.val,1:nrow(min.p.val))
  colnames(min.p.val) <- c("P.value","Beta","Position_SNP","Position_CpG")
  min.val <- min(min.p.val$Min,na.rm=T)
  return(min.p.val[min.val,])
}
