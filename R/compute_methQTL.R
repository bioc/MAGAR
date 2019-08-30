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
#' @return An object of type \code{\link{methQTLResult-class}} containing the called methQTL interactions.
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
#' @param sel.covariates Covariates as column names of the sample annotation sheet stored in \code{meth.qtl} to be
#'          used for covariate adjustment.
#' @param p.val.cutoff The p-value used for methQTL calling
#' @return A data frame with seven columns:
#'        \describe{
#'          \item{CpG}{The CpG ID (as cgNNNNNN) involved in the methQTL}
#'          \item{SNP}{The SNP ID (as rsNNNNNN) involved in the methQTL}
#'          \item{Beta}{The coefficient estimate of the linear model}
#'          \item{P.value}{The p-value associated with the coefficient estimate}
#'          \item{Chromosome}{The chromosome name}
#'          \item{Position.CpG}{The genomic position of the CpG}
#'          \item{Position.SNP}{The genomic position of the SNP}
#'          \item{Distance}{The distance between the CpG and the SNP}
#'        }
#' @author Michael Scherer
#' @noRd
do.methQTL.chromosome <- function(meth.qtl,chrom,sel.covariates,p.val.cutoff){
  anno <- getAnno(meth.qtl,"meth")
  sel.meth <- which(anno$Chromosome %in% chrom)
  sel.anno <- anno[sel.meth,]
  sel.meth <- getMethData(meth.qtl)[sel.meth,]
  cor.blocks <- compute.correlation.blocks(sel.meth,sel.anno)
  cor.blocks <- lapply(cor.blocks,as.numeric)
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
    }
    ph.dat <- ph.dat[,sel.covariates,drop=FALSE]
  }
  res.chr.p.val <- sapply(cor.blocks,call.methQTL.block,sel.meth,sel.geno,ph.dat)
  res.chr.p.val <- res.chr.p.val[res.chr.p.val$P.value<p.val.cutoff,]
  tests.performed <- length(cor.blocks)*nrow(sel.geno)
  chrom.frame <- data.frame(CpG=row.names(anno)[res.chr.p.val$Position_CpG],
                            SNP=row.names(anno.geno)[res.chr.p.val$Position_SNP],
                            Beta=res.chr.p.val$Beta,
                            P.value=res.chr.p.val$P.value,
                            Chromosome=anno$Chromosome[res.chr.p.val$Position_CpG],
                            Position.CpG=anno$Start[res.chr.p.val$Position_CpG],
                            Position.SNP=anno$Start[res.chr.p.val$Position_SNP])
  chrom.frame$Distance <- abs(chrom.frame$Position.CpG - chrom.frame$Position.SNP)
  return(chrom.frame)
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

#' call.methQTL.block
#'
#' This function computes a methQTL per correlation block.
#' @param cor.block A single correlation block as determined by \code{\link{compute.correlation.blocks}}
#' @param meth.data The methylation data matrix
#' @param geno.data The genotype data matrix
#' @param covs A \code{data.frame} specifying the covariates used for the analysis
#' @param model.type Type of the linear model to be used. Supported options are:
#'    \describe{
#'      \item{categorical.anova}{Performs linear regression after transforming the genotype states (0=homozygous
#'      reference allele, 1=heteroygous, 2=homoygous alternate allele) into categorical variables and then perform
#'      ANOVA}
#'      \item{classical.linear}{Performs linear regression with keeping the genotype calls as numeric values. This
#'      induces an ordering between heterozygous and homozygous reference. We then use the p-value from the t-statistic
#'      and the associated beta value.}
#'    }
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
call.methQTL.block <- function(cor.block,meth.data,geno.data,covs,model.type=qtl.getOption("linear.model.type")){
  sel.meth <- meth.data[cor.block,]
  min.p.val <- apply(sel.meth,1,function(cpg){
    all.snps <- apply(geno.data,1,function(snp){
      if(is.null(covs)){
        form <- as.formula("CpG~SNP")
      }else{
        form <- as.formula(paste("CpG~SNP",colnames(covs),sep="+"))
      }
      if(model.type == "categorical.anova"){
        in.mat <- data.frame(CpG=cpg,SNP=as.factor(snp),covs)
        lm.model <- lm(form,data=in.mat)
        an.model <- anova(lm.model)
        p.val <- an.model["SNP","Pr(>F)"]
        return(p.val=p.val,beta=NA)
      }else if(model.type == "classical.linear"){
        in.mat <- data.frame(CpG=cpg,SNP=snp,covs)
        lm.model <- lm(form,data=in.mat)
        p.val <- summary(lm.model)$coefficients["SNP","Pr(>|t|)"]
        beta <- summary(lm.model)$coefficients["SNP","Estimate"]
        return(c(p.val=p.val,beta=beta))
      }
    })
    all.snps <- t(all.snps)
    is.min <- which.min(all.snps[,'p.val'])
    c(all.snps[is.min,'p.val'],all.snps[is.min,'beta'],is.min)
  })
  min.p.val <- t(min.p.val)
  min.p.val <- data.frame(min.p.val,1:nrow(min.p.val))
  colnames(min.p.val) <- c("P.value","Beta","Position_SNP","Position_CpG")
  min.val <- which.min(min.p.val$P.value)
  return(min.p.val[min.val,])
}
