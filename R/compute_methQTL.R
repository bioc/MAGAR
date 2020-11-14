##########################################################################################
# compute_methQTL.R
# created: 2019-08-29
# creator: Michael Scherer
# ---------------------------------------------------------------------------------------
# Methods to call methQTL from DNA methylation and genotyping data.
##########################################################################################

#' doMethQTL
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
#' @param default.options Flag indicating if default options for \code{cluster.cor.threshold},
#'          \code{standard.deviation.gauss}, and \code{absolute.distance.cutoff} should be loaded for the
#'          data set used. See the option settings in \code{'inst/extdata'}.
#' @return An object of type \code{\link{methQTLResult-class}} containing the called methQTL interactions.
#' @details The process is split into 4 steps:
#'          \describe{
#'            \item{1}{First the two matrices are split according to the chromosomes.}
#'            \item{2}{We then compute correlations among the CpGs and compute CpG correlation blocks.}
#'            \item{3}{In each of the CpG correlation blocks, linear models according to the \code{linear.model.type}
#'              \code{\link{qtlSetOption}} with the CpG methylation state of the reference CpG specified by
#'              \code{representative.cpg.computation} as output and the SNP genotype state and all possible covariates
#'              as input are computed.}
#'            \item{4}{For each of the CpG correlation blocks, we report the p-value of the representative CpG.}
#'          }
#' @seealso doMethQTLChromosome
#' @import methods
#' @author Michael Scherer
#' @export
doMethQTL <- function(meth.qtl,
                       sel.covariates=NULL,
                       p.val.cutoff=1e-5,
                       ncores=1,
                       cluster.submit=F,
                       out.dir=getwd(),
                       default.options=TRUE){
  if(!inherits(meth.qtl,"methQTLInput")){
    stop("Invalid value for meth.qtl, needs to be of type methQTLInput")
  }
  if(default.options){
    logger.info("Loading default option setting")
    if(meth.qtl@platform %in% "CpG"){
      logger.info("Keeping default setting for sequencing based assays")
    }else{
      if(meth.qtl@platform %in% "probes27"){
        stop("This package does not support Illumina Infinium 27k arrays.")
      }
      qtlJSON2options(file.path(system.file("extdata/",package="MAGAR"),paste0("qtl_options_",meth.qtl@platform,".json")))
    }
  }
  if(!meth.qtl@imputed){
    meth.qtl <- imputeMeth(meth.qtl)
  }
  all.chroms <- unique(getAnno(meth.qtl)$Chromosome)
  res.all <- list()
  logger.start("Computing methQTLs")
  if(!cluster.submit){
#    if(ncores>1){
#      parallel.setup(ncores)
#      res.all <- foreach(chrom=all.chroms,.combine="c") %dopar%{
#        doMethQTLChromosome(meth.qtl,chrom,sel.covariates,p.val.cutoff)
#      }
#    }else{
      for(chrom in all.chroms){
        res.chrom <- doMethQTLChromosome(meth.qtl,chrom,sel.covariates,p.val.cutoff,out.dir,ncores=ncores)
        #res.all[[chrom]] <- res.chrom
	meth.qtl.path <- file.path(out.dir,paste0("methQTLResult_",chrom))
	saveMethQTLResult(res.chrom,meth.qtl.path)
	rm(res.chrom)
	gc()
	res.all[[chrom]] <- meth.qtl.path
      }
#    }
    res.all <- lapply(res.all,loadMethQTLResult)
    res.all <- joinMethQTLResult(res.all)
  }else{
    res.all <- submitClusterJobs(meth.qtl,sel.covariates,p.val.cutoff,out.dir,ncores = ncores)
  }
  logger.completed()
  return(res.all)
}

#' doMethQTLChromosome
#'
#' This functions computes the methQTL interactions for a single chromosome
#'
#' @param meth.qtl An Object of type \code{\link{methQTLInput-class}}.
#' @param chrom Character vector represeting the chromosome to be investigated.
#' @param sel.covariates Covariates as column names of the sample annotation sheet stored in \code{meth.qtl} to be
#'          used for covariate adjustment.
#' @param p.val.cutoff The p-value used for methQTL calling
#' @param out.dir Optional argument specifying the output directory
#' @param ncores The number of cores to be used
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
#' @seealso doMethQTL
#' @author Michael Scherer
#' @export
#' @import doParallel
doMethQTLChromosome <- function(meth.qtl,chrom,sel.covariates=NULL,p.val.cutoff=1e-5,out.dir=NULL,ncores=1){
  logger.start(paste("Computing methQTL for chromosome",chrom))
  anno <- getAnno(meth.qtl,"meth")
  sel.meth <- which(anno$Chromosome %in% chrom)
  sel.anno <- anno[sel.meth,]
  sel.meth <- getMethData(meth.qtl)[sel.meth,]
  if(qtlGetOption("hdf5dump")){
    sel.meth <- writeHDF5Array(sel.meth)
  }
  if(qtlGetOption('compute.cor.blocks')){
    cor.blocks <- computeCorrelationBlocks(sel.meth,sel.anno,assembly=meth.qtl@assembly,chromosome=chrom,segmentation=meth.qtl@segmentation)
    cor.blocks <- lapply(cor.blocks,as.numeric)
    if(!is.null(out.dir)){
      to.plot <- data.frame(Size=lengths(cor.blocks))
      plot <- ggplot(to.plot,aes(x=Size,y=..count..))+geom_histogram(binwidth = 1)+geom_vline(xintercept = mean(to.plot$Size,na.rm=T))+
        theme_bw()+theme(panel.grid=element_blank(),text=element_text(size=18,color="black"),
                                                                                         axis.ticks=element_line(color="black"),plot.title = element_text(size=18,color="black",hjust = .5),
                                                                                         axis.text = element_text(size=15,color="black"))
      ggsave(file.path(out.dir,paste0("CpG_cluster_sizes_",chrom,".pdf")),plot)
    }
  }else{
    cor.blocks <- as.list(1:nrow(sel.meth))
  }
  anno.geno <- getAnno(meth.qtl,"geno")
  sel.geno <- which(anno.geno$Chromosome %in% chrom)
  sel.anno.geno <- anno.geno[sel.geno,]
  sel.geno <- getGeno(meth.qtl)[sel.geno,]
  ph.dat <- getPheno(meth.qtl)
  n.comps <- qtlGetOption('n.prin.comp')
  if(!is.null(n.comps)){
    sel.covariates <- c(sel.covariates,paste0("PC",1:n.comps))
  }
  if(is.null(sel.covariates)){
    ph.dat <- NULL
  }else{
    logger.info(paste("Using covariates",paste(sel.covariates,collapse="")))
    if(!all(sel.covariates %in% colnames(ph.dat))){
      logger.warning("Not all the covariates are present in the sample annotation sheet.")
      sel.covariates <- sel.covariates[sel.covariates %in% colnames(ph.dat)]
      if(is.null(sel.covariates)){
        stop("No valid covariate present")
      }
    }
    ph.dat <- ph.dat[,sel.covariates,drop=FALSE]
  }
  if(qtlGetOption("meth.qtl.type")=="fastQTL"){
    logger.start("Running FastQTL")
    prep.fast.qtl <- generateFastQTLInput(meth.qtl,chrom,cor.blocks,ph.dat,out.dir=out.dir)
    res.chr.p.val <- runFastQTL(prep.fast.qtl,meth.qtl,chrom,out.dir)
    logger.completed()
  }else{
    logger.start("Compute methQTL per correlation block")
      parallel.setup(ncores)
      res.chr.p.val <- mclapply(cor.blocks,callMethQTLBlock,sel.meth,sel.geno,ph.dat,sel.anno,sel.anno.geno,mc.cores = ncores)
      res.all <- c()
      for(i in 1:length(res.chr.p.val)){
        res.all <- rbind(res.all,res.chr.p.val[[i]])
      }
      res.chr.p.val <- as.data.frame(res.all)
      rm(res.all)
    logger.completed()
  }
  if(qtlGetOption("p.value.correction")=="uncorrected.fdr"){
    res.chr.p.val <- res.chr.p.val[as.numeric(as.character(res.chr.p.val$P.value))<p.val.cutoff,]
  }
  if(is.null(res.chr.p.val)){
    logger.info(paste("No methQTL found for chromosome",chrom))
    chrom.frame <- data.frame()
  }else if(nrow(res.chr.p.val)==0){
    logger.info(paste("No methQTL found for chromosome",chrom))
    chrom.frame <- data.frame()
  }else{
    chrom.frame <- data.frame(CpG=as.character(res.chr.p.val$CpG),
                            SNP=as.character(res.chr.p.val$SNP),
                            Beta=as.numeric(as.character(res.chr.p.val$Beta)),
			    SE.Beta=as.numeric(as.character(res.chr.p.val$SE.Beta)),
                            P.value=as.numeric(as.character(res.chr.p.val$P.value)),
                            Chromosome=as.character(res.chr.p.val$Chromosome),
                            Position.CpG=as.numeric(as.character(res.chr.p.val$Position_CpG)),
                            Position.SNP=as.numeric(as.character(res.chr.p.val$Position_SNP)))
    chrom.frame$Distance <- chrom.frame$Position.CpG - chrom.frame$Position.SNP
    if(qtlGetOption("p.value.correction") == "uncorrected.fdr"){
      tests.performed <- length(cor.blocks)*nrow(sel.geno)
      if(is.na(tests.performed)){
        tests.performed <- .Machine$integer.max
      }
      chrom.frame$p.val.adj.fdr <- p.adjust(chrom.frame$P.value,method="fdr",n=tests.performed)
    }
    meth.qtl.id <- paste(chrom.frame$CpG,chrom.frame$SNP,sep="_")
    match.unique <- match(unique(meth.qtl.id),meth.qtl.id)
    chrom.frame <- chrom.frame[match.unique,]
  }
  methQTL.result <- new("methQTLResult",
                        result.frame=chrom.frame,
                        anno.meth=sel.anno,
                        anno.geno=sel.anno.geno,
                        correlation.blocks=cor.blocks,
                        method=qtlGetOption("linear.model.type"),
                        rep.type=qtlGetOption("representative.cpg.computation"),
                        chr=as.character(chrom))
  logger.completed()
  return(methQTL.result)
}

#' callMethQTLBlock
#'
#' This function computes a methQTL per correlation block.
#' @param cor.block A single correlation block as determined by \code{\link{computeCorrelationBlocks}}
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
#' @param n.permutations  Number of permutations used to correct the p-values. Only has an influence, if the parameter
#'            \code{"p.value.correction"="permutation".}
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
callMethQTLBlock <- function(cor.block,meth.data,geno.data,covs,anno.meth,anno.geno,
                               model.type=qtlGetOption("linear.model.type"),
                               n.permutations=qtlGetOption("n.permutations")){
  # computing CpG-wise medians across the samples
  # then order the CpGs by value and select the one in the middle
  reps <- computeRepresentativeCpG(cor.block,meth.data,anno.meth)
  sel.anno <- reps$anno
  reps <- as.vector(reps$meth)
  distances <- abs(anno.geno$Start - sel.anno$Start)
  if(all(distances>=qtlGetOption("absolute.distance.cutoff"))){
    logger.info(paste("No SNP closer than",qtlGetOption("absolute.distance.cutoff")))
    ret <- data.frame(as.character(row.names(sel.anno)),
                      NA,
                      NA,
		      NA,
                      NA,
                      as.character(sel.anno$Chromosome),
                      NA,
                      sel.anno$Start)
    colnames(ret) <- c("CpG","SNP","P.value","Beta","SE.Beta","Chromosome","Position_SNP","Position_CpG")
    return(ret)
  }
  geno.data <- geno.data[distances<qtlGetOption("absolute.distance.cutoff"),,drop=FALSE]
  if(is.null(geno.data) || all(is.na(geno.data))){
    logger.error("Geno data is null")
    ret <- data.frame(as.character(row.names(sel.anno)),
                      NA,
                      NA,
		      NA,
                      NA,
                      as.character(sel.anno$Chromosome),
                      NA,
                      sel.anno$Start)
    colnames(ret) <- c("CpG","SNP","P.value","Beta","SE.Beta","Chromosome","Position_SNP","Position_CpG")
    return(ret)
  }else if(nrow(geno.data)==0){
    logger.error("Geno data contains no rows")
    ret <- data.frame(as.character(row.names(sel.anno)),
                      NA,
                      NA,
		      NA,
                      NA,
                      as.character(sel.anno$Chromosome),
                      NA,
                      sel.anno$Start)
    colnames(ret) <- c("CpG","SNP","P.value","Beta","SE.Beta","Chromosome","Position_SNP","Position_CpG")
    return(ret)
  }
  anno.geno <- anno.geno[distances<qtlGetOption("absolute.distance.cutoff"),,drop=FALSE]
  all.snps <- apply(as.matrix(geno.data),1,function(snp){
    if(is.null(covs)){
      form <- as.formula("CpG~SNP")
    }else{
      form <- as.formula(paste0("CpG~SNP+",paste0(colnames(covs),collapse="+")))
    }
    if(model.type == "categorical.anova"){
      if(length(unique(snp))==1){
        return(c(p.val=NA,beta=NA,se.beta=NA))
      }
      in.mat <- data.frame(CpG=reps,SNP=as.factor(snp))
      if(!is.null(covs)){
        in.mat <- data.frame(in.mat,covs)
      }
      lm.model <- lm(form,data=in.mat)
      if(length(unique(snp))==2){
        if(is.na(coef(lm.model)["SNP"])){
          return(c(p.val=NA,beta=NA,se.beta=NA))
        }
      }else{
        if(any(is.na(coef(lm.model)["SNP1"]),is.na(coef(lm.model)["SNP2"]))){
          return(c(p.val=NA,beta=NA,se.beta=NA))
        }
      }
      an.model <- anova(lm.model)
      p.val <- an.model["SNP","Pr(>F)"]
      return(c(p.val=p.val,beta=NA,se.beta=NA))
    }else if(model.type == "classical.linear"){
      in.mat <- data.frame(CpG=reps,SNP=snp)
      if(!is.null(covs)){
        in.mat <- data.frame(in.mat,covs)
      }
      lm.model <- lm(form,data=in.mat)
      if(is.na(coef(lm.model)["SNP"])){
        return(c(p.val=NA,beta=NA,se.beta=NA))
      }
      p.val <- summary(lm.model)$coefficients["SNP","Pr(>|t|)"]
      beta <- summary(lm.model)$coefficients["SNP","Estimate"]
      se.beta <- summary(lm.model)$coefficients["SNP","Std. Error"]
      if(qtlGetOption("p.value.correction") == "corrected.fdr"){
        permuted.pvals <- matrix(nrow = n.permutations,ncol=2)
        # determine number of independent tests
        for(i in 1:n.permutations){
          in.mat[,1] <- in.mat[sample(1:nrow(in.mat),nrow(in.mat)),1]
          lm.model <- lm(form,data=in.mat)
          p.value <- summary(lm.model)$coefficients["SNP","Pr(>|t|)"]
          p.value.adj <- p.value*n.permutations
          if(p.value.adj>=1) p.value.adj <- 0.99999999999
          p.value.adj <- log10(1-(p.value.adj))
          permuted.pvals[i,] <- c(log10(1-p.value),p.value.adj)
        }
        mod <- lm(permuted.pvals[,1]~permuted.pvals[,2])
        if(is.na(coef(mod)[2])){
          n.tests <- 1
        }else{
          n.tests <- round(summary(mod)$coefficients[2,"Estimate"])
        }
        p.val <- p.adjust(p.val,"bonferroni",n=ifelse(n.tests<1,1,n.tests))
      }
     return(c(p.val=p.val,beta=beta,se.beta=se.beta))
    }
  })
  all.snps <- t(all.snps)
  if(qtlGetOption("meth.qtl.type")%in%"oneVSall"){
    is.min <- which.min(all.snps[,'p.val'])
  }else if(qtlGetOption("meth.qtl.type")%in%"allVSall"){
    is.min <- 1:nrow(all.snps)
  }else if(qtlGetOption("meth.qtl.type")%in%"twoVSall"){
    all.inds <- rep(FALSE,nrow(all.snps))
    all.inds[all.snps[,"beta"]>0][which.min(all.snps[all.snps[,"beta"]>0,"p.val"])] <- TRUE
    all.inds[all.snps[,"beta"]<0][which.min(all.snps[all.snps[,"beta"]<0,"p.val"])] <- TRUE
    is.min <- which(all.inds)
  }
  if(length(is.min)==0){
    #logger.info(paste("No methQTL found for block",cor.block))
    ret <- data.frame(as.character(row.names(sel.anno)),
                      NA,
                      NA,
		      NA,
                      NA,
                      as.character(sel.anno$Chromosome),
                      NA,
                      sel.anno$Start)
    colnames(ret) <- c("CpG","SNP","P.value","Beta","SE.Beta","Chromosome","Position_SNP","Position_CpG")
    return(ret)
  }
  min.p.val <- data.frame(as.character(row.names(sel.anno)),
                          as.character(row.names(anno.geno)[is.min]),
                          all.snps[is.min,'p.val'],
                          all.snps[is.min,'beta'],
			  all.snps[is.min,'se.beta'],
                          as.character(sel.anno$Chromosome),
                          anno.geno$Start[is.min],
                          sel.anno$Start)
  colnames(min.p.val) <- c("CpG","SNP","P.value","Beta","SE.Beta","Chromosome","Position_SNP","Position_CpG")
  return(min.p.val)
}

#' computeRepresentativeCpG
#'
#' This function computes respresentative CpGs per correlation blocks given and return the methylation data matrix
#' and genomic annotation of these representatives.
#'
#' @param cor.blocks The correlation blocks as a \code{list}
#' @param meth.data The originial DNA methylation data matrix
#' @param anno.meth The original genomic annotation of the CpGs
#' @return A list with two elements: \describe{
#'     \item{\code{meth:}}{The selected DNA methylation data matrix}
#'     \item{\code{anno:}}{The selexted genomic annotation}
#' }
#' @author Michael Scherer
#' @noRd
computeRepresentativeCpG <- function(cor.blocks,meth.data,annotation){
  res.meth <- c()
  res.anno <- c()
  repr.type <- qtlGetOption("representative.cpg.computation")
  if(is.list(cor.blocks)){
	  for(i in 1:length(cor.blocks)){
	    cor.block <- cor.blocks[[i]]
	    sel.meth <- meth.data[cor.block,,drop=F]
	    anno.meth <- annotation[cor.block,,drop=F]
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
	      row.names(sel.anno) <- paste0("mean_of_",nrow(sel.meth))
	    }else{
	      sel.anno <- anno.meth
	      reps <- t(as.matrix(sel.meth))
	    }
	    res.meth <- rbind(res.meth,reps)
	    res.anno <- rbind(res.anno,sel.anno)
	  }
 	 return(list(meth=res.meth,anno=res.anno))
  }else{
	    sel.meth <- meth.data[cor.blocks,,drop=FALSE]
	    anno.meth <- annotation[cor.blocks,,drop=FALSE]
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
	      row.names(sel.anno) <- paste0("mean_of_",nrow(sel.meth))
	    }else{
	      sel.anno <- anno.meth
	      reps <- t(as.matrix(sel.meth))
	    }
  	return(list(meth=reps,anno=sel.anno))
  }
}

#' runFastQTL
#'
#' This functions runs fastQTL on the specified input prepared using \code{generate.fastQTL.input}
#'
#' @param prepard.input The input as prepared through \code{generate.fastQTL.input}
#' @param meth.qtl The input object of type \code{\link{methQTLInput}}
#' @param chrom The chromosome to be analyzed
#' @param out.dir The output directory
#' @return A \code{data.frame} in analogy to \code{\link{callMethQTLBlock}}
#' @author Michael Scherer
#' @noRd
runFastQTL <- function(prepared.input,
                        meth.qtl,
                        chrom,
                        out.dir){
  fastQTL.path <- qtlGetOption("fast.qtl.path")
  n.permutations <- qtlGetOption("n.permutations")
  res.file <- file.path(out.dir,paste0("fastQTLresult_",chrom,".txt.gz"))
  cmd <- paste(fastQTL.path,
               "--vcf", prepared.input["genotypes"],
               "--bed", prepared.input["phenotypes"],
               "--permute", n.permutations,
               "--window",qtlGetOption("absolute.distance.cutoff"),
               "--out", res.file,
	       "--chunk 1 1",
	       "--silent")
  if(!is.na(prepared.input["covariates"])){
    cmd <- paste(cmd,"--cov",prepared.input["covariates"])
  }
  system(cmd)
  res.frame <- read.table(res.file,header = F,sep=" ")
  pos.cpg <- getAnno(meth.qtl)[as.character(res.frame[,1]),"Start"]
  pos.snp <- getAnno(meth.qtl,"geno")[as.character(res.frame[,6]),"Start"]
  ret <- data.frame(as.character(res.frame[,1]),
                    as.character(res.frame[,6]),
                    res.frame[,11],
                    res.frame[,9],
                    rep(NA,nrow(res.frame)),
                    rep(chrom,nrow(res.frame)),
                    pos.snp,
                    pos.cpg)
  colnames(ret) <- c("CpG","SNP","P.value","Beta","SE.Beta","Chromosome","Position_SNP","Position_CpG")
  return(ret)
}
