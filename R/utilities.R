##########################################################################################
# utilities.R
# created: 2020-01-31
# creator: Michael Scherer
# ---------------------------------------------------------------------------------------
# utitiliy functions
##########################################################################################

#' overlapQTLs
#'
#' This function overlaps a list of methQTLs to determine which interactions are common.
#'
#' @param meth.qtl.result.list A named list with each entry being an object of type \code{\link{methQTLResult-class}}.
#'                       The names are used in the visualization.
#' @param type Determines if either the SNP (default), the CpG, or the correlation block
#'     \code{'cor.block'} is to be visualized
#' @return A list with \code{length(meth.qtl.result.list)} elements, containing IDs of methQTL
#'     interactions according to the option \code{type}.
#' @author Michael Scherer
#' @export
#' @examples {
#'  meth.qtl.res.1 <- loadMethQTLResult(system.file("extdata/methQTLResult_chr18",package="MAGAR"))
#'  meth.qtl.res.2 <- meth.qtl.res.1
#'  res <- overlapQTLs(list(A=meth.qtl.res.1,B=meth.qtl.res.2),type="SNP")
#' }
overlapQTLs <- function(meth.qtl.result.list,type){
  if(!type %in% c("CpG","SNP",'cor.block')){
    stop("Invalid value for type, needs to be 'CpG', 'SNP', or 'cor.block'")
  }
  if(type%in%c("CpG","SNP")){
    res.all <- lapply(meth.qtl.result.list,function(res){
      getResult(res)[,type]
    })
  }else{
#    res.all.vec <- c()
#    res.all <- list()
#    logger.start("Constructing universe")
#    for(i in 1:length(meth.qtl.result.list)){
#      logger.start(paste("Obtaining correlation blocks for object",i))
#      cor.blocks <- getCorrelationBlocks(meth.qtl.result.list[[i]])
#      res <- getResult(meth.qtl.result.list[[i]],cor.blocks)
#      for(j in 1:nrow(res)){
#        snp <- res$SNP[j]
#        cor.block <- res$CorrelationBlock[j]
#        map.qtl <- unlist(sapply(cor.block[[1]],function(cpg){
#          matched <- unlist(grepl(cpg,res.all.vec) & grepl(snp,res.all.vec))
#          if(any(matched)){
#            return(which(matched))
#          }else{
#            return(NULL)
#          }
#        }))
#        if(!is.null(map.qtl)){
#          extended.qtl <- paste(snp,paste(cor.block[[1]],collapse = "_"),sep="_")
#	  if(extended.qtl==res.all.vec[map.qtl[1]]){
#		next
#	  }else{
#		old.cpgs <- unlist(strsplit(res.all.vec[map.qtl[1]],"_"))[-1]
#		res.all.vec[map.qtl[[1]]] <-  paste(snp,paste(union(old.cpgs,cor.block[[1]]),collapse = "_"),sep="_")
#	  }
#        }else{
#          #Caveat: only first interaction of a SNP with a correlation block will be considered, thus 'trans' effects
#          #are not considered
#          res.all.vec <- c(res.all.vec,paste(snp,paste(cor.block[[1]],collapse = "_"),sep="_"))
#        }
#      }
#      logger.completed()
#    }
#    logger.completed()
    res.all.vec <- list()
    res.all <- list()
    logger.start("Constructing universe")
    for(i in 1:length(meth.qtl.result.list)){
      logger.start(paste("Obtaining correlation blocks for object",i))
      cor.blocks <- getCorrelationBlocks(meth.qtl.result.list[[i]])
      res <- getResult(meth.qtl.result.list[[i]],cor.blocks)
      for(j in 1:nrow(res)){
        snp <- as.character(res$SNP[j])
        cor.block <- res$CorrelationBlock[j]
	if(!(snp%in%names(res.all.vec))){
		res.all.vec[[snp]] <- cor.block[[1]]
	}else{
		map.qtl <- res.all.vec[[snp]]
		matched <- sapply(cor.block[[1]],function(cpg)grepl(cpg,map.qtl))
	        if(any(matched)){
            		if(all(cor.block[[1]]%in%map.qtl)&all(map.qtl%in%cor.block[[1]])){
				next
			}else{
				res.all.vec[[snp]] <- union(map.qtl,cor.block[[1]])
			}
          	}
	}
      }
      logger.completed()
    }
    logger.start("Overlapping")
    for(i in 1:length(meth.qtl.result.list)){
      logger.start(paste("Obtaining correlation blocks for object",i))
      cor.blocks <- getCorrelationBlocks(meth.qtl.result.list[[i]])
      res <- getResult(meth.qtl.result.list[[i]],cor.blocks)
      res.all.class <- list()
      for(j in 1:nrow(res)){
        snp <- as.character(res$SNP[j])
        cor.block <- res$CorrelationBlock[j]
	if(snp%in%names(res.all.vec)){
		map.qtl <- res.all.vec[[snp]]
		matched <- sapply(cor.block[[1]],function(cpg)grepl(cpg,map.qtl))
		if(is.null(matched)){
	  		stop("Constructing the universe failed")
		}
		res.all.class[[snp]] <- map.qtl
	}else{
		stop("Constructing the universe failed")
	}
      }
      logger.completed()
      res.class.vec <- c()
      for(j in 1:length(res.all.class)){
        res.class.vec <- c(res.class.vec,
                           paste(names(res.all.class)[j],
                                 paste(res.all.class[[j]],collapse="_"),sep="_"))
      }
      res.all[[i]] <- res.class.vec
    }
    logger.completed()
  }
  names(res.all) <- names(meth.qtl.result.list)
  return(res.all)
}

#' overlapInputs
#'
#' Overlaps the input annotations
#'
#' @param meth.qtl.list A list of \code{\link{methQTLInput-class}} or \code{\link{methQTLResult-class}} objects to be overlapped
#' @param type The type of annotation to be overlapped. Needs to be \code{'SNP'}, \code{'CpG'} or \code{'cor.block'}
#' @return A data frame containing the annotations of the unique input values.
#' @author Michael Scherer
#' @export
#' @examples {
#'  meth.qtl.1 <- loadMethQTL(system.file("extdata/reduced_methQTL",package="MAGAR"))
#'  meth.qtl.2 <- meth.qtl.1
#'  res <- overlapInputs(list(A=meth.qtl.res.1,B=meth.qtl.res.2),type="SNP")
#' }
overlapInputs <- function(meth.qtl.list,type){
  if(!(inherits(meth.qtl.list[[1]],"methQTLResult")|inherits(meth.qtl.list[[1]],"methQTLInput"))){
    stop("Invalid value for meth.qtl.list, needs to be methQTLResult")
  }
  if(!(type%in%c("CpG","SNP","cor.block"))){
    stop("Invalid value for type, needs to be 'CpG', 'SNP', or 'cor.block'")
  }
  type <- ifelse(type%in%c("CpG","cor.block"),"meth","geno")
  all.input <- getAnno(meth.qtl.list[[1]],type)
  for(i in 2:length(meth.qtl.list)){
    anno <- getAnno(meth.qtl.list[[i]],type)
    all.input <- rbind(all.input[,intersect(colnames(anno),colnames(all.input))],
                       anno[!(row.names(anno)%in%row.names(all.input)),
                            intersect(colnames(anno),colnames(all.input))])
  }
  all.input <- as.data.frame(all.input)
#  if(type%in%"SNP"){
#    all.input$End <- as.numeric(as.character(all.input$Start))
#    input.ranges <- makeGRangesFromDataFrame(all.input)
#    sel.input <- rep(FALSE,nrow(all.input))
#    all.cpgs <- makeGRangesFromDataFrame(overlap.inputs(meth.qtl.list,"CpG"))
#    for(i in 1:length(all.cpgs)){
#      sel.input[distance(all.cpgs[i],input.ranges)<qtl.getOption("absolute.distance.cutoff")] <- TRUE
#    }
#    all.input <- all.input[sel.input,]
#  }
  return(all.input)
}

#' getOverlapUniverse
#'
#' This function overlaps results from a list of methQTLResults and returns the union of all
#' the input data points used.
#'
#' @param meth.qtl.res An object of type \code{\link{methQTLResult-class}} or a list of such objects
#' @param type The type of annotation to be overlapped. Needs to be \code{'SNP'}, \code{'CpG'} or \code{'cor.block'}
#' @return A list with two GRanges objects, one containing the overlapped set and the other the
#'    union of input data points as elements \code{'all.qtl'} and \code{'all.input'}
#' @author Michael Scherer
#' @export
#' @examples {
#'  meth.qtl.res.1 <- loadMethQTLResult(system.file("extdata/methQTLResult_chr18",package="MAGAR"))
#'  meth.qtl.res.2 <- meth.qtl.res.1
#'  res <- getOverlapUniverse(list(A=meth.qtl.res.1,B=meth.qtl.res.2),type="SNP")
#' }
getOverlapUniverse <- function(meth.qtl.res,type){
  if(inherits(meth.qtl.res,"methQTLResult")){
    type.anno <- ifelse(type%in%c("CpG","cor.block"),"meth","geno")
    all.input <- getAnno(meth.qtl.res,type.anno)
#    if(type%in%"SNP"){
#      all.input$End <- as.numeric(as.character(all.input$Start))
#      input.ranges <- makeGRangesFromDataFrame(all.input)
#      sel.input <- rep(TRUE,nrow(all.input))
#      all.cpgs <- makeGRangesFromDataFrame(getAnno(meth.qtl.res,"meth"))
#      print(Sys.time())
#      for(i in 1:length(all.cpgs)){
#        if((i%%10000)==0)print(i)
#        pot.input <- seqnames(input.ranges)%in%seqnames(all.cpgs[i])
#        sel.input[as.logical(pot.input)][distance(all.cpgs[i],input.ranges[as.logical(pot.input)])<qtl.getOption("absolute.distance.cutoff")] <- TRUE
#      }
#      print(Sys.time())
#      all.input <- all.input[sel.input,]
#    }
    if(type%in%"cor.block"){
      cor.blocks <- getCorrelationBlocks(meth.qtl.res)
      all.qtl <- apply(getResult(meth.qtl.res,cor.blocks),1,function(ro){
        paste(ro["SNP"],
              paste(ro["CorrelationBlock"][[1]],collapse = "_"),sep="_")
      })
    }else{
      all.qtl <- unique(as.character(getResult(meth.qtl.res)[,type]))
    }
  }
  if(is.list(meth.qtl.res)){
    if(!inherits(meth.qtl.res[[1]],"methQTLResult")){
      stop("Invalid value for meth.qtl.res, needs to be methQTLResult")
    }
    all.input <- overlapInputs(meth.qtl.res,type=type)
    all.qtl <- overlapQTLs(meth.qtl.res = meth.qtl.res,type=type)
    res <- all.qtl[[1]]
    for(i in 2:length(all.qtl)){
      res <- intersect(res,all.qtl[[i]])
    }
    all.qtl <- res
  }
  if(type%in%"cor.block"){
    all.qtl <- unlist(sapply(all.qtl,function(qtl){
      strsplit(qtl,"_")
    }))
    all.qtl <- unique(all.qtl[!grepl("rs|[:]",all.qtl)])
  }
  if(type%in%"SNP"){
    all.input$End <- all.input$Start+1
  }
  all.qtl <- all.input[all.qtl,]
  all.input <- makeGRangesFromDataFrame(all.input)
  all.qtl <- makeGRangesFromDataFrame(all.qtl)
  return(list("all.qtl"=all.qtl,"all.input"=all.input))
}

#' getSpecificQTL
#'
#' This function returns the methQTL interactions specific for a result
#' @param meth.qtl.res An object of type \code{\link{methQTLResult-class}} for which specific QTLs are to be obtained.
#' @param meth.qtl.background The background set as a list of \code{\link{methQTLResult-class}} objects.
#' @param type The type of annotation to be overlapped. Needs to be \code{'SNP'}, \code{'CpG'} or \code{'cor.block'}
#' @return A \code{data.frame} of methQTL interactions sorted by the effect size.
#' @author Michael Scherer
#' @export
#' @examples {
#'  meth.qtl.res.1 <- loadMethQTLResult(system.file("extdata/methQTLResult_chr18",package="MAGAR"))
#'  meth.qtl.res.2 <- meth.qtl.res.1
#'  res <- getSpecificQTL(meth.qtl.res.1,list(A=meth.qtl.res.1,B=meth.qtl.res.2),type="SNP")
#' }
getSpecificQTL <- function(meth.qtl.res,
                           meth.qtl.background,
                           type="SNP"){
  if(!inherits(meth.qtl.res,"methQTLResult")){
    stop("Invalid value for meth.qtl.res, needs to be methQTLResult")
  }
  #shared.qtl <- unique(unlist(lapply(meth.qtl.background,function(qtl.obj){
  #  ret <- overlap.QTLs(c(meth.qtl.res,qtl.obj),type=type)
  #  intersect(ret[[1]],ret[[2]])
  #})))
  #if(type%in%"cor.block"){
  #  cor.blocks <- getCorrelationBlocks(meth.qtl.res)
  #  res <- getResult(meth.qtl.res,cor.blocks)
  #  type <- "CorrelationBlock"
  #  sel.qtl <- !(sapply(res[,type],function(cori){
  #    any(sapply(cori,function(qtl){
  #      grepl(qtl,shared.qtl)}))
  #  }))
  #}else{
  #  res <- getResult(meth.qtl.res)
  #  sel.qtl <- !(res[,type] %in% shared.qtl)
  #}
  #res <- res[sel.qtl,]
  op.qtl <- overlapQTLs(c(meth.qtl.res,meth.qtl.background),type=type)
  op.qtl <- lapply(op.qtl,as.character)
  my.qtls <- op.qtl[[1]]
  if(length(op.qtl)>2){
    other.qtls <- unique(do.call(c,op.qtl[-1]))
  }else{
    other.qtls <- op.qtl[[2]]
  }
  spec.qtls <- setdiff(my.qtls,other.qtls)
  res <- getResult(meth.qtl.res)
  if(type%in%"cor.block"){
    snp.ids <- unlist(lapply(strsplit(spec.qtls,"_"),function(x)x[1]))
    cpg.ids <- lapply(strsplit(spec.qtls,"_"),function(x)x[-1])
    sel.qtls <- apply(res,1,function(r){
      pot.qtl <- which(snp.ids%in%r["SNP"])
      if(length(pot.qtl)>0){
        r["CpG"]%in%unlist(cpg.ids[pot.qtl])
      }else{
        FALSE
      }
    })
    res <- res[sel.qtls,]
  }else{
    res <- res[res[,type]%in%spec.qtls,]
  }
  return(res[order(abs(res$P.value),decreasing=FALSE),])
}

#' getOverlappingQTL
#'
#' This function merges the QTLs given and returns the methQTL table in a merged format.
#'
#' @param meth.qtl.list A list of \code{\link{methQTLResult-class}} objects to be merged
#' @param type The type of annotation to be overlapped. Needs to be \code{'SNP'}, \code{'CpG'} or \code{'cor.block'}
#' @return A \code{data.frame} with the methQTL interactions and an additional column specifying where the interaction
#'    displayed has been found. This value is generated from the \code{names()} argument of \code{meth.qtl.list}.
#' @author Michael Scherer
#' @export
#' @examples {
#'  meth.qtl.res.1 <- loadMethQTLResult(system.file("extdata/methQTLResult_chr18",package="MAGAR"))
#'  meth.qtl.res.2 <- meth.qtl.res.1
#'  res <- getOverlappingQTL(list(A=meth.qtl.res.1,B=meth.qtl.res.2),type="SNP")
#' }
getOverlappingQTL <- function(meth.qtl.list,type="SNP"){
  op.qtl <- overlapQTLs(meth.qtl.list,type=type)
  all.ops <- op.qtl[[1]]
  for(i in 1:length(op.qtl)){
    all.ops <- intersect(all.ops,op.qtl[[i]])
  }
  res.all <- c()
  for(i in 1:length(meth.qtl.list)){
    qtl <- meth.qtl.list[[i]]
    res <- getResult(qtl)
    if(type%in%"cor.block"){
      snp.ids <- unlist(lapply(strsplit(all.ops,"_"),function(x)x[1]))
      cpg.ids <- lapply(strsplit(all.ops,"_"),function(x)x[-1])
      sel.qtls <- apply(res,1,function(r){
        pot.qtl <- which(snp.ids%in%r["SNP"])
        if(length(pot.qtl)>0){
          r["CpG"]%in%unlist(cpg.ids[pot.qtl])
        }else{
          FALSE
        }
      })
      res <- res[sel.qtls,]
    }else{
      res <- res[res[,type]%in%all.ops,]
    }
    res$Type <- rep(names(meth.qtl.list)[i],nrow(res))
    res.all <- rbind(res.all,res)
  }
  return(res.all)
}

#' generateFastQTLInput
#'
#' This function generates the required input to FastQTL and stores the files on disk.
#'
#' @param meth.qtl An object of type \code{\link{methQTLInput-class}} with methylation and genotyping information
#' @param chrom The chromosome to be investigated.
#' @param correlation.block The correlation blocks generated.
#' @param sel.covariates The selected covariates as a \code{data.frame}
#' @param out.dir The target output directory
#' @return A vector comprising the following elements:\describe{
#'           \item{\code{genotypes:}}{The path to the genotypes file}
#'           \item{\code{phenotypes:}}{The path to the DNA methylation data file}
#'           \item{\code{covariates:}}{The path to the covariates file}
#' }
#' @noRd
generateFastQTLInput <- function(meth.qtl,
                                 chrom,
                                 correlation.block,
                                 sel.covariates,
                                 out.dir){
  geno.data <- getGeno(meth.qtl)
  anno.geno <- getAnno(meth.qtl,"geno")
  geno.data <- geno.data[anno.geno$Chromosome%in%chrom,]
  anno.geno <- anno.geno[anno.geno$Chromosome%in%chrom,]
  vcf.frame <- data.frame(CHROM=anno.geno$Chromosome,
                          POS=anno.geno$Start,
                          ID=row.names(anno.geno),
                          REF=anno.geno$Allele.1,
                          ALT=anno.geno$Allele.2,
                          QUAL=rep(100,nrow(anno.geno)),
                          FILTER=rep("PASS",nrow(anno.geno)),
                          INFO=rep("INFO",nrow(anno.geno)),
			  FORMAT=rep("DS",nrow(anno.geno)),
                          geno.data)
  f.name <- file.path(out.dir,paste0("genotypes_",chrom,".vcf"))
  write.table(vcf.frame,
              f.name,
              sep="\t",
              row.names = FALSE,
              col.names = FALSE,
              quote=FALSE)
  system(paste("sed -i '1i#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
               paste(colnames(geno.data),collapse="\t"),"'",f.name))
  system(paste("sed -i '1i##fileformat=VCFv4.1'",f.name))
  system(paste(qtlGetOption("bgzip.path"),
               f.name,"&&",
               qtlGetOption("tabix.path"),
               "-p vcf",
               paste0(f.name,".gz")))
  meth.data <- getMethData(meth.qtl)
  anno.meth <- getAnno(meth.qtl)
  meth.data <- meth.data[anno.meth$Chromosome%in%chrom,]
  anno.meth <- anno.meth[anno.meth$Chromosome%in%chrom,]
  reps <- computeRepresentativeCpG(correlation.block,meth.data,anno.meth)
  anno.meth <- reps$anno
  meth.data <- reps$meth
  pheno.frame <- data.frame(Chr=anno.meth$Chromosome,
                            start=anno.meth$Start,
                            end=anno.meth$End,
                            ID=row.names(anno.meth),
                            meth.data)
  pheno.frame <- pheno.frame[order(pheno.frame$start,decreasing=FALSE),]
  f.name <- file.path(out.dir,paste0("phenotypes_",chrom,".bed"))
  write.table(pheno.frame,f.name,
              sep="\t",
              row.names = FALSE,
              col.names = FALSE,
              quote=FALSE)
  system(paste("sed -i '1i#Chr\tstart\tend\tID",
               paste(colnames(meth.data),collapse="\t"),"'",f.name))
  system(paste(qtlGetOption("bgzip.path"),
               f.name,"&&",
               qtlGetOption("tabix.path"),
               "-p bed",
               paste0(f.name,".gz")))
  f.name <- file.path(out.dir,"covariates.txt")
  if(is.null(sel.covariates)){
	cov.file <- NULL
  }else{
	cov.file <- file.path(out.dir,"covariates.txt")
	  cov.frame <- data.frame(id=meth.qtl@samples,
		                  sel.covariates)
	  write.table(t(cov.frame),f.name,sep="\t",row.names = TRUE,col.names=FALSE,quote = FALSE)
  }
  return(c(genotypes=file.path(out.dir,paste0("genotypes_",chrom,".vcf.gz")),
           phenotypes=file.path(out.dir,paste0("phenotypes_",chrom,".bed.gz")),
           covariates=cov.file))
}
