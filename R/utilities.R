##########################################################################################
# utilities.R
# created: 2020-01-31
# creator: Michael Scherer
# ---------------------------------------------------------------------------------------
# utitiliy functions
##########################################################################################

#' overlap.QTLs
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
overlap.QTLs <- function(meth.qtl.result.list,type){
  if(!type %in% c("CpG","SNP",'cor.block')){
    stop("Invalid value for type, needs to be 'CpG', 'SNP', or 'cor.block'")
  }
  if(type%in%c("CpG","SNP")){
    res.all <- lapply(meth.qtl.result.list,function(res){
      getResult(res)[,type]
    })
  }else{
    res.all.vec <- c()
    res.all <- list()
    for(i in 1:length(meth.qtl.result.list)){
      logger.start(paste("Obtaining correlation blocks for object",i))
      cor.blocks <- getCorrelationBlocks(meth.qtl.result.list[[i]])
      res <- getResult(meth.qtl.result.list[[i]],cor.blocks)
      logger.completed()
      res.all.class <- c()
      for(j in 1:nrow(res)){
        snp <- res$SNP[j]
        cor.block <- res$CorrelationBlock[j]
        map.qtl <- unlist(sapply(cor.block[[1]],function(cpg){
          matched <- unlist(grepl(cpg,res.all.vec) & grepl(snp,res.all.vec))
          if(any(matched)){
            return(which(matched))
          }else{
            return(NULL)
          }
        }))
        if(!is.null(map.qtl)){
          res.all.class <- c(res.all.class,res.all.vec[map.qtl[1]])
          res.all.vec <- c(res.all.vec,res.all.vec[map.qtl[1]])
        }else{
          res.all.class <- c(res.all.class,paste(snp,paste(cor.block[[1]],collapse = "_"),sep="_"))
          res.all.vec <- c(res.all.vec,paste(snp,paste(cor.block[[1]],collapse = "_"),sep="_"))
        }
      }
      res.all[[i]] <- res.all.class
    }
  }
  names(res.all) <- names(meth.qtl.result.list)
  return(res.all)
}

#' overlap.inputs
#'
#' Overlaps the input annotations
#'
#' @param meth.qtl.list A list of \code{\link{methQTLInput}} or \code{\link{methQTLResult}} objects to be overlapped
#' @param type The type of annotation to be overlapped. Needs to be \code{'SNP'}, \code{'CpG'} or \code{'cor.block'}
#' @return A data frame containing the annotations of the unique input values.
#' @author Michael Scherer
#' @export
overlap.inputs <- function(meth.qtl.list,type){
  if(!(inherits(meth.qtl.list[[1]],"methQTLResult")|inherits(meth.qtl.list[[1]],"methQTLInput"))){
    stop("Invalid value for meth.qtl.list, needs to be methQTLResult")
  }
  if(!(type%in%c("CpG","SNP","cor.block"))){
    stop("Invalid value for type, needs to be 'CpG', 'SNP', or 'cor.block'")
  }
  type <- ifelse(type%in%c("CpG","cor.block"),"meth","geno")
  all.input <- c()
  for(i in 1:length(meth.qtl.list)){
    anno <- getAnno(meth.qtl.list[[1]],type)
    all.input <- cbind(all.input,anno[!(row.names(anno)%in%row.names(all.input)),])
  }
  all.input <- as.data.frame(all.input)
  if(type%in%"SNP"){
    all.input$End <- as.numeric(as.character(all.input$Start))
  }
  return(all.input)
}

#' get.overlap.universe
#'
#' This function overlaps results from a list of methQTLResults and returns the union of all
#' the inpute data points used.
#'
#' @param meth.qtl.res An objec of type \code{\link{methQTLResult-class}} or a list of such objects
#' @param type The type of annotation to be overlapped. Needs to be \code{'SNP'}, \code{'CpG'} or \code{'cor.block'}
#' @return A list with two GRanges objects, one containing the overlapped set and the other the
#'    union of input data points as elements \code{'all.qtl'} and \code{'all.input'}
#' @author Michael Scherer
#' @export
get.overlap.universe <- function(meth.qtl.res,type){
  if(inherits(meth.qtl.res,"methQTLResult")){
    type.anno <- ifelse(type%in%c("CpG","cor.block"),"meth","geno")
    all.input <- getAnno(meth.qtl.res,type.anno)
    if(type%in%"SNP"){
      all.input$End <- as.numeric(as.character(all.input$Start))
    }
    if(type%in%"cor.block"){
      cor.blocks <- getCorrelationBlocks(meth.qtl.res)
      all.qtl <- paste(getResult(meth.qtl.res)[,"SNP"],paste(getResult(meth.qtl.res,cor.blocks)[,"CorrelationBlock"],collapse = "_"),sep="_")
    }else{
      all.qtl <- as.character(getResult(meth.qtl.res)[,type])
    }
  }
  if(is.list(meth.qtl.res)){
    if(!inherits(meth.qtl.res[[1]],"methQTLResult")){
      stop("Invalid value for meth.qtl.res, needs to be methQTLResult")
    }
    all.input <- overlap.inputs(meth.qtl.res,type=type)
    all.qtl <- overlap.QTLs(meth.qtl.res = meth.qtl.res,type=type)
  }
  if(type%in%"cor.block"){
    all.qtl <- unlist(sapply(all.qtl,function(qtl){
      strsplit(qtl,"_")
    }))
    all.qtl <- all.qtl[!grepl("rs",all.qtl)]
  }
  all.qtl <- all.input[all.qtl,]
  all.input <- makeGRangesFromDataFrame(all.input)
  all.qtl <- makeGRangesFromDataFrame(all.qtl)
  return(list("all.qtl"=all.qtl,"all.input"=all.input))
}
