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
#' @param meth.qtl.list If \code{type}=\code{'cor.block'}, a list of \code{\link{methQTLInput-class}} object
#'     in the same order as \code{meth.qtl.result.list} to obtain correlation block annotations from.
#' @return A list with \code{length(meth.qtl.result.list)} elements, containing IDs of methQTL
#'     interactions according to the option \code{type}.
#' @author Michael Scherer
#' @export
overlap.QTLs <- function(meth.qtl.result.list,type,meth.qtl.list=NULL){
  if(!type %in% c("CpG","SNP",'cor.block')){
    stop("Invalid value for type, needs to be 'CpG', 'SNP', or 'cor.block'")
  }
  if(type%in%"cor.block"){
    if(is.null(meth.qtl.list)){
      stop("methQTLInput objects need to be specify to obtain correlation blocks")
    }
    if(length(meth.qtl.list)!=length(meth.qtl.result.list)){
      stop("Sizes of meth.qtl.list and meth.qtl.result.list do not match")
    }
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
      cor.blocks <- getCorrelationBlocks(meth.qtl.result.list[[i]],meth.qtl.list[[i]])
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
