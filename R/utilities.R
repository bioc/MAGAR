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
    logger.start("Constructing universe")
    for(i in 1:length(meth.qtl.result.list)){
      logger.start(paste("Obtaining correlation blocks for object",i))
      cor.blocks <- getCorrelationBlocks(meth.qtl.result.list[[i]])
      res <- getResult(meth.qtl.result.list[[i]],cor.blocks)
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
          extended.qtl <- paste(snp,paste(cor.block[[1]],collapse = "_"),sep="_")
	  if(extended.qtl==res.all.vec[map.qtl[1]]){
		next
	  }else{
		old.cpgs <- unlist(strsplit(res.all.vec[map.qtl[1]],"_"))[-1]
		res.all.vec[map.qtl[[1]]] <-  paste(snp,paste(union(old.cpgs,cor.block[[1]]),collapse = "_"),sep="_")
	  }
        }else{
          #Caveat: only first interaction of a SNP with a correlation block will be considered, thus 'trans' effects
          #are not considered
          res.all.vec <- c(res.all.vec,paste(snp,paste(cor.block[[1]],collapse = "_"),sep="_"))
        }
      }
      logger.completed()
    }
    logger.completed()
    logger.start("Overlapping")
    for(i in 1:length(meth.qtl.result.list)){
      logger.start(paste("Obtaining correlation blocks for object",i))
      cor.blocks <- getCorrelationBlocks(meth.qtl.result.list[[i]])
      res <- getResult(meth.qtl.result.list[[i]],cor.blocks)
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
        }else{
	  # You weren't supposed to be here
	  stop("Constructing the universe failed")
       }
      }
      logger.completed()
      res.all[[i]] <- res.all.class
    }

    logger.completed()
  }
  names(res.all) <- names(meth.qtl.result.list)
  return(res.all)
}

#' overlap.inputs
#'
#' Overlaps the input annotations
#'
#' @param meth.qtl.list A list of \code{\link{methQTLInput-class}} or \code{\link{methQTLResult-class}} objects to be overlapped
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
  all.input <- getAnno(meth.qtl.list[[1]],type)
  for(i in 2:length(meth.qtl.list)){
    anno <- getAnno(meth.qtl.list[[i]],type)
    all.input <- rbind(all.input[,intersect(colnames(anno),colnames(all.input))],anno[!(row.names(anno)%in%row.names(all.input)),intersect(colnames(anno),colnames(all.input))])
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

#' get.overlap.universe
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
get.overlap.universe <- function(meth.qtl.res,type){
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
        paste(ro["SNP"],paste(ro["CorrelationBlock"][[1]],collapse = "_"),sep="_")
      })
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
    all.qtl <- all.qtl[!grepl("rs",all.qtl)]
  }
  if(type%in%"SNP"){
    all.input$End <- all.input$Start+1
  }
  all.qtl <- all.input[all.qtl,]
  all.input <- makeGRangesFromDataFrame(all.input)
  all.qtl <- makeGRangesFromDataFrame(all.qtl)
  return(list("all.qtl"=all.qtl,"all.input"=all.input))
}

#' get.specific.qtl
#'
#' This function returns the methQTL interactions specific for a result
#' @param meth.qtl.res An object of type \code{\link{methQTLResult-class}} for which specific QTLs are to be obtained.
#' @param meth.qtl.background The background set as a list of \code{\link{methQTLResult-class}} objects.
#' @param type The type of annotation to be overlapped. Needs to be \code{'SNP'}, \code{'CpG'} or \code{'cor.block'}
#' @return A \code{data.frame} of methQTL interactions sorted by the effect size.
#' @author Michael Scherer
#' @export
get.specific.qtl <- function(meth.qtl.res,meth.qtl.background,type="SNP"){
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
  op.qtl <- overlap.QTLs(c(meth.qtl.res,meth.qtl.background),type=type)
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
  return(res[order(abs(res$P.value),decreasing=F),])
}

#' get.overlapping.qtl
#'
#' This function merges the QTLs given and returns the methQTL table in a merged format.
#'
#' @param meth.qtl.list A list of \code{\link{methQTLResult-class}} objects to be merged
#' @param type The type of annotation to be overlapped. Needs to be \code{'SNP'}, \code{'CpG'} or \code{'cor.block'}
#' @return A \code{data.frame} with the methQTL interactions and an additional column specifying where the interaction
#'    displayed has been found. This value is generated from the \code{names()} argument of \code{meth.qtl.list}.
#' @author Michael Scherer
#' @export
get.overlapping.qtl <- function(meth.qtl.list,type="SNP"){
  op.qtl <- overlap.QTLs(meth.qtl.list,type=type)
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
