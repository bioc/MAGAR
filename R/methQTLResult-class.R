##########################################################################################
# methQTLResult-class.R
# created: 2019-08-27
# creator: Michael Scherer
# ---------------------------------------------------------------------------------------
# methQTLResult class definition
##########################################################################################

#' methQTLResult-class
#'
#' Class storing methQTL analysis results and the associated genomic annotations
#'
#' @details
#'   This class stores the results of the methQTL analysis. It stores a \code{data.frame} with the methQTL results,
#'   and associated genomic annotations for both the methylation sites and SNPs.
#'
#' @section Slots:
#' \describe{
#'   \item{\code{result.frame}}{The methQTL results as a \code{data.frame}}
#'   \item{\code{anno.meth}}{Genomic annotation of the methylation sites as a \code{data.frame}.}
#'   \item{\code{anno.geno}}{Genomic annotation of the SNPs as a \code{data.frame}.}
#'   \item{\code{correlation.blocks}}{Correlation blocks determined from the methylation matrix.}
#'   \item{\code{method}}{The method used to call methQTL.}
#'   \item{\code{rep.type}}{Method used to determine representative CpGs from correlation blocks.}
#'   \item{\code{chr}}{Optional argument specifying if methQTL were called on a single chromosome.}
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{\link[=getResult,methQTLResult-method]{getResult}}}{Returns the methQTL results.}
#'   \item{\code{\link[=getAnno,methQTLResult-method]{getAnno}}}{Returns the genomic annotation.}
#' }
#'
#' @name methQTLResult-class
#' @rdname methQTLResult-class
#' @author Michael Scherer
#' @exportClass methQTLResult

setClass("methQTLResult",
         representation(
           result.frame="data.frame",
           anno.meth="data.frame",
           anno.geno="data.frame",
           correlation.blocks="list",
           method="character",
           rep.type="character",
           chr="characterOrNULL"
         ),
         prototype(
           result.frame=data.frame(),
           anno.meth=data.frame(),
           anno.geno=data.frame(),
           correlation.blocks=list(),
           method="classical.linear",
           rep.type="row.medians",
           chr=NULL
         ),
         package="methQTL")

# CONSTRUCTOR
setMethod("initialize","methQTLResult",
          function(.Object,
            result.frame=data.frame(),
            anno.meth=data.frame(),
            anno.geno=data.frame(),
            correlation.blocks=list(),
            method="classical.linear",
            rep.type="row.medians",
            chr=NULL
          ){
            .Object@result.frame <- result.frame
            .Object@anno.meth <- anno.meth
            .Object@anno.geno <- anno.geno
            .Object@correlation.blocks=correlation.blocks
            .Object@method <- method
            .Object@rep.type <- rep.type
            .Object@chr <- chr

            .Object
          })

##########################################################################################
# GETTERS
##########################################################################################
if(!isGeneric("getResult")) setGeneric("getResult",function(object,...) standardGeneric("getResult"))

#' getResult
#'
#' Returns the methQTL results stores in the object.
#'
#' @param object An of type \code{\link{methQTLResult-class}}.
#' @param cor.blocks Correlation blocks as obtained using \code{getCorrelationBlocks}. Please note that the
#'          correlation blocks need to contain the CpG identifiers, so the \code{\link{methQTLInput-class}} object
#'          needs to be provided to \code{getCorrelationBlocks}.
#' @return The methQTL results as a \code{data.frame} with each row being a methQTL.
#' @rdname getResult
#' @docType methods
#' @aliases getResult,methQTLResult-method
#' @export
setMethod("getResult",signature(object="methQTLResult"),
          function(object,cor.blocks=NULL){
            ret <- object@result.frame
            keep.lines <- apply(ret,1,function(line){
              any(!is.na(line))
            })
            ret <- ret[keep.lines,]
            if(!is.null(cor.blocks)){
              if(is.list(cor.blocks[[1]])){
                cor.blocks.assigned <- list()
                for(i in 1:length(cor.blocks)){
                  chr <- cor.blocks[[i]]
                  for(j in 1:length(chr)){
                    block <- chr[[j]]
                    cpg <- as.character(ret$CpG)[as.character(ret$CpG) %in% block]
                    if(length(cpg)>0){
                      cor.blocks.assigned[[cpg]] <- block
                    }
                  }
                }
              }else{
                cor.blocks.assigned <- list()
                for(i in 1:length(cor.blocks)){
                  block <- cor.blocks[[i]]
                  cpg <- as.character(ret$CpG)[as.character(ret$CpG) %in% block]
                  if(length(cpg)>0){
                    cor.blocks.assigned[[cpg]] <- block
                  }
                }
              }
              ret$CorrelationBlock <- cor.blocks.assigned[as.character(ret$CpG)]
            }
            return(ret)
          }
)

if(!isGeneric("getAnno")) setGeneric("getAnno",function(object,...) standardGeneric("getAnno"))

#' getAnno
#' @rdname getAnno
#' @docType methods
#' @aliases getAnno,methQTL-method
#' @export
setMethod("getAnno",signature(object="methQTLResult"),
          function(object,type="meth"){
            if(type=="meth"){
              return(object@anno.meth)
            }else if(type=="geno"){
              return(object@anno.geno)
            }else{
              stop("Invalid value for type: needs to be 'meth' or 'geno'")
            }
          }
)

if(!isGeneric("getCorrelationBlocks")) setGeneric("getCorrelationBlocks",function(object) standardGeneric("getCorrelationBlocks"))

#' getCorrelationBlocks
#'
#' Returns the correlation blocks defined for the given dataset
#'
#' @param object An object of class \code{\link{methQTLResult-class}}.
#' @return A \code{list} object containing the correlation blocks.
#' @rdname getCorrelationBlocks
#' @docType methods
#' @aliases getCorrelationBlocks,methQTLResult-method
#' @return
setMethod("getCorrelationBlocks",signature(object="methQTLResult"),
          function(object){
            cor.blocks <- object@correlation.blocks
            if(is.list(cor.blocks[[1]])){
              anno <- getAnno(object)
              ret <- list()
              if("chr1"%in%names(cor.blocks)){
                chr.ids <- names(cor.blocks)
              }else{
                chr.ids <- paste0("chr:",1:length(cor.blocks))
              }
              for(chr.id in chr.ids){
                chr <- cor.blocks[[chr.id]]
                chr.id <- ifelse(chr.id=="chr23","X",ifelse(chr.id=="chr24","Y",chr.id))
                anno.chr <- anno[anno$Chromosome%in%chr.id,]
                res <- lapply(chr,function(cg){
                  row.names(anno.chr)[cg]
                })
                ret[[chr.id]] <- res
              }
              ret <- ret[order(as.numeric(gsub("chr","",names(ret))))]
            }else{
              anno.chr <- getAnno(object)
              anno.chr <- anno.chr[anno.chr$Chromosome%in%object@chr,]
              ret <- lapply(cor.blocks,function(cg){
                row.names(anno.chr)[cg]
              })
            }
            return(ret)
          }
)

setMethod("show","methQTLResult",
          function(object){
            ret.str <- list()
            ret.str[1] <- "Object of class methQTLResult\n"
            ret.str[2] <- paste("\t Contains",nrow(getResult(object)),"methQTL\n")
            if(length(object@correlation.blocks)>0){
              if(is.list(object@correlation.blocks[[1]])){
                ret.str[3] <- paste("\t Contains",sum(lengths(object@correlation.blocks)),"correlation blocks\n")
              }else{
                ret.str[3] <- paste("\t Contains",length(object@correlation.blocks),"correlation blocks\n")
              }
            }else{
              ret.str[3] <- "\t Contains 0 correlation blocks\n"
            }
            ret.str[4] <- paste("\t methQTL called using",object@method,"\n")
            ret.str[5] <- paste("\t representative CpGs computed with",object@rep.type,"\n")
            if(!is.null(object@chr)){
              ret.str[6] <- paste("\t methQTL called for chromosome",object@chr,"\n")
            }
            cat(do.call(paste0,ret.str))
          }
)

if(!isGeneric("filter.pval")) setGeneric("filter.pval", function(object,...)standardGeneric("filter.pval"))

#' filter.pval
#'
#' This functions filters the methQTL results according to a given p-value cutoff
#'
#' @param object The \code{\link{methQTLResult-class}} object to be filtered
#' @param p.val.cutoff The p-value cutoff to be employed
#' @return The filtered \code{\link{methQTLResult-class}} object
#' @rdname filter.pval
#' @docType methods
#' @aliases filter.pval,methQTLResult-method
#' @author Michael Scherer
setMethod("filter.pval","methQTLResult",
          function(object,p.val.cutoff=0.01){
            res <- object@result.frame
            res <- res[res$p.val.adj.fdr <= p.val.cutoff,]
            object@result.frame <- res
            return(object)
          }
)

if(!isGeneric("save.methQTLResult")) setGeneric("save.methQTLResult", function(object,...)standardGeneric("save.methQTLResult"))

#' save.methQTLResult
#'
#' This functions stores a methQTLInputResult object in disk.
#'
#' @param object The \code{\link{methQTLResult-class}} object to be stored on disk.
#' @param path A path to a non-existing directory for files to be stored.
#'
#' @rdname save.methQTLResult
#' @docType methods
#' @aliases save.methQTLResult,methQTL-method
#' @author Michael Scherer
setMethod("save.methQTLResult","methQTLResult",
          function(object,path){
            if(file.exists(path)){
              if(dir.exists(path)){
                path <- file.path(path,"methQTLResult")
                if(file.exists(path)){
                  stop("Will not overwrite existing data")
                }
                dir.create(path)
              }else{
                stop("Will not overwrite existing data")
              }
            }else{
              dir.create(path)
            }
            saveRDS(object@result.frame,file=file.path(path,"result_frame.RDS"))
            saveRDS(object@anno.meth,file=file.path(path,"anno_meth.RDS"))
            saveRDS(object@anno.geno,file=file.path(path,"anno_geno.RDS"))
            saveRDS(object@correlation.blocks,file=file.path(path,"correlation_blocks.RDS"))
            object@result.frame <- data.frame()
            object@anno.meth <- data.frame()
            object@anno.geno <- data.frame()
            object@correlation.blocks <- list()
            save(object,file=file.path(path,"methQTLResult.RData"))
          }
)

#' load.methQTLResult
#'
#' This functions load a \code{\link{methQTLResult-class}} object from disk.
#'
#' @param path Path to the directory that has been created by \code{save.methQTLResult,methQTLInput-method}.
#' @return The object of type \code{\link{methQTLResult-class}} that has been stored on disk.
#' @author Michael Scherer
#' @export
load.methQTLResult <- function(path){
  if(any(!(file.exists(file.path(path,"result_frame.RDS"))),
         !file.exists(file.path(path,"anno_meth.RDS")),
         !file.exists(file.path(path,"anno_geno.RDS")),
         !file.exists(file.path(path,"correlation_blocks.RDS")))){
    stop("Invalid value for path. Potentially not a directory saved with save.methQTLResult")
  }
  load_env<-new.env(parent=emptyenv())
  load(file.path(path, "methQTLResult.RData"),envir=load_env)
  object <- get("object",load_env)
  result.frame <- readRDS(file.path(path,"result_frame.RDS"))
  anno.meth <- readRDS(file.path(path,"anno_meth.RDS"))
  anno.geno <- readRDS(file.path(path,"anno_geno.RDS"))
  correlation.blocks <- readRDS(file.path(path,"correlation_blocks.RDS"))
  object@result.frame <- result.frame
  object@anno.meth <- anno.meth
  object@anno.geno <- anno.geno
  object@correlation.blocks <- correlation.blocks
  return(object)
}

#' join.methQTLResult
#'
#' This function combines a list of \code{\link{methQTLResult-class}} objects.
#'
#' @param obj.list A list of \code{\link{methQTLResult-class}} objects to be joined
#' @return An object of type \code{\link{methQTLResult-class}} containing the combined information
#' @author Michael Scherer
#' @export
join.methQTLResult <- function(obj.list){
  if(any(!unlist(lapply(obj.list,function(x)inherits(x,"methQTLResult"))))){
    logger.error("Objects needs to be of type methQTLResult")
  }
  result.frame <- c()
  anno.meth <- c()
  anno.geno <- c()
  correlation.blocks <- list()
  methods <- c()
  rep.types <- c()
  for(obj in obj.list){
    result.frame <- rbind(result.frame,getResult(obj))
    anno.meth <- rbind(anno.meth,getAnno(obj))
    anno.geno <- rbind(anno.geno,getAnno(obj,"geno"))
    if(!is.null(obj@chr)){
      correlation.blocks[[obj@chr]] <- obj@correlation.blocks
    }else{
      correlation.blocks <- c(correlation.blocks,obj@correlation.blocks)
    }
    methods <- c(methods,obj@method)
    if(any(methods != obj@method)){
      logger.error("Incompatible methQTL calling methods")
    }
    rep.types <- c(rep.types,obj@rep.type)
    if(any(rep.types != obj@rep.type)){
      logger.error("Incompatible representative CpG computation methods")
    }
  }
  result.frame <- data.frame(result.frame[order(result.frame[,1]),])
  anno.meth <- data.frame(anno.meth[order(anno.meth[,1]),])
  anno.geno <- data.frame(anno.geno[order(anno.geno[,1]),])
  ret.obj <- new("methQTLResult",
                 result.frame=result.frame,
                 anno.meth=anno.meth,
                 anno.geno=anno.geno,
                 correlation.blocks=correlation.blocks,
                 method=methods[1],
                 rep.type=rep.types[1],
                 chr=NULL)
  return(ret.obj)
}
