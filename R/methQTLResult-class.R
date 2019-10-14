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
#'   \item{\code{method}}{The method used to call methQTL.}
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
           method="character"
         ),
         prototype(
           result.frame=data.frame(),
           anno.meth=data.frame(),
           anno.geno=data.frame(),
           method="classical.linear"
         ),
         package="methQTL")

# CONSTRUCTOR
setMethod("initialize","methQTLResult",
          function(.Object,
            result.frame=data.frame(),
            anno.meth=data.frame(),
            anno.geno=data.frame(),
            method="classical.linear"
          ){
            .Object@result.frame <- result.frame
            .Object@anno.meth <- anno.meth
            .Object@anno.geno <- anno.geno
            .Object@method <- method

            .Object
          })

##########################################################################################
# GETTERS
##########################################################################################
if(!isGeneric("getResult")) setGeneric("getResult",function(object) standardGeneric("getResult"))

#' getResult
#'
#' Returns the methQTL results stores in the object.
#'
#' @param object An of type \code{\link{methQTLResult-class}}.
#' @return The methQTL results as a \code{data.frame} with each row being a methQTL.
#' @rdname getResult
#' @docType methods
#' @aliases getResult,methQTLResult-method
#' @export
setMethod("getResult",signature(object="methQTLResult"),
          function(object){
            return(object@result.frame)
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

setMethod("show","methQTLResult",
          function(object){
            ret.str <- list()
            ret.str[1] <- "Object of class methQTLResult\n"
            ret.str[2] <- paste("\t Contains",nrow(object@result.frame),"methQTL\n")
            ret.str[3] <- paste("\t methQTL called using",object@method,"\n")
            cat(do.call(paste0,ret.str))
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
            object@result.frame <- data.frame()
            object@anno.meth <- data.frame()
            object@anno.geno <- data.frame()
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
         !file.exists(file.path(path,"anno_geno.RDS")))){
    stop("Invalid value for path. Potentially not a directory saved with save.methQTLResult")
  }
  load_env<-new.env(parent=emptyenv())
  load(file.path(path, "methQTLResult.RData"),envir=load_env)
  object <- get("object",load_env)
  result.frame <- readRDS(file.path(path,"result_frame.RDS"))
  anno.meth <- readRDS(file.path(path,"anno_meth.RDS"))
  anno.geno <- readRDS(file.path(path,"anno_geno.RDS"))
  object@result.frame <- result.frame
  object@anno.meth <- anno.meth
  object@anno.geno <- anno.geno
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
join.methQTL <- function(obj.list){
  if(any(!unlist(lapply(obj.list,function(x)inherits(x,"methQTLResult"))))){
    logger.error("Objects needs to be of type methQTLResult")
  }
  result.frame <- c()
  anno.meth <- c()
  anno.geno <- c()
  methods <- c()
  for(obj in obj.list){
    result.frame <- rbind(result.frame,getResult(obj))
    anno.meth <- rbind(anno.meth,getAnno(obj))
    anno.geno <- rbind(anno.geno,getAnno(obj,"geno"))
    methods <- c(methods,obj@method)
    if(any(methods != obj@method)){
      logger.error("Incompatible methQTL calling methods")
    }
  }
  result.frame <- data.frame(result.frame[order(result.frame[,1]),])
  anno.meth <- data.frame(anno.meth[order(anno.meth[,1]),])
  anno.geno <- data.frame(anno.geno[order(anno.geno[,1]),])
  ret.obj <- new("methQTLResult",
                 result.frame=result.frame,
                 anno.meth=anno.meth,
                 anno.geno=anno.geno,
                 method=methods[1])
  return(ret.obj)
}
