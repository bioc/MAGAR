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
           method="lm"
         ),
         package="methQTL")

# CONSTRUCTOR
setMethod("initialize","methQTLResult",
          function(.Object,
            result.frame=data.frame(),
            anno.meth=data.frame(),
            anno.geno=data.frame(),
            method="lm"
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
#'
#' Returns genomic annotation information for the given dataset.
#'
#' @param object An object of class \code{\link{methQTLResult-class}}.
#' @param type The type of annotation to be returned. Can either be \code{'meth'} or \code{'geno'} for methylation,
#'    and genotyping information, respectively.
#' @return The genomic annotation as a \code{data.frame}.
#'
#' @rdname getAnno
#' @docType methods
#' @aliases getAnno,methQTLResult-method
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
