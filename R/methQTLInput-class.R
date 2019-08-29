##########################################################################################
# methQTLInput-class.R
# created: 2019-08-27
# creator: Michael Scherer
# ---------------------------------------------------------------------------------------
# methQTLInput class definition
##########################################################################################

##########################################################################################
# CLASS DEFINITIONS
##########################################################################################
setClassUnion("matrixOrHDF",c("matrix","HDF5Matrix"))
setClassUnion("characterOrNULL",c("character","NULL"))

#' methQTLInput-class
#'
#' Class storing methQTL input data, such as DNA methylation and genotyping data, as well as sample metadata
#'
#' @details
#'   This class is the basis for computing methQTLs in the methQTL-package. It stores all the relevant information
#'   including methylation data and genotype data as a matrix or HDF5Matrix, the phenotypic data as a data frame
#'   and the genomic annotation of both the methylation sites and the SNP data.
#'
#' @section Slots:
#' \describe{
#'   \item{\code{meth.data}}{The methylation data as a numeric matrix of beta values or as an object of type \code{\link{HDF5Matrix}}}
#'   \item{\code{geno.data}}{The genotyping data as a numeric matrix of SNP genotypes (0=homozygote reference,
#'       1=heterozygote, 2=homozygote alternative allele) or as an object of type \code{\link{HDF5Matrix}}}
#'   \item{\code{pheno.data}}{Phenotypic data describing the samples used in the study. Matches the dimensions of
#'       both \code{meth.data} and \code{geno.data}}
#'   \item{\code{anno.meth}}{Genomic annotation of the methylation sites as a \code{data.frame}. Has the same number
#'       of rows as \code{meth.data}.}
#'   \item{\code{anno.geno}}{Genomic annotation of the SNPs as a \code{data.frame}. Has the same number of rows as
#'       \code{geno.data}.}
#'   \item{\code{samples}}{The sample identifiers used both for \code{meth.data} and \code{geno.data}, and as the rownames of
#'       \code{pheno.data}.}
#'   \item{\code{assembly}}{The genome assembly used.}
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{\link[=getMeth,methQTL-method]{getMeth}}}{Returns the methylation matrix.}
#'   \item{\code{\link[=getGeno,methQTL-method]{getGeno}}}{Returns the genotyping matrix.}
#'   \item{\code{\link[=getPheno,methQTL-method]{getPheno}}}{Returns the phenotypic information.}
#'   \item{\code{\link[=getAnno,methQTL-method]{getAnno}}}{Returns the genomic annotation.}
#' }
#'
#' @name methQTLInput-class
#' @rdname methQTLInput-class
#' @author Michael Scherer
#' @exportClass methQTLInput
setClass("methQTLInput",
         representation(
           meth.data="matrixOrHDF",
           geno.data="matrixOrHDF",
           pheno.data="data.frame",
           anno.meth="data.frame",
           anno.geno="data.frame",
           samples="characterOrNULL",
           assembly="character"
         ),
         prototype(
           meth.data=matrix(nrow=0,ncol=0),
           geno.data=matrix(nrow=0,ncol=0),
           pheno.data=data.frame(),
           anno.meth=data.frame(),
           anno.geno=data.frame(),
           samples=c(),
           assembly="hg19"
         ),
         package="methQTL")

# CONSTRUCTOR
setMethod("initialize","methQTLInput",
          function(.Object,
            meth.data=matrix(nrow=0,ncol=0),
            geno.data=matrix(nrow=0,ncol=0),
            pheno.data=data.frame(),
            anno.meth=data.frame(),
            anno.geno=data.frame(),
            samples=c(),
            assembly="hg19"
          ){
            if(length(samples) != ncol(meth.data) | length(samples) != ncol(geno.data) | length(samples) != nrow(pheno.data)){
              stop("Samples do not match dimension of the matrices.")
            }
            .Object@meth.data <- meth.data
            .Object@geno.data <- geno.data
            .Object@pheno.data <- pheno.data
            .Object@anno.meth <- anno.meth
            .Object@anno.geno <- anno.geno
            .Object@samples <- samples
            .Object@assembly <- assembly

            .Object
          })

#' get.value
#'
#' This functions returns the values for a particular sample and position from the matrix given in the first argument
#'
#' @param mat The matrix from which data is to be extracted. Can be either methylation of genotyping data.
#' @param site The sites to be selected either as a numeric or logical vector. If \code{NULL} all sites are returned.
#' @param sample The samples to be selected either as a numeric or logical vector. If \code{NULL} all samples are returned.
#' @return The selected values from the matrix either as a matrix or \code{HDF5Matrix}.
#' @author Michael Scherer
#' @noRd
get.value <- function(mat,site=NULL,sample=NULL){
  if(!is.element(class(site),c("NULL","integer","numeric","logical"))){
    stop("Invalid value for site, needs to be numeric, logical or NULL")
  }
  if(!is.element(class(sample),c("NULL","integer","numeric","logical"))){
    stop("Invalid value for site, needs to be numeric, logical or NULL")
  }
  if(is.null(site)){
    if(is.null(sample)){
      res <- mat
    }else{
      res <- mat[,sample]
    }
  }else{
    if(is.null(sample)){
      res <- mat[site,]
    }else{
      res <- mat[site,sample]
    }
  }
  return(res)
}
##########################################################################################
# GETTERS
##########################################################################################

if(!isGeneric("getMethData")) setGeneric("getMethData",function(object,...) standardGeneric("getMethData"))

#' getMethData
#'
#' Returns methylation information for the given dataset.
#'
#' @param object An object of class \code{\link{methQTLInput-class}}.
#' @param site The sites to be selected either as a numeric or logical vector. If \code{NULL} all sites are returned.
#' @param sample The samples to be selected either as a numeric or logical vector. If \code{NULL} all samples are returned.
#' @return The methylation matrix either as a matrix of \code{\link{HDF5Matrix}}.
#'
#' @rdname getMethData
#' @docType methods
#' @aliases getMethData,methQTL-method
#' @export
setMethod("getMethData",signature(object="methQTLInput"),
          function(object,site=NULL,sample=NULL){
              return(get.value(object@meth.data,site=site,sample=sample))

  }
)

if(!isGeneric("getGeno")) setGeneric("getGeno",function(object,...) standardGeneric("getGeno"))

#' getGeno
#'
#' Returns genotyping information for the given dataset.
#'
#' @param object An object of class \code{\link{methQTLInput-class}}.
#' @param site The sites to be selected either as a numeric or logical vector. If \code{NULL} all sites are returned.
#' @param sample The samples to be selected either as a numeric or logical vector. If \code{NULL} all samples are returned.
#' @return The genotyping matrix either as a matrix of \code{\link{HDF5Matrix}}.
#'
#' @rdname getGeno
#' @docType methods
#' @aliases getGeno,methQTL-method
#' @export
setMethod("getGeno",signature(object="methQTLInput"),
          function(object,site=NULL,sample=NULL){
            return(get.value(object@geno.data,site=site,sample=sample))
          }
)

if(!isGeneric("getPheno")) setGeneric("getPheno",function(object) standardGeneric("getPheno"))

#' getPheno
#'
#' Returns phenotypic information for the given dataset.
#'
#' @param object An object of class \code{\link{methQTLInput-class}}.
#' @return The phenotypic data either as a \code{data.frame}.
#'
#' @rdname getPheno
#' @docType methods
#' @aliases getPheno,methQTL-method
#' @export
setMethod("getGeno",signature(object="methQTLInput"),
          function(object){
            return(object@pheno.data)
          }
)

if(!isGeneric("getAnno")) setGeneric("getAnno",function(object,...) standardGeneric("getAnno"))

#' getAnno
#'
#' Returns genomic annotation information for the given dataset.
#'
#' @param object An object of class \code{\link{methQTLInput-class}}.
#' @param type The type of annotation to be returned. Can either be \code{'meth'} or \code{'geno'} for methylation,
#'    and genotyping information, respectively.
#' @return The genomic annotation as a \code{data.frame}.
#'
#' @rdname getAnno
#' @docType methods
#' @aliases getAnno,methQTL-method
#' @export
setMethod("getAnno",signature(object="methQTLInput"),
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

if(!isGeneric("getSamples")) setGeneric("getSamples",function(object) standardGeneric("getSamples"))

#' getSamples
#'
#' Returns the samples of the given dataset.
#'
#' @param object An object of class \code{\link{methQTLInput-class}}.
#' @return The samples of the dataset as a character vector.
#'
#' @rdname getSamples
#' @docType methods
#' @aliases getSamples,methQTL-method
#' @export
setMethod("getSamples",signature(object="methQTLInput"),
          function(object){
            return(object@samples)
          }
)

setMethod("show","methQTLInput",
  function(object){
    ret.str <- list()
    ret.str[1] <- "Object of class methQTLInput\n"
    ret.str[2] <- paste("\t Contains",length(object@samples),"samples\n")
    ret.str[3] <- paste("\t Methylation data for",nrow(object@meth.data),"CpGs\n")
    ret.str[4] <- paste("\t Genotyping data for",nrow(object@geno.data),"SNPs\n")
    ret.str[5] <- paste("\t Genome assembly:",object@assembly,"\n")
    cat(do.call(paste0,ret.str))
  }
)

if(!isGeneric("save.methQTL")) setGeneric("save.methQTL", function(object,...)standardGeneric("save.methQTL"))

setMethod("save.methQTL","methQTLInput",
          function(object,path){
            if(file.exists(path)){
              if(dir.exists(path)){
                path <- file.path(path,"methQTL")
              }else{
                stop("Will not overwrite existing data")
              }
            }
            save(object,file=file.path(path,"methQTLInput.RData"))
            saveRDS(object@meth.data,file=file.path(path,"meth_data.RDS"))
            saveRDS(object@geno.data,file=file.path(path,"meth_data.RDS"))
          }
)

load.methQTL <- function(path){

}
