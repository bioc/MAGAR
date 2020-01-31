##########################################################################################
# interpretation.R
# created: 2020-01-31
# creator: Michael Scherer
# ---------------------------------------------------------------------------------------
# functions to interpret methQTL results including LOLA and GO enrichments
##########################################################################################

#' qtl.lola.enrichment
#'
#' This functions performn LOLA enrichment analysis for the methQTL sites or the sites that are
#' shared across all methQTLs in the input list.
#'
#' @param meth.qtl.res An object of type \code{\link{methQTLResult-class}} or a list of such objects
#' @param type The type of methQTL to be visualized. Can be either \code{'SNP'}, \code{'CpG'},
#'     or \code{'cor.block'}
#' @param lola.db The location of a LOLA DB already downloaded
#' @return The LOLA enrichment result
#' @author Michael Scherer
#' @export
qtl.lola.enrichment <- function(meth.qtl.res,type="SNP",lola.db=NULL){
  if(!requireNamespace("LOLA")){
    stop("Please install the 'LOLA' R package")
  }
  stats <- get.overlap.universe(meth.qtl.res,type)
  all.input <- stats$all.input
  all.qtl <- stats$all.qtl
  if(is.null(lola.db)){
    assembly <- ifelse(inherits(meth.qtl.res,"methQTLResult"),meth.qtl.res@assmbly,meth.qtl.res[[1]]@assembly)
    lola.db <- downloadLolaDbs(tempdir())[[assembly]]
  }
  lola.db <- loadRegionDB(lola.db)
  lola.res <- runLOLA(all.input,all.qtl,lola.db)
  return(lola.res)
}

#' qtl.annotation.enrichment
#'
#' This functions performs enrichment analysis using the Fisher's test for the methQTLs detected
#' with respect to different genomic annotations.
#'
#' @param meth.qtl.res An object of type \code{\link{methQTLResult}} or a list of such objects.
#' @param type The type of methQTL to be visualized. Can be either \code{'SNP'}, \code{'CpG'},
#'     or \code{'cor.block'}
#' @param annotation The genomic annotation to be used. Can be the ones available in \code{\link{region.types()}}
#' @return The p-value of the Fisher exact test for enrichment
#' @details We use all data points that have been used to calculate methQTLs as the background
#'    and compare the overlaps with the annotation of interest in comparison to the methQTLs that
#'    have been computed.
#' @author Michael Scherer
#' @export
qtl.annotation.enrichment <- function(meth.qtl.res,
                                      type="SNP",
                                      annotation="cpgislands"){
  stats <- get.overlap.universe(meth.qtl.res,type)
  all.input <- stats$all.input
  all.qtl <- stats$all.qtl
  annotation <- unlist(rnb.get.annotation("cpgislands"))
  tps <- length(findOverlaps(all.qtl,annotation))
  fps <- length(all.qtl)-tps
  tns <- length(findOverlaps(all.input,annotation))
  fns <- length(all.input)-tns
  return(fisher.test(matrix(c(tps,fps,fns,tns),2,2),alternative="greater")$p.value)
}
