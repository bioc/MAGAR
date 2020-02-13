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
#' @param assembly The assembly used (defaul \code{'hg19'})
#' @return The LOLA enrichment result
#' @details We use all data points that have been used to calculate methQTLs as the background
#'  and compare the overlaps with the annotation of interest in comparison to the methQTLs that
#'  have been computed in case a \code{\link{methQTLResult-class}} is provided. If a list of \code{\link{methQTLResult-class}} objects
#'  is provided, the intersection between the methQTLs from all objects in the list is compared with the union of all interactions
#'  that have been tested.
#' @author Michael Scherer
#' @export
qtl.lola.enrichment <- function(meth.qtl.res,type="SNP",lola.db=NULL,assembly="hg19"){
  if(!requireNamespace("LOLA")){
    stop("Please install the 'LOLA' R package")
  }else{
    library(LOLA)
  }
  stats <- get.overlap.universe(meth.qtl.res,type)
  all.input <- stats$all.input
  all.qtl <- stats$all.qtl
  if(is.null(lola.db)){
    lola.db <- downloadLolaDbs(tempdir())[[assembly]]
  }
  lola.db <- loadRegionDB(lola.db)
  lola.res <- runLOLA(userSets=all.qtl,userUniverse = all.qtl,lola.db)
  return(list(lola.res=lola.res,lola.db=lola.db))
}

#' qtl.annotation.enrichment
#'
#' This functions performs enrichment analysis using the Fisher's test for the methQTLs detected
#' with respect to different genomic annotations.
#'
#' @param meth.qtl.res An object of type \code{\link{methQTLResult-class}} or a list of such objects.
#' @param type The type of methQTL to be visualized. Can be either \code{'SNP'}, \code{'CpG'},
#'     or \code{'cor.block'}
#' @param annotation The genomic annotation to be used. Can be the ones available in \code{\link{region.types()}} or
#'     \code{"ctcf", "distal", "dnase", "proximal", "tfbs", "tss"}
#' @return A list of two p-values named \code{'enrichment'} for overrepresentation and \code{'depletion'} for underrepresentation
#' @details We use all data points that have been used to calculate methQTLs as the background
#'    and compare the overlaps with the annotation of interest in comparison to the methQTLs that
#'    have been computed in case a \code{\link{methQTLResult-class}} is provided. If a list of \code{\link{methQTLResult-class}} objects
#'    is provided, the intersection between the methQTLs from all objects in the list is compared with the union of all interactions
#'    that have been tested.
#' @author Michael Scherer
#' @export
qtl.annotation.enrichment <- function(meth.qtl.res,
                                      type="SNP",
                                      annotation="cpgislands"){
  ensembl.anno <- c("ctcf","distal","dnase","proximal","tfbs","tss")
  if(!(annotation%in%c("cpgislands","genes","promoters",ensembl.anno))){
    stop(paste("Invalid value for annotation, needs to be",c("cpgislands","genes","promoters",ensembl.anno)))
  }
  if(annotation%in%c("cpgislands","genes","promoters")){
    anno.type <- annotation
  }else{
    anno.type <- paste0("ensembleRegBuildBP",annotation)
    if(!anno.type%in%rnb.region.types()){
      rnb.load.annotation.from.db(anno.type)
    }
  }
  if(!(type%in%c("SNP","CpG","cor.block"))){
    stop("Invalif value for type, needs to be 'SNP','CpG', or 'cor.block'")
  }
  stats <- get.overlap.universe(meth.qtl.res,type)
  all.input <- stats$all.input
  all.qtl <- stats$all.qtl
  all.input <- all.input[!(names(all.input)%in%names(all.qtl))]
  annotation <- unlist(rnb.get.annotation(anno.type))
  tps <- length(findOverlaps(all.qtl,annotation))
  fps <- length(all.qtl)-tps
  fns <- length(findOverlaps(all.input,annotation))
  tns <- length(all.input)-fns
  gr <- fisher.test(matrix(c(tps,fns,fps,tns),2,2),alternative="greater")$p.value
  le <- fisher.test(matrix(c(tps,fns,fps,tns),2,2),alternative="less")$p.value
  or <- (tps/fps)/(fns/tns)
  return(list(enrichment=gr,depletion=le,OddsRatio=or))
}

#' qtl.base.substitution.enrichment
#'
#' This function tests for enrichment of a specific base substitution in the methQTL interactions.
#'
#' @param meth.qtl.res An object of type \code{\link{methQTLResult-class}} or a list of such objects.
#' @return A list with one element for each potential base substitution containing the enrichment p-value.
#' @details The names of the list are e.g. \code{A_G}, which refers to a replacement of the reference base \code{A}
#'   with an \code{A}. Enrichment is computed using Fisher's exact test, using all SNP that have been used
#'   as input as the background.
#' @author Michael Scherer
#' @export
qtl.base.substitution.enrichment <- function(meth.qtl.res){
  stats <- get.overlap.universe(meth.qtl.res,type="SNP")
  if(is.list(meth.qtl.res)){
    anno.cpgs <- overlap.inputs(meth.qtl.res,type="SNP")
  }else{
    anno.cpgs <- getAnno(meth.qtl.res,"geno")
  }
  sel.anno <- anno.cpgs[names(stats$all.qtl),]
  subs.qtl <- paste(sel.anno[,"Allele.1"],sel.anno[,"Allele.2"],sep="_")
  cu.qtl <- plyr::count(subs.qtl)
  sel.anno <- anno.cpgs[names(stats$all.input),]
  subs.input <- paste(sel.anno[,"Allele.1"],sel.anno[,"Allele.2"],sep="_")
  cu.input <- plyr::count(subs.input)
  all.interactions <- as.character(sapply(c("A","C","G","T"),function(b){
    paste(c("A","C","G","T")[!c("A","C","G","T")%in%b],b,sep ="_")
  }))
  enr.all <- sapply(all.interactions, function(int){
    tps <- cu.qtl[cu.qtl$x%in%int,"freq"]
    fps <- sum(cu.qtl[!(cu.qtl$x%in%int),"freq"])
    fns <- cu.input[cu.input$x%in%int,"freq"]
    tns <- sum(cu.input[!(cu.input$x%in%int),"freq"])
    p.val <- fisher.test(matrix(c(tps,fns,fps,tns),2,2),alternative="greater")$p.value
    or <- (tps/fps)/(fns/tns)
    return(c(p.value=p.val,OddsRatio=or))
  })
  return(enr.all)
}
