#' options.R
#'
#' This files contains code to generate the options of the methQTL package.
#'
#'
## G L O B A L S #######################################################################################################

QTL.OPTIONS <- new.env()
assign('ALL',c('rnbeads.options','meth.data.type','rnbeads.report','rnbeads.qc','HDF5dump','hardy.weinberg.p',
               'minor.allele.frequency','missing.values.samples','plink.path',
               'cluster.cor.threshold','standard.deviation.gauss','absolute.distance.cutoff',
               'linear.model.type'),QTL.OPTIONS)
assign('RNBEADS.OPTIONS',system.file("extdata/rnbeads_options.xml",package="methQTL"),QTL.OPTIONS)
assign('METH.DATA.TYPE',"idat.dir",QTL.OPTIONS)
assign('RNBEADS.REPORT',"temp",QTL.OPTIONS)
assign('RNBEADS.QC',FALSE,QTL.OPTIONS)
assign('HDF5DUMP',FALSE,QTL.OPTIONS)
assign("HARDY.WEINBERG.P",0.001,QTL.OPTIONS)
assign("MINOR.ALLELE.FREQUENCY",0.05,QTL.OPTIONS)
assign("MISSING.VALUES.SAMPLES",0.05,QTL.OPTIONS)
assign("PLINK.PATH",system.file("bin/plink",package="methQTL"),QTL.OPTIONS)
assign("CLUSTER.COR.THRESHOLD",0.2,QTL.OPTIONS)
assign("STANDARD.DEVIATION.GAUSS",5000,QTL.OPTIONS)
assign("ABSOLUTE.DISTANCE.CUTOFF",1e6,QTL.OPTIONS)
assign("LINEAR.MODEL.TYPE","classical.linear",QTL.OPTIONS)

#' qtl.setOption
#'
#' Change global options for methQTL calculation
#'
#' @param rnbeads.options Path to an XML file specifying the RnBeads options used for data import. The default options
#'            are suitable for Illumina Beads Array data sets.
#' @param meth.data.type Type of DNA methylation data used. Choices are listed in \code{\link{rnb.execute.import}}.
#' @param rnbeads.report Path to an existing directory, in which the preprocessing report of RnBeads is to be stored.
#'            Defaults to the temporary file.
#' @param rnbeads.qc Flag indicating if the quality control module of RnBeads is to be executed.
#' @param HDF5dump Flag indicating, if large matrices are to be stored on disk rather than in main memory using the
#'            \code{\link{HDF5Array}} package.
#' @param hardy.weinberg.p P-value used for the markers to be excluded if the do not follow the
#'            Hardy-Weinberg equilibrium as implemented in \code{PLINK}.
#' @param minor.allele.frequency Threshold for the minor allele frequency of the SNPs to be used in the analysis.
#' @param missing.values.samples Threshold specifying how much missing values per SNP are allowed across the samples
#'            to be included in the analyis.
#' @param plink.path Path to an installation of PLINK (also comes with the package)
#' @param cluster.cor.threshold Threshold for CpG methylatin state correlation to be considered as connected in
#'            the distance graph used to compute the correlation clustering.
#' @param standard.deviation.gauss Standard deviation of the Gauss distribution used to weight the correlation
#'            according to its distance.
#' @param absolute.distance.cutoff Distance cutoff after which a CpG correlation is not considered anymore.
#' @param linear.model.type Linear model type to be used. Can be either \code{"categorical.anova"} or \code{"classical.linear"}.
#'            see \code{\link{call.methQTL.block}} for more informations.
#' @export
#' @author Michael Scherer
#' @examples
#' \donttest{
#' qtl.getOption("rnbeads.report")
#' qtl.setOption(rnbeads.report=getwd())
#' qtl.getOption("rnbeads.report")
#' }
qtl.setOption <- function(rnbeads.options=system.file("extdata/rnbeads_options.xml",package="methQTL"),
                       meth.data.type="idat.dir",
                       rnbeads.report="temp",
                       rnbeads.qc=F,
                       HDF5dump=F,
                       hardy.weinberg.p=0.001,
                       minor.allele.frequency=0.05,
                       missing.values.samples=0.05,
                       plink.path=system.file("bin/plink",package="methQTL"),
                       cluster.cor.threshold=0.2,
                       standard.deviation.gauss=1000,
                       absolute.distance.cutoff=1e6,
                       linear.model.type="classial.linear"){
  if(length(rnbeads.options)!=1){
    stop("Please specify the options one by one, not as a vector or list.")
  }
  if(!missing(rnbeads.options)){
    if(!grepl(".xml",rnbeads.options)){
      stop("Invalid value for rnbeads.options: needs to be a path to a XML configuration file")
    }
    QTL.OPTIONS[['RNBEADS.OPTIONS']] <- rnbeads.options
  }
  if(!missing(rnbeads.report)){
    if(!(meth.data.type %in% c("idat.dir","data.dir","data.files","GS.report","GEO","rnb.set"))){
      stop("Invalid value for meth.data.type, see rnb.execute.import for options.")
    }
    QTL.OPTIONS[['METH.DATA.TYPE']] <- meth.data.type
  }
  if(!missing(rnbeads.report)){
    if(!rnbeads.report =="temp" && !dir.exists(rnbeads.report)){
      stop("Invalid value for rnbeads.report, needs to be a path to an existing directory")
    }
    QTL.OPTIONS[['RNBEADS.REPORT']] <- rnbeads.report
  }
  if(!missing(rnbeads.qc)){
    if(!is.logical(rnbeads.qc)){
      stop("Invalid value for rnbeads.qc, needs to be TRUE/FALSE")
    }
    QTL.OPTIONS[['RNBEADS.QC']] <- rnbeads.qc
  }
  if(!missing(HDF5dump)){
    if(!is.logical(HDF5dump)){
      stop("Invalid value for HDF5dump, needs to be TRUE/FALSE")
    }
    QTL.OPTIONS[['HDF5DUMP']] <- HDF5dump
  }
  if(!missing(hardy.weinberg.p)){
    if(!is.numeric(hardy.weinberg.p) && hardy.weinberg.p > 1){
      stop("Invalid value for hardy.weinberg.p, needs to be numeric < 1")
    }
    QTL.OPTIONS[['HARDY.WEINBERG.P']] <- hardy.weinberg.p
  }
  if(!missing(minor.allele.frequency)){
    if(!is.numeric(minor.allele.frequency) && minor.allele.frequency > 1){
      stop("Invalid value for minor.allele.frequency, needs to be numeric < 1")
    }
    QTL.OPTIONS[['MINOR.ALLELE.FREQUENCY']] <- minor.allele.frequency
  }
  if(!missing(missing.values.samples)){
    if(!is.numeric(missing.values.samples) && missing.values.samples > 1){
      stop("Invalid value for missing.values.samples, needs to be numeric < 1")
    }
    QTL.OPTIONS[['MISSING.VALUES.SAMPLES']] <- missing.values.samples
  }
  if(!missing(plink.path)){
    er <- tryCatch(system(plink.path),function(x)x)
    if(inherits(er,"error")){
      stop("Invalid value for plink.path, needs to be path to an executable")
    }
    QTL.OPTIONS[['PLINK.PATH']] <- plink.path
  }
  if(!missing(cluster.cor.threshold)){
    if(!is.numeric(cluster.cor.threshold) & cluster.cor.threshold >1){
      stop("Invalid value for cluster.cor.threshold, needs to be a numeric value between 0 and 1")
    }
    QTL.OPTIONS['CLUSTER.COR.THRESHOLD'] <- cluster.cor.threshold
  }
  if(!missing(standard.deviation.gauss)){
    if(!is.numeric(standard.deviation.gauss)){
      stop("Invalid value for standard.deviation.gauss, needs to be numeric")
    }
    QTL.OPTIONS['STANDARD.DEVIATION.GAUSS'] <- standard.deviation.gauss
  }
  if(!missing(absolute.distance.cutoff)){
    if(!is.numeric(absolute.distance.cutoff)){
      stop("Invalid value for absolute.distance.cutoff, needs to be numeric")
    }
    QTL.OPTIONS['ABSOLUTE.DISTANCE.CUTOFF'] <- absolute.distance.cutoff
  }
  if(!missing(linear.model.type)){
    if(!linear.model.type %in% c("classical.linear","categorical.anova")){
      stop("Invalid value for linear.model.type. Needs to be classical linear or categorical.anova.")
    }
    QTL.OPTIONS['LINEAR.MODEL.TYPE'] <- linear.model.type
  }
}

#' qtl.getOption
#' Print the value of the global option
#'
#' @param names string or character vector containing the names of the options to be printed
#'
#' @return the option for the specified option
#' @author Michael Scherer
#' @export
qtl.getOption <- function(names){
  if(!all(names %in% QTL.OPTIONS[['ALL']])){
    stop(paste0('No option(s) available named: ',names[!(names%in%QTL.OPTIONS[['ALL']])]))
  }
  ret <- c()
  if('rnbeads.options'%in%names){
    ret <- c(ret,rnbeads.options=QTL.OPTIONS[['RNBEADS.OPTIONS']])
  }
  if('meth.data.type'%in%names){
    ret <- c(ret,meth.data.type=QTL.OPTIONS[['METH.DATA.TYPE']])
  }
  if('rnbeads.report'%in%names){
    ret <- c(ret,rnbeads.report=QTL.OPTIONS[['RNBEADS.REPORT']])
  }
  if('rnbeads.qc'%in%names){
    ret <- c(ret,rnbeads.qc=QTL.OPTIONS[['RNBEADS.QC']])
  }
  if('HDF5dump'%in%names){
    ret <- c(ret,HDF5dump=QTL.OPTIONS[['HDF5DUMP']])
  }
  if('hardy.weinberg.p'%in%names){
    ret <- c(ret,hardy.weinberg.p=QTL.OPTIONS[['HARDY.WEINBERG.P']])
  }
  if('minor.allele.frequency'%in%names){
    ret <- c(ret,minor.allele.frequency=QTL.OPTIONS[['MINOR.ALLELE.FREQUENCY']])
  }
  if('missing.values.samples'%in%names){
    ret <- c(ret,missing.values.samples=QTL.OPTIONS[['MISSING.VALUES.SAMPLES']])
  }
  if('plink.path'%in%names){
    ret <- c(ret,plink.path=QTL.OPTIONS[['PLINK.PATH']])
  }
  if('cluster.cor.threshold'%in%names){
    ret <- c(ret,cluster.cor.threshold=QTL.OPTIONS[['CLUSTER.COR.THRESHOLD']])
  }
  if('standard.deviation.gauss'%in%names){
    ret <- c(ret,standard.deviation.gauss=QTL.OPTIONS[['STANDARD.DEVIATION.GAUSS']])
  }
  if('absolute.distance.cutoff'%in%names){
    ret <- c(ret,absolute.distance.cutoff=QTL.OPTIONS[['ABSOLUTE.DISTANCE.CUTOFF']])
  }
  if('linear.model.type'%in%names){
    ret <- c(ret,linear.model.type=QTL.OPTIONS[['LINEAR.MODEL.TYPE']])
  }
  return(ret[names])
}
