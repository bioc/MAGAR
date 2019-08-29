##########################################################################################
# compute_methQTL.R
# created: 2019-08-29
# creator: Michael Scherer
# ---------------------------------------------------------------------------------------
# Methods to call methQTL from DNA methylation and genotyping data.
##########################################################################################

#' do.methQTL
#'
#' Function to compute methQTL given DNA methylation and genotyping data.
#'
#' @param meth.qtl An object of type \code{\link{methQTLInput-class}} on which methQTL computation is to be performed
#' @param sel.covariates Covariates as column names of the sample annotation sheet stored in \code{meth.qtl} to be
#'          used for covariate adjustment.
#' @param p.val.cutoff The p-value cutoff used for methQTL calling
#' @param ncores The number of cores used.
do.methQTL <- function(meth.qtl,sel.covariates=NULL,p.val.cutoff=1e-5,ncores=1){
  if(!inherits(meth.qtl,"methQTLInput")){
    stop("Invalid value fpr meth.qtl, needs to be of type methQTLInput")
  }
}
