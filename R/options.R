#' options.R
#'
#' This files contains code to generate the options of the methQTL package.
#'
#'
## G L O B A L S #######################################################################################################

QTL.OPTIONS <- new.env()

assign('ALL',c('rnbeads.options',
               'meth.data.type',
               'geno.data.type',
               'rnbeads.report',
               'rnbeads.qc',
               'hdf5dump',
               'hardy.weinberg.p',
             	 'db.snp.ref',
               'minor.allele.frequency',
               'missing.values.samples',
               'plink.geno',
               'plink.path',
		           'fast.qtl.path',
		           'bgzip.path',
		           'tabix.path',
		           'n.prin.comp',
		           'correlation.type',
		           'cluster.cor.threshold',
		           'standard.deviation.gauss',
		           'absolute.distance.cutoff',
               'linear.model.type',
		           'representative.cpg.computation',
		           'meth.qtl.type',
               'max.cpgs',
		           'rscript.path',
		           'cluster.config',
		           'recode.allele.frequencies',
               'n.permutations',
		           'p.value.correction',
		           'compute.cor.blocks',
               'use.segmentation',
		           'use.functional.annotation',
		           'functional.annotation.weight',
		           'impute.geno.data',
		           'vcftools.path',
		           'cluster.architecture',
		           'imputation.user.token',
		           'imputation.reference.panel',
		           'imputation.phasing.method',
		           'imputation.population'),QTL.OPTIONS)
assign('RNBEADS.OPTIONS',NULL,QTL.OPTIONS)
assign('METH.DATA.TYPE',"idat.dir",QTL.OPTIONS)
assign('GENO.DATA.TYPE',"plink",QTL.OPTIONS)
assign('RNBEADS.REPORT',"temp",QTL.OPTIONS)
assign('RNBEADS.QC',FALSE,QTL.OPTIONS)
assign('HDF5DUMP',FALSE,QTL.OPTIONS)
assign("HARDY.WEINBERG.P",0.001,QTL.OPTIONS)
assign("DB.SNP.REF",NULL,QTL.OPTIONS)
assign("MINOR.ALLELE.FREQUENCY",0.05,QTL.OPTIONS)
assign("MISSING.VALUES.SAMPLES",0.05,QTL.OPTIONS)
assign("PLINK.GENO",0.1,QTL.OPTIONS)
assign("N.PRIN.COMP",NULL,QTL.OPTIONS)
assign("CORRELATION.TYPE","pearson",QTL.OPTIONS)
assign("PLINK.PATH",NULL,QTL.OPTIONS)
assign("FAST.QTL.PATH",NULL,QTL.OPTIONS)
assign('BGZIP.PATH',NULL,QTL.OPTIONS)
assign("TABIX.PATH",NULL,QTL.OPTIONS)
assign("CLUSTER.COR.THRESHOLD",0.25,QTL.OPTIONS)
assign("STANDARD.DEVIATION.GAUSS",250,QTL.OPTIONS)
assign("ABSOLUTE.DISTANCE.CUTOFF",5e5,QTL.OPTIONS)
assign("LINEAR.MODEL.TYPE","classical.linear",QTL.OPTIONS)
assign("REPRESENTATIVE.CPG.COMPUTATION","row.medians",QTL.OPTIONS)
assign("METH.QTL.TYPE","oneVSall",QTL.OPTIONS)
assign("MAX.CPGS",40000,QTL.OPTIONS)
assign("RSCRIPT.PATH","/usr/bin/Rscript",QTL.OPTIONS)
assign("CLUSTER.CONFIG",c(h_vmem="5G",mem_free="5G"),QTL.OPTIONS)
assign("N.PERMUTATIONS",100,QTL.OPTIONS)
assign("P.VALUE.CORRECTION","uncorrected.fdr",QTL.OPTIONS)
assign("COMPUTE.COR.BLOCKS",TRUE,QTL.OPTIONS)
assign("RECODE.ALLELE.FREQUENCIES",TRUE,QTL.OPTIONS)
assign("USE.SEGMENTATION",FALSE,QTL.OPTIONS)
assign("USE.FUNCTIONAL.ANNOTATION",FALSE,QTL.OPTIONS)
assign("FUNCTIONAL.ANNOTATION.WEIGHT",1.1,QTL.OPTIONS)
assign("IMPUTE.GENO.DATA",FALSE,QTL.OPTIONS)
assign("VCFTOOLS.PATH",NULL,QTL.OPTIONS)
assign("IMPUTATION.USER.TOKEN",NULL,QTL.OPTIONS)
assign("IMPUTATION.REFERENCE.PANEL","apps@hrc-r1.1",QTL.OPTIONS)
assign("IMPUTATION.PHASING.METHOD","shapeit",QTL.OPTIONS)
assign("IMPUTATION.POPULATION","eur",QTL.OPTIONS)
assign("CLUSTER.ARCHITECTURE","sge",QTL.OPTIONS)

#' qtlSetOption
#'
#' Change global options for methQTL calculation
#'
#' @param rnbeads.options Path to an XML file specifying the RnBeads options used for data import. The default options
#'            are suitable for Illumina Beads Array data sets.
#' @param meth.data.type Type of DNA methylation data used. Choices are listed in \code{\link{rnb.execute.import}}.
#' @param geno.data.type The type of data to be imported. Can be either \code{'plink'} for \code{'.bed', '.bim',} and \code{'.fam'} or
#'   \code{'.dos'} and \code{'txt'} files or \code{'idat'} for raw IDAT files.
#' @param rnbeads.report Path to an existing directory, in which the preprocessing report of RnBeads is to be stored.
#'            Defaults to the temporary file.
#' @param rnbeads.qc Flag indicating if the quality control module of RnBeads is to be executed.
#' @param hdf5dump Flag indicating, if large matrices are to be stored on disk rather than in main memory using the
#'            \code{\link{HDF5Array}} package.
#' @param hardy.weinberg.p P-value used for the markers to be excluded if they do not follow the
#'            Hardy-Weinberg equilibrium as implemented in \code{PLINK}.
#' @param db.snp.ref Path to a locally stored version of dbSNP[3]. If this option is specified, the reference allele
#'             is determined from this file instead of from the allele frequencies of the dataset. This circumvents problems
#'	       with some imputation methods. If \code{NULL}(default), recoding will not be performed.
#' @param minor.allele.frequency Threshold for the minor allele frequency of the SNPs to be used in the analysis.
#' @param missing.values.samples Threshold specifying how much missing values per SNP are allowed across the samples
#'            to be included in the analyis.
#' @param plink.geno Threshold for missing values per SNP
#' @param impute.geno.data Flag indicating if imputation of genotyping data is to be perfomed using the Michigan imputation
#'            server (https://imputationserver.sph.umich.edu/index.html)[2].
#' @param n.prin.comp Number of principal components of the genetic data to be used as covariates
#'            in the methQTL calling. \code{NULL} means that no adjustment is conducted.
#' @param plink.path Path to an installation of PLINK (also comes with the package)
#' @param fast.qtl.path Path to an installation of fastQTL (comes with the package for Linux)
#' @param bgzip.path Path to an installation of BGZIP (comes with the package for Linux)
#' @param tabix.path Path to an installation of TABIX (comes with the package for Linux)
#' @param correlation.type The type of correlation to be used. Please note that for \code{type='pearson'} (default) the more efficient
#'          implementation of correlation in the \code{bigstatsr} is used. Further available options are \code{'spearman'} and
#'          \code{'kendall'}.
#' @param cluster.cor.threshold Threshold for CpG methylatin state correlation to be considered as connected in
#'            the distance graph used to compute the correlation clustering.
#' @param standard.deviation.gauss Standard deviation of the Gauss distribution used to weight the correlation
#'            according to its distance.
#' @param absolute.distance.cutoff Distance cutoff after which a CpG correlation is not considered anymore.
#' @param linear.model.type Linear model type to be used. Can be either \code{"categorical.anova"} or \code{"classical.linear"}. If \code{'meth.qtl.type'='fastQTL'}, this option is automatically set to \code{'fastQTL'}
#'            see \code{\link{callMethQTLBlock}} for more informations.
#' @param representative.cpg.computation Option specifying how reference CpGs per correlation block are to be computed. Available
#'            options are \code{"row.medians"} for the site that is the row median across the samples within the
#'            correlation block (for ties a random selection is performed), \code{"mean.center"} for an artifical site in the geometric center of the block with
#'            the average methylation level or \code{"best.all"} for the CpG with the best p-value across all of the
#'            CpGs in the correlation block.
#' @param meth.qtl.type Option specifying how a methQTL interaction is computed. Since the package is based on correlation
#'            blocks, a single correlation block can be associated with either one SNP (\code{meth.qtl.type='oneVSall'}),
#'            with multiple SNPs (\code{meth.qtl.type='allVSall'}), or each correlation block can once be positively and once
#'            negatively correlated with a SNP genotype (\code{meth.qtl.type='twoVSall'}). Additionally, we provide the option
#'            to use (\code{FastQTL})[1] as a methQTL mapping tool (option \code{'fastQTL'}).
#' @param max.cpgs Maximum number of CpGs used in the computation (used to save memory). 40,000 is a reasonable
#'             default for machines with ~128GB of main memory. Should be smaller for smaller machines and larger
#'             for larger ones.
#' @param cluster.architecture The type of HPC cluster architecture present. Currently supported are \code{'sge'} and \code{'slurm'}
#' @param cluster.config Resource parameters needed to setup an SGE or SLURM cluster job. Includes \code{h_vmem} and \code{mem_free} for SGE and \code{clock.limit} and \code{mem.size} for SLURM.
#' An example configuration for SLURM would be \code{c("clock.limit"="1-0","mem.size"="10G")} for 1 day of running time (format days:hours) and 10 GB of maximum memory usage. Additionally, \code{'n.cpus'} can be specified as the SLURM option \code{cpus-per-task}
#' @param rscript.path Path to an executable version of Rscript needed for submitting batch jobs to a cluster
#' @param n.permutations The number of permutations used to correct the p-values for multiple testing. See
#'              (http://fastqtl.sourceforge.net/) for further information.
#' @param p.value.correction The p-value correction method for multiple testing correction. Can be one of
#'              \code{"uncorrected.fdr"} or \code{"corrected.fdr}. \code{"uncorrected.fdr"} uses nominal p-values
#'              per correlation block as a lenient filtering and then uses FDR with the number of tests performed.
#'              \code{"corrected.fdr} corrects the p-values per correlation block for multiple testing, accounting
#'              for the correlation structure of the p-values and then uses FDR-correction for all the interactions computed.
#' @param compute.cor.blocks Flag indicating if correlation blocks are to be called. If \code{FALSE}, each CpG is considered
#'              separately.
#' @param use.segmentation Flag indicating if segmenation into partically methylated domains (PMDs) and non-PMDs is to
#'              be used for computing correlation block. If \code{TRUE}, two CpGs cannot be connected in the graph if they
#'              are in different segments.
#' @param use.functional.annotation Flag indicating if functional genome annotation according to the ENSEMBL regulatory
#'              build is to be incorportated into the correlation block computation. If \code{TRUE}, samples within the
#'              same annotation are prioritized and in different ones penalized.
#' @param functional.annotation.weight Numeric value specifying how strong CpGs in the same functional annotation should
#'              be prioritized.
#' @param recode.allele.frequencies Flag indicating if the reference allele is to be redefined according to the frequenciess
#'              found in the cohort investigated.
#' @param vcftools.path Path to the installation of VCFtools. Necessary is the vcf-sort function in this folder.
#' @param imputation.user.token The user token that is required for authorization with the Michigan imputation server. Please
#'              have a look at https://imputationserver.sph.umich.edu, create a user account and request a user token for access
#'              in your user profile.
#' @param imputation.reference.panel The reference panel used for imputation. Please see https://imputationserver.readthedocs.io/en/latest/reference-panels/
#'              for further information which panels are supported by the Michigan imputation server.
#' @param imputation.phasing.method The phasing method employed by the Michigan imputation server. See
#'              https://imputationserver.readthedocs.io/en/latest/api/ for further information.
#' @param imputation.population The population for the phasing method required by the Michigan imputation server. See
#'              https://imputationserver.readthedocs.io/en/latest/api/ for further information.
#' @return None
#' @export
#' @author Michael Scherer
#' @examples
#' {
#' qtlGetOption("rnbeads.report")
#' qtlSetOption(rnbeads.report=getwd())
#' qtlGetOption("rnbeads.report")
#' }
#' @references
#'   1. Ongen, H., Buil, A., Brown, A. A., Dermitzakis, E. T., & Delaneau, O. (2016).
#'    Fast and efficient QTL mapper for thousands of molecular phenotypes. Bioinformatics, 32(10),
#'    1479–1485. https://doi.org/10.1093/bioinformatics/btv722
#'   2. Das S, Forer L, Schönherr S, Sidore C, Locke AE, et al. (2016).
#'    Next-generation genotype imputation service and methods. Nature Genetics 48, 1284–1287,
#'    https://doi.org/10.1038/ng.3656
#'   3. Sherry, S. T. et al. (2001). dbSNP: the NCBI database of genetic variation.
#'    Nucleic Acids Res. 29, 308–311, https://doi.org/10.1093/nar/29.1.308.
qtlSetOption <- function(rnbeads.options=NULL,
                       meth.data.type="idat.dir",
                       geno.data.type="plink",
                       rnbeads.report="temp",
                       rnbeads.qc=FALSE,
                       hdf5dump=FALSE,
                       hardy.weinberg.p=0.001,
		                   db.snp.ref=NULL,
                       minor.allele.frequency=0.05,
                       missing.values.samples=0.05,
		                   plink.geno=0.1,
                       impute.geno.data=FALSE,
		                   n.prin.comp=NULL,
                       plink.path=NULL,
                       fast.qtl.path=NULL,
                       bgzip.path=NULL,
                       tabix.path=NULL,
		                   correlation.type="pearson",
                       cluster.cor.threshold=0.25,
                       standard.deviation.gauss=250,
                       absolute.distance.cutoff=5e5,
                       linear.model.type="classial.linear",
                       representative.cpg.computation="row.medians",
                       meth.qtl.type="oneVSall",
                       max.cpgs=40000,
                       rscript.path="/usr/bin/Rscript",
		                   cluster.architecture='sge',
                       cluster.config=c(h_vmem="5G",mem_free="5G"),
                       n.permutations=1000,
                       p.value.correction="uncorrected.fdr",
                       compute.cor.blocks=TRUE,
                       recode.allele.frequencies=FALSE,
                       use.segmentation=FALSE,
                       use.functional.annotation=FALSE,
                       functional.annotation.weight=1.1,
                       vcftools.path=NULL,
                       imputation.user.token=NULL,
                       imputation.reference.panel="apps@hrc-r1.1",
                       imputation.phasing.method="shapeit",
		                   imputation.population="eur"){
  if(length(rnbeads.options)>1 & !is.null(rnbeads.options)){
    stop("Please specify the options one by one, not as a vector or list.")
  }
  if(!missing(rnbeads.options)){
    if(is.null(rnbeads.options)){
      logger.info("Loading system default for option 'rnbeads.options'")
      rnbeads.options=system.file("extdata/rnbeads_options.xml",package="MAGAR")
    }
    if(!grepl(".xml",rnbeads.options)){
      stop("Invalid value for rnbeads.options: needs to be a path to a XML configuration file")
    }
    QTL.OPTIONS[['RNBEADS.OPTIONS']] <- rnbeads.options
  }
  if(!missing(meth.data.type)){
    if(!(meth.data.type %in% c("idat.dir",
                               "data.dir",
                               "data.files",
                               "GS.report",
                               "GEO",
                               "rnb.set"))){
      stop("Invalid value for meth.data.type, see rnb.execute.import for options.")
    }
    QTL.OPTIONS[['METH.DATA.TYPE']] <- meth.data.type
  }
  if(!missing(geno.data.type)){
    if(!(geno.data.type %in% c("idat","plink"))){
      stop("Invalid value for geno.data.type, needs to be 'idat' or 'plink'.")
    }
    QTL.OPTIONS[['GENO.DATA.TYPE']] <- geno.data.type
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
  if(!missing(hdf5dump)){
    if(!is.logical(hdf5dump)){
      stop("Invalid value for hdf5dump, needs to be TRUE/FALSE")
    }
    QTL.OPTIONS[['HDF5DUMP']] <- hdf5dump
  }
  if(!missing(hardy.weinberg.p)){
    if(!is.numeric(hardy.weinberg.p) && hardy.weinberg.p > 1){
      stop("Invalid value for hardy.weinberg.p, needs to be numeric < 1")
    }
    QTL.OPTIONS[['HARDY.WEINBERG.P']] <- hardy.weinberg.p
  }
  if(!missing(db.snp.ref)){
    if(!is.null(db.snp.ref) && !file.exists(db.snp.ref)){
      stop("Please download dbSNP from UCSC (https://genome.ucsc.edu/), and specify the path here")
    }
    QTL.OPTIONS[['DB.SNP.REF']] <- db.snp.ref
  }
  if(!missing(minor.allele.frequency)){
    if(!is.numeric(minor.allele.frequency) || minor.allele.frequency > 1){
      stop("Invalid value for minor.allele.frequency, needs to be numeric < 1")
    }
    QTL.OPTIONS[['MINOR.ALLELE.FREQUENCY']] <- minor.allele.frequency
  }
  if(!missing(missing.values.samples)){
    if(!is.numeric(missing.values.samples) || missing.values.samples > 1){
      stop("Invalid value for missing.values.samples, needs to be numeric < 1")
    }
    QTL.OPTIONS[['MISSING.VALUES.SAMPLES']] <- missing.values.samples
  }
  if(!missing(plink.geno)){
    if(!is.numeric(plink.geno) || plink.geno > 1){
      stop("Invalid value for plink.geno, needs to be numeric < 1")
    }
    QTL.OPTIONS[['PLINK.GENO']] <- plink.geno
  }
  if(!missing(impute.geno.data)){
    if(!is.logical(impute.geno.data)){
      stop("Invalid value for impute.geno.data, needs to be logical")
    }
    QTL.OPTIONS[['IMPUTE.GENO.DATA']] <- impute.geno.data
  }
  if(!missing(n.prin.comp)){
    if(!is.numeric(n.prin.comp) && !is.null(n.prin.comp)){
      stop("Invalid value for n.prin.comp, needs to be an integer or NULL")
    }
    QTL.OPTIONS[['N.PRIN.COMP']] <- n.prin.comp
  }
  if(!missing(plink.path)){
    if(!is.null(plink.path)){
      er <- tryCatch(system(plink.path,intern=TRUE),error=function(x)x)
      if(inherits(er,"error")){
        stop("Invalid value for plink.path, needs to be path to an executable")
      }
    }
    QTL.OPTIONS[['PLINK.PATH']] <- plink.path
  }
  if(!missing(fast.qtl.path)){
    if(!is.null(fast.qtl.path)){
      er <- tryCatch(system(fast.qtl.path,intern=TRUE),error=function(x)x)
      if(inherits(er,"error")){
        stop("Invalid value for fast.qtl.path, needs to be path to an executable")
      }
    }
    QTL.OPTIONS[['FAST.QTL.PATH']] <- fast.qtl.path
  }
  if(!missing(bgzip.path)){
    if(!is.null(bgzip.path)){
      er <- tryCatch(system(bgzip.path,timeout = 1, intern = TRUE),error=function(x)x)
      if(inherits(er,"error")){
        stop("Invalid value for bgzip.path, needs to be path to an executable")
      }
    }
    QTL.OPTIONS[['BGZIP.PATH']] <- bgzip.path
  }
  if(!missing(tabix.path)){
    if(!is.null(tabix.path)){
      er <- tryCatch(system(tabix.path,timeout = 1, intern = TRUE),error=function(x)x)
      if(inherits(er,"error")){
        stop("Invalid value for tabix.path, needs to be path to an executable")
      }
    }
    QTL.OPTIONS[['TABIX.PATH']] <- tabix.path
  }
  if(!missing(correlation.type)){
    if(!(correlation.type %in% c("pearson","spearman","kendall"))){
      stop("Invalid value for correlation.type, needs to be 'pearson', 'spearman', or 'kendall'")
    }
    QTL.OPTIONS[['CORRELATION.TYPE']] <- correlation.type
  }
  if(!missing(cluster.cor.threshold)){
    if(!is.numeric(cluster.cor.threshold) & cluster.cor.threshold >1){
      stop("Invalid value for cluster.cor.threshold, needs to be a numeric value between 0 and 1")
    }
    QTL.OPTIONS[['CLUSTER.COR.THRESHOLD']] <- cluster.cor.threshold
  }
  if(!missing(standard.deviation.gauss)){
    if(!is.numeric(standard.deviation.gauss)){
      stop("Invalid value for standard.deviation.gauss, needs to be numeric")
    }
    QTL.OPTIONS[['STANDARD.DEVIATION.GAUSS']] <- standard.deviation.gauss
  }
  if(!missing(absolute.distance.cutoff)){
    if(!is.numeric(absolute.distance.cutoff)){
      stop("Invalid value for absolute.distance.cutoff, needs to be numeric")
    }
    QTL.OPTIONS[['ABSOLUTE.DISTANCE.CUTOFF']] <- absolute.distance.cutoff
  }
  if(!missing(linear.model.type)){
    if(!linear.model.type %in% c("classical.linear",
                                 "categorical.anova",
                                 "fastQTL")){
      stop("Invalid value for linear.model.type. Needs to be classical linear or categorical.anova.")
    }
    QTL.OPTIONS[['LINEAR.MODEL.TYPE']] <- linear.model.type
  }
  if(!missing(representative.cpg.computation)){
    if(!representative.cpg.computation %in% c("row.medians",
                                              "mean.center",
                                              "best.all")){
      stop("Invalid value for representative.cpg.computation. Needs to be 'row.medians', 'mean.center' or 'best.all'.")
    }
    QTL.OPTIONS[['REPRESENTATIVE.CPG.COMPUTATION']] <- representative.cpg.computation
  }
  if(!missing(meth.qtl.type)){
    if(!meth.qtl.type%in%c("oneVSall",
                           "allVSall",
                           "twoVSall",
                           "fastQTL")){
      stop("Invalid value for meth.qtl.type. Needs to be 'oneVSall', 'allVSall', 'twoVSall', or 'fastQTL'")
    }
    QTL.OPTIONS[['METH.QTL.TYPE']] <- meth.qtl.type
    if(meth.qtl.type=="fastQTL"){
      QTL.OPTIONS[['LINEAR.MODEL.TYPE']] <- "fastQTL"
    }
  }
  if(!missing(max.cpgs)){
    if(!is.numeric(max.cpgs)){
      stop("Invalid value for max.cpgs. Needs to be numeric.")
    }
    QTL.OPTIONS[['MAX.CPGS']] <- max.cpgs
  }
  if(!missing(rscript.path)){
    if(!is.character(rscript.path)){
      stop("Invalid value for rscript.path, needs to be character")
    }
    tryCatch(o <- system(paste(rscript.path,"--version"),intern = TRUE),error=function(e){
      logger.warning("Invalid value for rscript.path. Needs to be a path to an executable version of Rscript")})
    QTL.OPTIONS[['RSCRIPT.PATH']] <- rscript.path
  }
  if(!missing(cluster.architecture)){
	if(!is.character(cluster.architecture) || !(cluster.architecture%in%c('sge','slurm'))){
		stop("Invalid value for cluster.architecture, needs to be 'sge' or 'slurm'")
	}
	QTL.OPTIONS[['CLUSTER.ARCHITECTURE']] <- cluster.architecture
  }
  if(!missing(cluster.config)){
    cluster.config <- unlist(cluster.config)
    if(!is.character(cluster.config)){
      stop("Invalid value for cluster.config, needs to be character")
    }else if(qtlGetOption("cluster.architecture")=="sge" && any(!(c("h_vmem","mem_free") %in% names(cluster.config)))){
		stop("h_vmem and mem_free required for cluster.architecture='sge'")
    }else if(qtlGetOption("cluster.architecture")=="slurm" && any(!(c("clock.limit","mem.size") %in% names(cluster.config)))){
		stop("clock.limit and mem.size required for cluster.architecture='slurm'")
    }
    QTL.OPTIONS[['CLUSTER.CONFIG']] <- cluster.config
  }
  if(!missing(n.permutations)){
    if(!is.numeric(n.permutations)){
      stop("Invalid value for n.permutations, needs to be numeric")
    }
    QTL.OPTIONS[['N.PERMUTATIONS']] <- n.permutations
  }
  if(!missing(p.value.correction)){
    if(!is.character(p.value.correction) || !p.value.correction %in% c("uncorrected.fdr","corrected.fdr")){
      stop("Invalid value for p.value.correction, needs to be either 'uncorrected.fdr' or 'corrected.fdr'")
    }
    QTL.OPTIONS[['P.VALUE.CORRECTION']] <- p.value.correction
  }
  if(!missing(compute.cor.blocks)){
    if(!is.logical(compute.cor.blocks)){
      stop("Invalid value for compute.cor.blocks, needs to be logical.")
    }
    QTL.OPTIONS[['COMPUTE.COR.BLOCKS']] <- compute.cor.blocks
  }
  if(!missing(recode.allele.frequencies)){
    if(!is.logical(compute.cor.blocks)){
      stop("Invalid value for recode.allele.frequencies, needs to be logical.")
    }
    QTL.OPTIONS[['RECODE.ALLELE.FREQUENCIES']] <- recode.allele.frequencies
  }
  if(!missing(use.segmentation)){
    if(!is.logical(use.segmentation)){
      stop("Invalid value for use.segmentation, needs to be logical.")
    }
    QTL.OPTIONS[['USE.SEGMENTATION']] <- use.segmentation
  }
  if(!missing(use.functional.annotation)){
    if(!is.logical(use.functional.annotation)){
      stop("Invalid value for use.functional.annotation, needs to be logical.")
    }
    QTL.OPTIONS[['USE.FUNCTIONAL.ANNOTATION']] <- use.functional.annotation
  }
  if(!missing(functional.annotation.weight)){
    if(!is.numeric(functional.annotation.weight)){
      stop("Invalid value for functional.annotation.weight, needs to be numeric.")
    }
    QTL.OPTIONS[['FUNCTIONAL.ANNOTATION.WEIGHT']] <- functional.annotation.weight
  }
  if(!missing(vcftools.path)){
    if(!is.null(vcftools.path)){
	    if(!is.character(vcftools.path) || !file.exists(file.path(vcftools.path,"vcf-sort"))){
	      stop("Invalid value for option 'vcftools.path', needs to point to a folder with the program 'vcf-sort'")
	    }
	    QTL.OPTIONS[['VCFTOOLS.PATH']] <- vcftools.path
    }
  }
  if(!missing(imputation.user.token)){
    if(!is.null(imputation.user.token)){
	    if(!is.character(imputation.user.token)){
	      stop("Invalid value for option 'imputation.user.token', needs to character")
	    }
	    QTL.OPTIONS[['IMPUTATION.USER.TOKEN']] <- imputation.user.token
    }
  }
  if(!missing(imputation.reference.panel)){
    if(!is.character(imputation.reference.panel)){
      stop("Invalid value for option 'imputation.reference.panel', needs to character")
    }
    QTL.OPTIONS[['IMPUTATION.REFERENCE.PANEL']] <- imputation.reference.panel
  }
  if(!missing(imputation.phasing.method)){
    if(!is.character(imputation.phasing.method)){
      stop("Invalid value for option 'imputation.phasing.method', needs to character")
    }
    QTL.OPTIONS[['IMPUTATION.PHASING.METHOD']] <- imputation.phasing.method
  }
  if(!missing(imputation.population)){
    if(!is.character(imputation.population)){
      stop("Invalid value for option 'imputation.population', needs to character")
    }
    QTL.OPTIONS[['IMPUTATION.POPULATION']] <- imputation.population
  }
}

#' qtlGetOption
#' Print the value of the global option
#'
#' @param names string or character vector containing the names of the options to be printed. All options are listed in \code{\link{qtlSetOption}}
#'
#' @return the option for the specified option
#' @author Michael Scherer
#' @export
#' @examples {
#' qtlGetOption("cluster.cor.threshold")
#' }
qtlGetOption <- function(names){
  if(!all(names %in% QTL.OPTIONS[['ALL']])){
    stop(paste0('No option(s) available named: ',names[!(names%in%QTL.OPTIONS[['ALL']])]))
  }
  ret <- c()
  if('rnbeads.options'%in%names){
    if(is.null(QTL.OPTIONS[['RNBEADS.OPTIONS']])){
      logger.info("Loading system default for option 'rnbeads.options'")
      qtlSetOption('rnbeads.options'=system.file("extdata/rnbeads_options.xml",package="MAGAR"))
    }
    ret <- c(ret,rnbeads.options=QTL.OPTIONS[['RNBEADS.OPTIONS']])
  }
  if('meth.data.type'%in%names){
    ret <- c(ret,meth.data.type=QTL.OPTIONS[['METH.DATA.TYPE']])
  }
  if('geno.data.type'%in%names){
    ret <- c(ret,geno.data.type=QTL.OPTIONS[['GENO.DATA.TYPE']])
  }
  if('rnbeads.report'%in%names){
    ret <- c(ret,rnbeads.report=QTL.OPTIONS[['RNBEADS.REPORT']])
  }
  if('rnbeads.qc'%in%names){
    ret <- c(ret,rnbeads.qc=QTL.OPTIONS[['RNBEADS.QC']])
  }
  if('hdf5dump'%in%names){
    ret <- c(ret,hdf5dump=QTL.OPTIONS[['HDF5DUMP']])
  }
  if('hardy.weinberg.p'%in%names){
    ret <- c(ret,hardy.weinberg.p=QTL.OPTIONS[['HARDY.WEINBERG.P']])
  }
  if('db.snp.ref'%in%names){
    ret <- c(ret,db.snp.ref=QTL.OPTIONS[['DB.SNP.REF']])
  }
  if('minor.allele.frequency'%in%names){
    ret <- c(ret,minor.allele.frequency=QTL.OPTIONS[['MINOR.ALLELE.FREQUENCY']])
  }
  if('missing.values.samples'%in%names){
    ret <- c(ret,missing.values.samples=QTL.OPTIONS[['MISSING.VALUES.SAMPLES']])
  }
  if('plink.geno'%in%names){
    ret <- c(ret,plink.geno=QTL.OPTIONS[['PLINK.GENO']])
  }
  if('impute.geno.data'%in%names){
    ret <- c(ret,impute.geno.data=QTL.OPTIONS[['IMPUTE.GENO.DATA']])
  }
  if('fast.qtl.path'%in%names){
    fast.qtl.path <- QTL.OPTIONS[['FAST.QTL.PATH']]
    if(is.null(fast.qtl.path)){
      logger.info("Loading system default for option 'fast.qtl.path'")
      fast.qtl.path=system.file("bin/fastQTL.static",package="MAGAR")
      er <- tryCatch(system(fast.qtl.path,timeout = 1, intern = TRUE),error=function(x)x)
      if(inherits(er,"error")){
        stop("Non-functional default version of fastQTL, please install it manually and specify it with 'fast.qtl.path'")
      }
    }
    ret <- c(ret,fast.qtl.path=fast.qtl.path)
  }
  if('n.prin.comp'%in%names){
    ret <- c(ret,n.prin.comp=QTL.OPTIONS[['N.PRIN.COMP']])
  }
  if('plink.path'%in%names){
    plink.path <- QTL.OPTIONS[['PLINK.PATH']]
    if(is.null(plink.path)){
      logger.info("Loading system default for option 'plink.path'")
      plink.path=system.file("bin/plink",package="MAGAR")
      er <- tryCatch(system(plink.path,timeout = 1, intern = TRUE),error=function(x)x)
      if(inherits(er,"error")){
        stop("Non-functional default version of plink, please install it manually and specify it with 'plink.path'")
      }
    }
    ret <- c(ret,plink.path=plink.path)
  }
  if('bgzip.path'%in%names){
    bgzip.path <- QTL.OPTIONS[['BGZIP.PATH']]
    if(is.null(bgzip.path)){
      logger.info("Loading system default for option 'bgzip.path'")
      tabix.path=system.file("bin/bgzip",package="MAGAR")
      er <- tryCatch(system(bgzip.path,timeout = 1, intern = TRUE),error=function(x)x)
      if(inherits(er,"error")){
        stop("Non-functional default version of bgzip, please install it manually and specify it with 'bgzip.path'")
      }
    }
    ret <- c(ret,bgzip.path=bgzip.path)
  }
  if('tabix.path'%in%names){
    tabix.path <- QTL.OPTIONS[['TABIX.PATH']]
    if(is.null(tabix.path)){
      logger.info("Loading system default for option 'tabix.path'")
      tabix.path=system.file("bin/tabix",package="MAGAR")
      er <- tryCatch(system(tabix.path,timeout = 1, intern = TRUE),error=function(x)x)
      if(inherits(er,"error")){
        stop("Non-functional default version of tabix, please install it manually and specify it with 'tabix.path'")
      }
    }
    ret <- c(ret,tabix.path=tabix.path)
  }
   if('correlation.type'%in%names){
    ret <- c(ret,correlation.type=QTL.OPTIONS[['CORRELATION.TYPE']])
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
  if('representative.cpg.computation'%in%names){
    ret <- c(ret,representative.cpg.computation=QTL.OPTIONS[['REPRESENTATIVE.CPG.COMPUTATION']])
  }
  if('meth.qtl.type'%in%names){
    ret <- c(ret,meth.qtl.type=QTL.OPTIONS[['METH.QTL.TYPE']])
  }
  if('max.cpgs'%in%names){
    ret <- c(ret,max.cpgs=QTL.OPTIONS[['MAX.CPGS']])
  }
  if('rscript.path'%in%names){
    ret <- c(ret,rscript.path=QTL.OPTIONS[['RSCRIPT.PATH']])
  }
  if('cluster.architecture'%in%names){
    ret <- c(ret,cluster.architecture=QTL.OPTIONS[['CLUSTER.ARCHITECTURE']])
  }
  if('cluster.config'%in%names){
    ret <- c(ret,cluster.config=list(QTL.OPTIONS[['CLUSTER.CONFIG']]))
  }
  if('n.permutations'%in%names){
    ret <- c(ret,n.permutations=QTL.OPTIONS[['N.PERMUTATIONS']])
  }
  if('p.value.correction'%in%names){
    ret <- c(ret,p.value.correction=QTL.OPTIONS[['P.VALUE.CORRECTION']])
  }
  if('compute.cor.blocks'%in%names){
    ret <- c(ret,compute.cor.blocks=QTL.OPTIONS[['COMPUTE.COR.BLOCKS']])
  }
  if('recode.allele.frequencies'%in%names){
    ret <- c(ret,recode.allele.frequencies=QTL.OPTIONS[['RECODE.ALLELE.FREQUENCIES']])
  }
  if('use.segmentation'%in%names){
    ret <- c(ret,use.segmentation=QTL.OPTIONS[['USE.SEGMENTATION']])
  }
  if('use.functional.annotation'%in%names){
    ret <- c(ret,use.functional.annotation=QTL.OPTIONS[['USE.FUNCTIONAL.ANNOTATION']])
  }
  if('functional.annotation.weight'%in%names){
    ret <- c(ret,functional.annotation.weight=QTL.OPTIONS[['FUNCTIONAL.ANNOTATION.WEIGHT']])
  }
  if('vcftools.path'%in%names){
    ret <- c(ret,vcftools.path=QTL.OPTIONS[['VCFTOOLS.PATH']])
  }
  if('imputation.user.token'%in%names){
    ret <- c(ret,imputation.user.token=QTL.OPTIONS[['IMPUTATION.USER.TOKEN']])
  }
  if('imputation.reference.panel'%in%names){
    ret <- c(ret,imputation.reference.panel=QTL.OPTIONS[['IMPUTATION.REFERENCE.PANEL']])
  }
  if('imputation.phasing.method'%in%names){
    ret <- c(ret,imputation.phasing.method=QTL.OPTIONS[['IMPUTATION.PHASING.METHOD']])
  }
  if('imputation.population'%in%names){
    ret <- c(ret,imputation.population=QTL.OPTIONS[['IMPUTATION.POPULATION']])
  }
  return(ret[names])
}

#' qtlOptions2JSON
#'
#' This function stores the current options setting as a JSON file at the specified path
#'
#' @param path A filename, to which the option setting is to be saved
#' @author Michael Scherer
#' @return None
#' @export
#' @examples {
#'   qtlSetOption('cluster.cor.threshold'=0.5)
#'   qtlOptions2JSON("my_opts.json")
#'   qtlJSON2options("my_opts.json")
#' }
#' @import jsonlite
qtlOptions2JSON <- function(path=file.path(getwd(),"methQTL_options.json")){
  all.options <- as.list(QTL.OPTIONS)
  all.options <- all.options[!(names(all.options) %in% "ALL")]
  names(all.options) <- sapply(names(all.options),tolower)
  all.options <- rjson::toJSON(all.options)#,null = "null")
  write(all.options,path)
}

#' qtlJSON2options
#'
#' This function reads an option setting from a JSON file and applies them to the current session
#'
#' @param path Path to a JSON file containing the options to be specified
#' @author Michael Scherer
#' @return None
#' @export
#' @import jsonlite
#' @examples {
#' qtlJSON2options(system.file("extdata/qtl_options_probesEPIC.json",package="MAGAR"))
#' }
qtlJSON2options <- function(path){
  if(!file.exists(path) || !grepl(".json",path,ignore.case = TRUE)){
    logger.error("Invalid value for path, needs to be a JSON file")
  }
  all.options <- fromJSON(path)
  all.options <- lapply(all.options,function(opt){
  	if(is(opt,"data.frame")){
             unlist(opt)
  	}else{
  	   opt
  	}
  })
  do.call(qtlSetOption,all.options)
}
