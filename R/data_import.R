##########################################################################################
# data_import.R
# created: 2019-08-28
# creator: Michael Scherer
# ---------------------------------------------------------------------------------------
# Methods to perform data import for methQTL analysis.
##########################################################################################

#' do.import
#'
#' Performs input for the given DNA methylation and genotyping data.
#'
#' @param data.location Named character vector specifying the data location. The names correspond to:
#'          \describe{
#'            \item{idat.dir}{Path to the idat-folder for raw DNA methylation data.}
#'            \item{plink.dir}{Path to the PLINK file directory (i.e. with \code{bed}, \code{bim} and \code{fam} files)}
#'          }
#' @param s.anno Path to the sample annotation sheet. If \code{NULL}, the program searches for potential sample
#'          annotation sheets in the data location directories.
#' @param assembly.meth Assembly used for the DNA methylation data. Typically is \code{"hg19"} for Illumina BeadArray
#'          data.
#' @param assembly.geno Assembly used for the genotyping data. If it is not the same as \code{assembly.geno}, the
#'          positions will be matched using liftOver.
#' @param tab.sep The table separator used for the sample annotation sheet.
#' @param s.id.col The column name of the sample annotation sheet that specifies the sample identifier.
#' @return An object of type \code{\link{methQTLInput-class}} with the methylation and genotyping information added.
#' @details Import of DNA methylation and genotyping data is done separately:
#'          \describe{
#'            \item{DNA methylation data}{DNA methylation data is imported using the \code{\link{RnBeads}} package. We
#'            use a default option setting commonly used for DNA methylation data obtained from the Illumina BeadArray
#'            series. If you want to specify further options, we refer to the \code{\link{rnb.options}}.}
#'            \item{Genotyping data}{Genotyping data is processed using PLINK. We focus on genotyping data generated
#'            with the Illumina BeadArray series and use default options. For further option settings, you should consult
#'            the \code{\link{qtl.setOption}} documentation.}
#'          }
#'
#'          If \code{assembly.meth} and \code{assembly.geno} do not match, we use the liftOver function to match the
#'          older assembly to the newer one and discard those sites, that are only covered by the newer version.
#'
#' @author Michael Scherer
#' @export

do.import <- function(data.location,s.anno=NULL,assembly.meth="hg19",assembly.geno="hg19",tab.sep=",",s.id.col="sample_id"){
  logger.start("Import methQTL data")
  if(is.null(s.anno)){
    for(i in 1:length(data.location)){
      all.files <- list.files(data.location,full.names=T)
      s.anno <- all.files[grepl("*annotation.csv",all.files)]
      if(is.null(s.anno)){
        s.anno <- all.files[grepl("*annotation.tsv",all.files)]
        tab.sep <- "\t"
        if(is.null(s.anno) && i==length(data.location)){
          stop("Could not find sample annotation sheet")
        }else if(!is.null(anno)){
          break
        }
      }else{
        break
      }
    }
  }
  pheno.data <- read.table(s.anno,sep=tab.sep,header = T)
  geno.import <- do.geno.import(data.location,pheno.data,s.id.col)
  pheno.data <- geno.import$pheno.data
  meth.import <- do.meth.import(data.location,assembly.meth,s.anno,s.id.col,tab.sep)
  s.names <- as.character(pheno.data[,s.id.col])
  if(is.null(s.names) || (length(unique(s.names)) < length(s.names))){
    stop("Invalid value for s.id.col, needs to specify unique identfiers in the sample annotation sheet")
  }
  if(s.names != meth.import$samples || s.names != geno.import$samples){
    stop("Samples are not present in all of the datasets")
  }
  row.names(pheno.data) <- s.names
  pheno.data <- pheno.data[,!(colnames(pheno.data) %in% s.id.col)]
  dataset.import <- new("methQTLInput",
    meth.data=meth.import$data,
    geno.data=geno.import$data,
    anno.meth=meth.import$annotation,
    anno.geno=geno.import$annotation,
    pheno.data=pheno.data,
    samples=s.names,
    assembly=assembly.meth,
    disk.dump=qtl.getOption("HDF5dump")
  )
  if(assembly.meth != assembly.geno){
    dataset.import <- match.assemblies(dataset.import)
  }
  logger.completed()
  return(dataset.import)
}

#' do.meth.import
#'
#' This function performs import of DNA methylation data and basic preprocessing steps implemented in the \code{\link{RnBeads}} package.
#'
#' @param data.location The location of the data files. See \code{\link{do.import}} for further details.
#' @param assembly The assembly used.
#' @param s.anno Path to the sample annotation sheet
#' @param s.id.col Column name in the sample annotation sheet specifying the sample identifiers.
#' @param tab.sep Table separator used.
#' @return A list with three elements:
#'          \describe{
#'            \item{sample}{The samples used in the dataset as a character vector}
#'            \item{data}{The DNA methylation data matrix as the results of import and preprocessing}
#'            \item{annotation}{The genomic annotation of the sites present in \code{data}}
#'          }
#' @details This function execute the import and preprocessing modules of RnBeads. In the default setting, a common
#'          option setting for Illumina BeadArray data is used (described in \code{inst/extdata/rnbeads_options.xml}).
#'          For further option setting, you should use the \code{rnbeads.options} option of the methQTL package to
#'          specify another XML file with custom RnBeads options.
#' @author Michael Scherer
#' @noRd
do.meth.import <- function(data.location,assembly="hg19",s.anno,s.id.col,tab.sep=","){
  require("RnBeads")
  logger.start("Processing DNA methylation data")
  rnb.xml2options(qtl.getOption("rnbeads.options"))
  rnb.options(identifiers.column=s.id.col,
              assembly=assembly,
              import.table.separator=tab.sep)
  data.s <- c(data.location["idat.dir"],s.anno)
  rnb.imp <- rnb.execute.import(data.source = data.s, data.type = qtl.getOption("meth.data.type"))
  if(qtl.getOption("rnbeads.qc")){
    if(dir.exists(qtl.getOption("rnbeads.report"))){
      rnb.report <- file.path(qtl.getOption("rnbeads.report"),"rnbeads_QC")
    }else{
      rnb.report <- file.path(tempdir(),"rnbeads_QC")
    }
    rnb.run.qc(rnb.imp,rnb.report)
  }
  if(dir.exists(qtl.getOption("rnbeads.report"))){
    rnb.report <- file.path(qtl.getOption("rnbeads.report"),"rnbeads_preprocessing")
  }else{
    rnb.report <- file.path(tempdir(),"rnbeads_preprocessing")
  }
  rnb.imp <- rnb.run.preprocessing(rnb.imp,rnb.report)$rnb.set
  s.names <- samples(rnb.imp)
  meth.data <- meth(rnb.imp)
  if(qtl.getOption("HDF5dump")){
    meth.data <- writeHDF5Array(meth.data)
  }
  anno.meth <- annotation(rnb.imp)
  logger.completed()
  return(list(samples=s.names,data=meth.data,annotation=anno.meth))
}

#' do.geno.import
#'
#' This routine performs data import and processing of genotyping data using the \code{PLINK} program.
#'
#' @param data.location Path to the directory, where the PLINK files are stored. For more information, see \code{\link{do.import}}
#' @param s.anno The sample annotation sheet.
#' @param s.id.col The column name of the sample annotation sheet specifying the sample identifiers.
#' @return A list of three elements:
#'        \describe{
#'          \item{data}{The processed genotyping data as a data.frame}
#'          \item{annotation}{The genomic annotation of the SNPs in \code{data}}
#'          \item{pheno.data}{The, potentially modified, sample annotation sheet}
#'        }
#' @author Michael Scherer
#' @noRd
do.geno.import <- function(data.location,s.anno,s.id.col){
  require(snpStats)
  logger.start("Processing genotyping data")
  snp.loc <- data.location["plink.dir"]
  all.files <- list.files(snp.loc,full.names=T)
  bed.file <- all.files[grepl(".bed",all.files)]
  bim.file <- all.files[grepl(".bim",all.files)]
  fam.file <- all.files[grepl(".fam",all.files)]
  if(any(c(is.null(bed.file),is.null(bim.file),is.null(fam.file)))){
    stop("Incompatible input to genotyping processing, needs to be in PLINK format (i.e. .bed, .bim and .fam files)")
  }
  snp.dat <- read.plink(bed = bed.file,bim = bim.file,fam = fam.file)
  fam <- snp.dat$fam
  s.anno <- s.anno[s.anno[,s.id.col] %in% row.names(fam),]
  if(is.null(s.anno)){
    stop("Sample ids and IDs of PLINK files do not match.")
  }
  keep.frame <- fam[as.character(s.anno[,s.id.col]),c("pedigree","member")]
  keep.file <- file.path(tempdir(),"kept_samples.txt")
  proc.data <- file.path(tempdir(),"processed_snp_data")
  write.table(keep.frame,keep.file,sep="\t",row.names=F,col.names=F,quote=F)
  plink.file <- gsub(".bed","",bed.file)
  cmd <- paste(qtl.getOption("plink.path"),"--bfile",plink.file,"--keep",keep.file,"--hwe",qtl.getOption("hardy.weinberg.p"),
               "--maf",qtl.getOption("minor.allele.frequency"),"--mind",qtl.getOption("missing.values.samples"),
               "--make-bed --out",proc.data)
  system(cmd)
  snp.dat <- read.plink(bed=paste0(proc.data,".bed"),bim=paste0(proc.data,".bim"),fam=paste0(proc.data,".fam"))
  snp.mat <- t(as(snp.dat$genotypes,"numeric"))
  snp.mat <- snp.mat[,s.anno[,s.id.col]]
  if(qtl.getOption("HDF5dump")){
    snp.mat <- writeHDF5Array(snp.mat)
  }
  anno.geno <- snp.dat$map
  colnames(anno.geno) <- c("Chromosome","Name","cM","Start","Allele.1","Allele.2")
  row.names(anno.geno) <- as.character(anno.geno$Name)
  anno.geno <- anno.geno[,c("Chromosome","Start","cM","Allele.1","Allele.2")]
  if(!any(grepl("chr*",anno.geno$Chromosome))){
    anno.geno$Chromosome <- paste0("chr",anno.geno$Chromosome)
  }
  logger.completed()
  return(list(data=snp.mat,annotation=anno.geno,pheno.data=s.anno,samples=s.anno[,s.id.col]))
}

match.assemblies <- function(meth.qtl){
  require(rtracklayer)
  print("Not yet implemented")
  return(meth.qtl)
}
