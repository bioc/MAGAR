  #  test.import <- function(){
#    idat.files <- "/DEEP_fhgfs/projects/mscherer/data/EPIC/CEDAR/idats/"
#    s.anno <- "/DEEP_fhgfs/projects/mscherer/data/EPIC/CEDAR/annotation/sample_annotation_CD4_IPC.tsv"
#    plink.files <- "/DEEP_fhgfs/projects/mscherer/data/450K/CEDAR/publication/139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GENO//"
#    data.loc <- c(idat.dir=idat.files,geno.dir=plink.files)
#    qtl.setOption(hdf5dump=T)
#    res <- do.import(data.location = data.loc,s.anno = s.anno,s.id.col = "ind_IPC",tab.sep = "\t",
#                     out.folder="/TL/deep/projects/work/mscherer/projects/methQTL/test/")
#    checkTrue(inherits(res,"methQTLInput"))
#  }

#  test.import.imputed <- function(){
#    idat.files <- "/DEEP_fhgfs/projects/mscherer/data/EPIC/CEDAR/idats/"
#    s.anno <- "/DEEP_fhgfs/projects/mscherer/data/EPIC/CEDAR/annotation/sample_annotation_CD4_IPC.tsv"
#    imputed.files <- "/DEEP_fhgfs/projects/mscherer/data/450K/CEDAR/publication/139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GENO_IMPUTED/"
#    data.loc <- c(idat.dir=idat.files,geno.dir=imputed.files)
#    qtl.setOption(hdf5dump=T)
#    res <- do.import(data.location = data.loc,s.anno = s.anno,s.id.col = "ind_IPC",tab.sep = "\t",
#                     out.folder="/TL/deep/projects/work/mscherer/projects/methQTL/test/")
#    checkTrue(inherits(res,"methQTLInput"))
#  }

#  test.import.segmentation <- function(){
#    idat.files <- "/DEEP_fhgfs/projects/mscherer/data/EPIC/CEDAR/idats/"
#    s.anno <- "/DEEP_fhgfs/projects/mscherer/data/EPIC/CEDAR/annotation/sample_annotation_CD4_IPC.tsv"
#    plink.files <- "/DEEP_fhgfs/projects/mscherer/data/450K/CEDAR/publication/139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GENO//"
#    data.loc <- c(idat.dir=idat.files,geno.dir=plink.files)
#    qtl.setOption(hdf5dump=T)
#    qtl.setOption(use.segmentation=TRUE)
#    res <- do.import(data.location = data.loc,s.anno = s.anno,s.id.col = "ind_IPC",tab.sep = "\t",
#                     out.folder="/TL/deep/projects/work/mscherer/projects/methQTL/test/")
#    checkTrue(inherits(res,"methQTLInput"))
#  }

#  test.methQTL.segmentation <- function(){
#    idat.files <- "/DEEP_fhgfs/projects/mscherer/data/EPIC/CEDAR/idats/"
#    s.anno <- "/DEEP_fhgfs/projects/mscherer/data/EPIC/CEDAR/annotation/sample_annotation_CD4_IPC.tsv"
#    plink.files <- "/DEEP_fhgfs/projects/mscherer/data/450K/CEDAR/publication/139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GENO//"
#    data.loc <- c(idat.dir=idat.files,geno.dir=plink.files)
#    qtl.setOption(hdf5dump=T)
#    qtl.setOption(use.segmentation=TRUE)
#    res <- do.import(data.location = data.loc,s.anno = s.anno,s.id.col = "ind_IPC",tab.sep = "\t",
#                     out.folder="/TL/deep/projects/work/mscherer/projects/methQTL/test/")
#    meth.qtl <- do.methQTL.chromosome(res,
#	chrom="ch19",
#	sel.covariates=NULL,
#	p.val.cutoff=1e-5,
#	out.dir=getwd(),
#	ncores=10)
#    checkTrue(inherits(res,"methQTLResult"))
#  }

#   test.fastQTL <- function(){
#     meth.qtl <- load.methQTL("/DEEP_fhgfs/projects/mscherer/data/EPIC/CEDAR/methQTL/package/imputation/stringent/allVSall/RE/methQTLInput/")
#     qtl.setOption(meth.qtl.type="fastQTL",
#		bgzip.path="/TL/deep-share/archive00/software/packages/htslib/htslib-1.3.2/bgzip",
#		tabix.path="/TL/deep-share/archive00/software/packages/htslib/htslib-1.3.2/tabix")
#    meth.qtl <- do.methQTL.chromosome(meth.qtl,
#	chrom="chr20",
#	sel.covariates=NULL,
#	p.val.cutoff=1e-5,
#	out.dir=getwd(),
#	ncores=10)
#   }

test.idat <- function(){
 idat.files <- "/DEEP_fhgfs/projects/mscherer/data/450K/methQTLDo2016Tcells/idat/"
 s.anno <- "/DEEP_fhgfs/projects/mscherer/data/450K/methQTLDo2016Tcells/annotation/sample_annotation_genotypes_red.tsv"
 s.id.col <- "title"
 out.dir <- "/DEEP_fhgfs/projects/mscherer/data/450K/methQTLDo2016Tcells/test_package/"
 idat.platform="humanomni258v1p1b"
 gender.col=NULL
 qtl.setOption(meth.data.type="GEO")
 qtl.setOption(missing.values.samples = 0.5)
 data.loc <- c(idat.dir="GSE79144",geno.dir=idat.files)
 rnbeads.report = out.dir
 res <- do.import(data.loc,
                  data.type.geno="idat",
                  s.anno=s.anno,
                  tab.sep="\t",
                  s.id.col=s.id.col,
                  out.folder=out.dir,
                  idat.platform=idat.platform)
 res <- do.methQTL(res,
                   sel.covariates = NULL,
                   out.dir=out.dir,
                   ncores=10)
}
#  require("RUnit")
  require("methQTL")
#  test.fastQTL()
test.idat()

test.imputation <- function(){
  plink.files <- paste0("/DEEP_fhgfs/projects/mscherer/data/450K/methQTLDo2016Tcells/test_package/","processed_snp_data",c(".bed",".bim",".fam"))
  qtl.setOption("imputation.user.token"="eyJjdHkiOiJ0ZXh0XC9wbGFpbiIsImFsZyI6IkhTMjU2In0.eyJtYWlsIjoibXNjaGVyZXJAbXBpLWluZi5tcGcuZGUiLCJleHBpcmUiOjE1OTAxNTIzMjkwNzIsIm5hbWUiOiJNaWNoYWVsIFNjaGVyZXIiLCJhcGkiOnRydWUsInVzZXJuYW1lIjoibXNjaGVyZXIifQ.UJDTTrK-poRQbfRKhelvOB4QQQ4Xo13q7Q5g3Z1P-EU")
  qtl.setOption(vcftools.path="/TL/deep-share/archive00/software/packages/vcftools/vcftools/perl/")
  qtl.setOption(missing.values.samples = 0.75,
                bgzip.path = "/TL/deep-share/archive00/software/packages/htslib/htslib-1.3.2/bgzip",
                tabix.path = "/TL/deep-share/archive00/software/packages/htslib/htslib-1.3.2/tabix",
                imputation.reference.panel="apps@hrc-r1.1"
  )
  res <- do.imputation(plink.files[1],
                       plink.files[2],
                       plink.files[3],
                       out.dir="/DEEP_fhgfs/projects/mscherer/data/450K/methQTLDo2016Tcells/test_package/")
  bed.file <- plink.files[1]
  bim.file <- plink.files[2]
  fam.file <- plink.files[3]
  out.dir <- "/DEEP_fhgfs/projects/mscherer/data/450K/methQTLDo2016Tcells/test_package/"
}
library(methQTL)
test.imputation()
