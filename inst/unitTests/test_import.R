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

#  require("RUnit")
#  require("methQTL")
#  test.methQTL.segmentation()

