test.import <- function(){
  idat.files <- "~/Documents/nfs/mscherer/data/EPIC/CEDAR/idats/"
  s.anno <- "~/Documents/nfs/mscherer/data/EPIC/CEDAR/annotation/sample_annotation_CD4_IPC.tsv"
  plink.files <- "~/Documents/nfs/mscherer/data/450K/CEDAR/publication/139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GENO//"
  data.loc <- c(idat.dir=idat.files,plink.dir=plink.files)
  res <- do.import(data.location = data.loc,s.anno = s.anno,s.id.col = "ind_IPC",tab.sep = "\t")
  checkTrue(inherits(res,"methQTLInput"))
}

require("RUnit")
require("methQTL")
test.import()

