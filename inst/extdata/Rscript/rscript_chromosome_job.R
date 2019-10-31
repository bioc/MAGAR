suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(methQTL))

ap <- ArgumentParser()
ap$add_argument("-m","--methQTL",action="store",help="methQTL object to be used")
ap$add_argument("-j","--json",action="store",help="Configuration JSON file")
ap$add_argument("-c","--chr",action="store",help="chromosome to be used")
ap$add_argument("-u","--covariates",action="store",default=NULL,help="covariates to be included")
ap$add_argument("-p","--p.val",action="store",default=1e-5,help="p-value cutoff")
ap$add_argument("-o","--output",action="store",help="Output directory")
cmd.args <- ap$parse_args()

logger.start(paste("Running on:",Sys.info()["nodename"]))

logger.start("Configuring job")
qtl.json2options(cmd.args$json)
logger.completed()

logger.start("Loading methQTL object")
meth.qtl <- load.methQTL(cmd.args$methQTL)
logger.completed()

if(!is.null(cmd.args$covariates)){
  logger.start("Reading covariates")
  covs <- unlist(readLines(cmd.args$covariates))
  logger.completed()
}else{
  covs <- NULL
}

p.val <- as.numeric(cmd.args$p.val)

methQTL.res <- do.methQTL.chromosome(meth.qtl,cmd.args$chr,sel.covariates = covs,p.val.cutoff = p.val, out.dir=cmd.args$output)

logger.start("Saving results")
path.save <- file.path(cmd.args$output,paste0("methQTLResult_",cmd.args$chr))
save.methQTLResult(methQTL.res,path.save)
logger.completed()
