suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(methQTL))

ap <- ArgumentParser()
ap$add_argument("-o","--output",action="store",help="Output directory")
cmd.args <- ap$parse_args()

logger.start("Combining results")
methQTL.result.files <- list.files(cmd.args$output,pattern = "methQTLResult_")
methQTL.results <- lapply(methQTL.result.files,load.methQTLResult)
methQTL.results <- join.methQTL(methQTL.results)
unlink(methQTL.result.files,recursive = T)
logger.completed()

logger.start("Saving results")
path.save <- file.path(cmd.args$output,paste0("methQTLResult"))
save.methQTLResult(methQTL.results,path.save)
logger.completed()
