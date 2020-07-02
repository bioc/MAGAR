##########################################################################################
# cluster.R
# created: 2019-10-01
# creator: Michael Scherer
# ---------------------------------------------------------------------------------------
# Methods to submit methQTL jobs to a SGE cluster.
##########################################################################################

#' submitClusterJobs
#'
#' This functions distributes chromosome-wise methQTL jobs across a SGE high performance computing cluster
#'
#' @param methQTL.input An object of type \code{\link{methQTLInput-class}} on which the methQTL are to be
#'         computed
#' @param covariates The selected covariates as a character vector
#' @param p.val.cutoff The p-value cutoff employed
#' @param out.dir The output directory
#' @param ncores The number of cores to be used
#' @return A list containing \code{\link{methQTLResult-class}} objects for each chromosome
#' @details The function was specifically created for a Sun Grid Engine (SGE) cluster, but can be extendend, to
#'  e.g. a SLURM architecture.
#' @author Michael Scherer
#' @noRd
submitClusterJobs <- function(methQTL.input,covariates,p.val.cutoff,out.dir,ncores=1){
  logger.start("Prepare cluster submission")
  json.path <- file.path(out.dir,"methQTL_configuration.json")
  qtl.options2json(json.path)
  if(!is.null(covariates)){
    cov.file <- file.path(out.dir,"methQTL_covariates.txt")
    writeLines(covariates,cov.file)
  }
  methQTL.file <- file.path(out.dir,"methQTLInput")
  if(!file.exists(methQTL.file)){
    save.methQTL(methQTL.input,methQTL.file)
  }
  logger.completed()
  all.chroms <- unique(getAnno(methQTL.input)$Chromosome)
  set.seed(Sys.time())
  id <- sample(1:10000,1)
  dep.tok <- ""
  req.res <- qtl.getOption("cluster.config")$cluster.config
  for(i in 1:length(names(req.res))){
    dep.tok <- paste(dep.tok,"-l",paste0(names(req.res)[i],"=",req.res[i]))
  }
  job.names <- sapply(all.chroms,function(chr){
    cmd.tok <- paste("qsub -V",
                     "-N",paste0("methQTL_",id,"_",chr),
                     "-o",file.path(out.dir,paste0("methQTL_",id,"_",chr,".log")),
                     dep.tok,
                     "-j y",
                     "-b y",
                     paste0("'",qtl.getOption("rscript.path")," ",system.file("extdata/Rscript/rscript_chromosome_job.R",package="methQTL")),
                     "-m",methQTL.file,
                     "-j",json.path,
                     "-c",chr,
                     "-p",p.val.cutoff,
                     "-o",paste0(out.dir,"'"),
                     "-n",ncores
                     )
    if(!is.null(covariates)){
      cmd.tok <- paste(cmd.tok,"-u",cov.file)
    }
    system(cmd.tok)
    paste0("methQTL_",id,"_",chr)
  })
  cmd.tok <- paste("qsub -V",
                   "-N",paste0("methQTL_",id,"_summary"),
                   "-o",file.path(out.dir,paste0("methQTL_",id,"_summary.log")),
                   dep.tok,
                   "-j y",
                   "-hold_jid",paste0(job.names,collapse = ","),
                   "-b y",
                   paste0("'",qtl.getOption("rscript.path")," ",system.file("extdata/Rscript/rscript_summary.R",package="methQTL")),
                   paste0("-o ",out.dir,"'")
                  )
  system(cmd.tok)
  logger.start("Waiting for jobs to finish")
  finished <- F
  while(!finished){
    Sys.sleep(100)
    qstat.res <- system(paste("qstat -j",paste0("methQTL_",id,"_summary")),ignore.stdout = T, ignore.stderr = T)
    # 0 for running, 1 for finished
    if(qstat.res == 1){
      finished <- T
    }
  }
  methQTL.result <- loadMethQTLResult(file.path(out.dir,"methQTLResult"))
  logger.completed()
  return(methQTL.result)
}
