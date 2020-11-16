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
#' @import argparse
#' @author Michael Scherer
#' @noRd
submitClusterJobs <- function(methQTL.input,covariates,p.val.cutoff,out.dir,ncores=1){
  logger.start("Prepare cluster submission")
  json.path <- file.path(out.dir,"methQTL_configuration.json")
  qtlOptions2JSON(json.path)
  if(!is.null(covariates)){
    cov.file <- file.path(out.dir,"methQTL_covariates.txt")
    writeLines(covariates,cov.file)
  }
  methQTL.file <- file.path(out.dir,"methQTLInput")
  if(!file.exists(methQTL.file)){
    saveMethQTL(methQTL.input,methQTL.file)
  }
  logger.completed()
  if(qtlGetOption("cluster.architecture")=='sge'){
    methQTL.result <- submitClusterJobsSGE(methQTL.input,out.dir,json.path,p.val.cutoff,ncores,covariates,cov.file,methQTL.file)
  }else if(qtlGetOption("cluster.architecture")=='slurm'){
submitClusterJobsSLURM(methQTL.input,out.dir,json.path,p.val.cutoff,ncores,covariates,cov.file,methQTL.file)
  }else{
	stop("You weren't supposed to be here...")
  }

  logger.completed()
  return(methQTL.result)
}

#' submitClusterJobsSGE
#' Implementation of submitClusterJobs for Sun Grid Engine
#' @noRd
submitClusterJobsSGE <- function(methQTL.input,
				out.dir,
				json.path,
				p.val.cutoff,
				ncores,
				covariates,
				cov.file,
				methQTL.file){
  all.chroms <- unique(getAnno(methQTL.input)$Chromosome)
  set.seed(Sys.time())
  id <- sample(1:10000,1)
  dep.tok <- ""
  req.res <- qtlGetOption("cluster.config")$cluster.config
  for(i in 1:length(names(req.res))){
    dep.tok <- paste(dep.tok,"-l",paste0(names(req.res)[i],"=",req.res[i]))
  }
  hdf.dir <- getHDF5DumpDir()
  ff.dir <- getOption("fftempdir")
  job.names <- sapply(all.chroms,function(chr){
    cmd.tok <- paste("qsub -V",
                     "-N",paste0("methQTL_",id,"_",chr),
                     "-o",file.path(out.dir,paste0("methQTL_",id,"_",chr,".log")),
                     dep.tok,
                     "-j y",
                     "-b y",
                     paste0("'",qtlGetOption("rscript.path")," ",system.file("extdata/Rscript/rscript_chromosome_job.R",package="MAGAR")),
                     "-m",methQTL.file,
                     "-j",json.path,
                     "-c",chr,
                     "-p",p.val.cutoff,
		     "-d",hdf.dir,
		     "-f",ff.dir,
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
                   paste0("'",qtlGetOption("rscript.path")," ",system.file("extdata/Rscript/rscript_summary.R",package="MAGAR")),
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
  return(methQTL.result)
}

#' submitClusterJobsSLURM
#' Implementation of submitClusterJobs for SLURM
#' @noRd
submitClusterJobsSLURM <- function(methQTL.input,
				out.dir,
				json.path,
				p.val.cutoff,
				ncores,
				covariates,
				cov.file,
				methQTL.file){
  all.chroms <- unique(getAnno(methQTL.input)$Chromosome)
  set.seed(Sys.time())
  id <- sample(1:10000,1)
  dep.tok <- ""
  req.res <- qtlGetOption("cluster.config")$cluster.config
  if(any(!(names(req.res)%in%c("clock.limit","mem.size","n.cpus")))){
	stop("Only 'clock.limit', 'mem.size', and 'n.cpus' currently supported for SLURM")
  }
  hdf.dir <- getHDF5DumpDir()
  ff.dir <- getOption("fftempdir")
  dep.tok <- paste(dep.tok,ifelse(is.null(req.res["clock.limit"]),"",paste("-t",req.res["clock.limit"])),
	ifelse(is.null(req.res["mem.size"]),"",paste0("--mem=",req.res["mem.size"])),
	ifelse(is.null(req.res["n.cpus"]),"",paste0("--cpus-per-task=",req.res["n.cpus"])))
  job.names <- sapply(all.chroms,function(chr){
    cmd.tok <- paste("sbatch --export=ALL",
                     paste0("--job-name=","methQTL_",id,"_",chr),
                     "-o",file.path(out.dir,paste0("methQTL_",id,"_",chr,".log")),
                     dep.tok,
                     paste0("--wrap='",qtlGetOption("rscript.path")," ",system.file("extdata/Rscript/rscript_chromosome_job.R",package="MAGAR")),
                     "-m",methQTL.file,
                     "-j",json.path,
                     "-c",chr,
                     "-p",p.val.cutoff,
                     "-n",ncores,
		     "-d",hdf.dir,
		     "-f",ff.dir,
                     "-o",paste0(out.dir,"'")
                     )
    if(!is.null(covariates)){
      cmd.tok <- paste(cmd.tok,"-u",cov.file)
    }
    as.numeric(gsub("Submitted batch job ","",system(cmd.tok,intern=T)))
    #paste0("methQTL_",id,"_",chr)
  })
  cmd.tok <- paste("sbatch --export=ALL",
                   paste0("--job-name=","methQTL_",id,"_summary"),
                   "-o",file.path(out.dir,paste0("methQTL_",id,"_summary.log")),
                   dep.tok,
                   paste0("--depend=",paste0(job.names,collapse = ",")),
                   paste0("--wrap='",qtlGetOption("rscript.path")," ",system.file("extdata/Rscript/rscript_summary.R",package="MAGAR")),
                   paste0("-o ",out.dir,"'")
                  )
  system(cmd.tok)
  logger.start("Waiting for jobs to finish")
  finished <- F
  while(!finished){
    Sys.sleep(100)
    qstat.res <- system(paste("squeue --name",paste0("methQTL_",id,"_summary")),ignore.stdout = T, ignore.stderr = T)
    # 0 for running, 1 for finished
    if(qstat.res == 1){
      finished <- T
    }
  }
  methQTL.result <- loadMethQTLResult(file.path(out.dir,"methQTLResult"))
  return(methQTL.result)
}
