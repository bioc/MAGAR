#' Unit tests for methQTL R package

test.constructors <- function(){
  meth.matrix <- cbind(sample(1:10,5),sample(1:10,5),sample(1:10,5))
  geno.matrix <- cbind(sample(1:10,4),sample(1:10,4),sample(1:10,4))
  anno.meth <- as.data.frame(cbind(sample(1:10,5),sample(1:10,5)))
  anno.geno <- as.data.frame(cbind(sample(1:10,4),sample(1:10,4)))
  pheno.dat <- as.data.frame(rbind(c("A","B"),c("A","C"),c("B","D")))
  s.names <- paste0("S",1:3)
  colnames(meth.matrix) <- s.names
  colnames(geno.matrix) <- s.names
  row.names(pheno.dat) <- s.names
  obj <- new("methQTLInput",
             meth.data=meth.matrix,
             geno.data=geno.matrix,
             anno.meth=anno.meth,
             anno.geno=anno.geno,
             pheno.data=pheno.dat,
             samples=s.names)
  print(obj)
  res.frame <- as.data.frame(cbind(sample(1:10,10),sample(1:10,10)))
  obj <- new("methQTLResult",
             result.frame=res.frame,
             anno.meth=anno.meth,
             anno.geno=anno.geno)
  print(obj)
}

test.options <- function(){
  rnb.op <- "options.xml"
  qtlSetOption(rnbeads.options=rnb.op)
  res <- rnb.op == qtlGetOption("rnbeads.options")
  qtlSetOption(rnbeads.report=getwd())
  res <- res & qtlGetOption("rnbeads.report") == getwd()
  checkTrue(res)
}

test_calling <- function(){
    meth.qtl <- loadMethQTL(system.file("extdata","reduced_methQTL",package="MAGAR"))
    meth.qtl@meth.data <- as.matrix(meth.qtl@meth.data[1:10,])
    meth.qtl@geno.data <- as.matrix(meth.qtl@geno.data[1:10,])
    meth.qtl@anno.meth <- meth.qtl@anno.meth[1:10,]
    meth.qtl@anno.geno <- meth.qtl@anno.geno[1:10,]
    qtlSetOption(compute.cor.blocks=FALSE,
        p.value.correction="corrected.fdr",
        representative.cpg.computation="mean.center",
        rnbeads.options=system.file("extdata","rnbeads_options.xml",package="MAGAR"))
    meth.qtl.res <- doMethQTL(meth.qtl,p.val.cutoff=1)
    checkTrue(inherits(meth.qtl.res,"methQTLResult"))
}

test_cor_blocks <- function(){
    meth.qtl <- loadMethQTL(system.file("extdata","reduced_methQTL",package="MAGAR"))
    qtlSetOption(cluster.cor.threshold=0.1,
        standard.deviation.gauss=5000,
        absolute.distance.cutoff=100000)
    cor.blocks <- computeCorrelationBlocks(getMethData(meth.qtl),annotation=getAnno(meth.qtl))
    checkTrue(is(cor.blocks,"list"))
}

execute.unit.tests <- function(){
  if(requireNamespace("RUnit")){
    logger.start("Unit Testing")
      logger.start("Testing constructors")
        test.constructors()
      logger.completed()
      logger.start("Testing options")
        test.options()
      logger.completed()
      logger.start("Testing cor blocks")
        test_cor_blocks()
      logger.completed()
      logger.start("Test methQTL calling")
        test_calling()
      logger.completed()
    logger.completed()
  }
}

execute.unit.tests()
