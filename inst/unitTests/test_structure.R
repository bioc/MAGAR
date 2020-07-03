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

execute.unit.tests <- function(){
  if(requireNamespace("RUnit")){
    logger.start("Unit Testing")
      logger.start("Testing constructors")
        test.constructors()
      logger.completed()
      logger.start("Testing options")
        test.options()
      logger.completed()
    logger.completed()
  }
}

execute.unit.tests()
