##########################################################################################
# plotting.R
# created: 2019-10-16
# creator: Michael Scherer
# ---------------------------------------------------------------------------------------
# plotting functions
##########################################################################################

my_theme <- theme_bw()+theme(panel.grid=element_blank(),text=element_text(size=18,color="black"),
                  axis.ticks=element_line(color="black"),plot.title = element_text(size=18,color="black",hjust = .5),
                  axis.text = element_text(size=15,color="black"))


#' qtl.plot.SNP.CpG.interaction
#'
#' Compares the methylation states of a given CpG for the genotype states availabe at the given SNP
#'
#' @param meth.qtl An object of type \code{\link{methQTLInput-class}} containing the methylation and genotype information
#'                  for the given CpG and the given SNP
#' @param cpg The CpG identifier as a character (e.g. cg12345678)
#' @param snp The SNP identifier as a character (e.g. rs12345678)
#' @param out.dir If specified, the plot is stored as a pdf in this directory
#' @return An object of type \code{ggplot} comparing the CpG methylation states as boxplots across the different genotype states
#' @author Michael Scherer
#' @export
qtl.plot.SNP.CpG.interaction <- function(meth.qtl,cpg,snp,out.dir=NULL){
  if(meth.qtl@rep.type == "mean.center"){
    logger.error("Interaction plot not available for representative CpG computation 'mean.center;. Pseudo CpG was created.")
  }
  sel.cpg <- which(row.names(getAnno(meth.qtl,"meth")) %in% cpg)
  sel.snp <- which(row.names(getAnno(meth.qtl,"geno")) %in% snp)
  sel.cpg <- getMethData(meth.qtl)[sel.cpg,]
  sel.snp <- getGeno(meth.qtl)[sel.snp,]
  count.geno <- c(0,0,0)
  for(gen in c(0,1,2)){
    count.geno[gen+1] <- sum(sel.snp == gen,na.rm=T)
  }
  to.plot <- data.frame(SNP=factor(ifelse(sel.snp==0,paste0("ref/ref (",count.geno[1],")"),ifelse(sel.snp=="1",paste0("ref/alt (",count.geno[2],")"),paste0("alt/alt (",count.geno[3],")"))),c(paste0("ref/ref (",count.geno[1],")"),paste0("ref/alt (",count.geno[2],")"),paste0("alt/alt (",count.geno[3],")"))),CpG=sel.cpg)
  plot <- ggplot(to.plot,aes(x=SNP,y=CpG))+geom_boxplot()+theme_bw()+
    my_theme+
    ylab(paste(cpg,"methylation"))+xlab(paste(snp,"genotype"))
  if(!is.null(out.dir)){
    if(file.exists(out.dir)){
      ggsave(file.path(out.dir,paste0("SNP_CpG_interaction_",cpg,"_",snp,".pdf")),plot)
    }
  }
  return(plot)
}

#' distance.scatterplot
#'
#' Computes a scatterplot between CpG-SNP distance with both effect size and p-value
#'
#' @param meth.qtl.result An object of type \code{\link{methQTLResult-class}} containing called methQTL
#' @param out.dir If specified, the plot is stored as a pdf in this directory
#' @return An object of type \code{ggplot} comparing the distance between CpG and SNP. Negative values indicate that the
#'          SNP is downstream of the CpG.
#' @author Michael Scherer
#' @export
distance.scatterplot <- function(meth.qtl.result,out.dir=NULL){
  require("gridExtra")
  res <- getResult(meth.qtl.result)
  cori.pval <- round(cor(abs(res$Distance),res$P.value),2)
  cori.pval.pval <- round(cor.test(abs(res$Distance),res$P.value)$p.value,3)
  g1 <- ggplot(res,aes(x=Distance,y=-log10(P.value)))+geom_point()+
    ggtitle(paste("Correlation btw absolute distance and p-value:",cori.pval,"p-value:",cori.pval.pval))+
    my_theme+xlab("Distance CpG - SNP")+ylab("-log10(methQTL p-value)")
  if(meth.qtl.result@rep.type != "mean.center"){
    cori.beta <- round(cor(abs(res$Distance),abs(res$Beta)),2)
    cori.beta.pval <- round(cor.test(abs(res$Distance),abs(res$Beta))$p.value,3)
    g2 <- ggplot(res,aes(x=Distance,y=abs(Beta),color=Beta))+geom_point()+
      ggtitle(paste("Correlation btw absolute distance and absolute beta:",cori.beta,"p-value:",cori.beta.pval))+
      my_theme+xlab("Distance CpG - SNP")+ylab("absolute beta")+labs(color="beta")+
      scale_color_gradient2(mid="#19547b",low = "#ffd89b",high = "chartreuse3")
    g3 <- ggplot(res,aes(x=Distance,y=Beta,color=-log10(P.value)))+geom_point()+
      my_theme+xlab("Distance CpG - SNP")+ylab("beta")+labs(color="-log10(p.val)")+
      scale_color_continuous(low="#19547b",high = "#ffd89b")
  }else{
    g2 <- ggplot()+my_theme+annotate("text",y=0,x=0,label="Not beta available")
    g3 <- ggplot()+my_theme+annotate("text",y=0,x=0,label="Not beta available")
  }
  ret <- grid.arrange(g1,g2,g3,ncol=1)
  if(!is.null(out.dir)){
    if(file.exists(out.dir)){
      ggsave(file.path(out.dir,"distance_correlations.pdf"),ret,width = 10,height = 12)
    }
  }
}
