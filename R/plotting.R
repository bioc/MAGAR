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
#' @param meth.qtl.res An optional argument of type \code{\link{methQTLResult-class}} containing information on the results.
#'            If either \code{cpg} or \code{snp} are NULL, this function sorts the results by increasing p-value and the uses the
#'            best results for plotting.
#' @param out.name Optional name for the resulting plot
#' @return An object of type \code{ggplot} comparing the CpG methylation states as boxplots across the different genotype states
#' @author Michael Scherer
#' @export
qtl.plot.SNP.CpG.interaction <- function(meth.qtl,cpg=NULL,snp=NULL,out.dir=NULL,meth.qtl.res=NULL,out.name=NULL){
  if(!is.null(meth.qtl.res)){
    if(meth.qtl.res@rep.type == "mean.center"){
      logger.error("Interaction plot not available for representative CpG computation 'mean.center'. Pseudo CpG was created.")
    }
    res <- getResult(meth.qtl.res)
    res <- res[order(res$p.val.adj.fdr),]
    if(is.null(cpg)){
      cpg <- as.character(res$CpG[1])
    }
    if(is.null(snp)){
      snp <- as.character(res$SNP[1])
    }
  }
  sel.cpg <- which(row.names(getAnno(meth.qtl,"meth")) %in% cpg)
  sel.snp <- which(row.names(getAnno(meth.qtl,"geno")) %in% snp)
  sel.cpg <- getMethData(meth.qtl)[sel.cpg,]
  sel.snp <- getGeno(meth.qtl)[sel.snp,]
  if(is.integer(sel.snp)){
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
        if(is.null(out.name)){
          ggsave(file.path(out.dir,paste0("SNP_CpG_interaction_",cpg,"_",snp,".pdf")),plot)
        }else{
          ggsave(file.path(out.dir,out.name),plot)
        }
      }
    }
    return(plot)
  }else{
    to.plot <- data.frame(SNPDosage=sel.snp,CpG=sel.cpg)
    plot <- ggplot(to.plot,aes(x=SNP,y=CpG))+geom_point()+geom_smooth(method="lm",se=F)+theme_bw()+
      my_theme+
      ylab(paste(cpg,"methylation"))+xlab(paste(snp,"genotype"))
    if(!is.null(out.dir)){
      if(file.exists(out.dir)){
        if(is.null(out.name)){
          ggsave(file.path(out.dir,paste0("SNP_CpG_interaction_",cpg,"_",snp,".pdf")),plot)
        }else{
          ggsave(file.path(out.dir,out.name),plot)
        }
      }
    }
    return(plot)
  }
}

#' qtl.distance.scatterplot
#'
#' Computes a scatterplot between CpG-SNP distance with both effect size and p-value
#'
#' @param meth.qtl.result An object of type \code{\link{methQTLResult-class}} containing called methQTL
#' @param out.dir If specified, the plot is stored as a pdf in this directory
#' @param out.name Optional name for the resulting plot
#' @return An object of type \code{ggplot} comparing the distance between CpG and SNP. Negative values indicate that the
#'          SNP is downstream of the CpG.
#' @author Michael Scherer
#' @export
qtl.distance.scatterplot <- function(meth.qtl.result,out.dir=NULL,out.name=NULL){
  if(!requireNamespace("gridExtra")){
    stop("Please install the 'gridExtra' package")
  }
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
    g2 <- ggplot()+my_theme+annotate("text",y=0,x=0,label="No beta available")
    g3 <- ggplot()+my_theme+annotate("text",y=0,x=0,label="No beta available")
  }
  ret <- gridExtra::grid.arrange(g1,g2,g3,ncol=1)
  if(!is.null(out.dir)){
    if(file.exists(out.dir)){
      if(is.null(out.name)){
        ggsave(file.path(out.dir,"distance_correlations.pdf"),ret,width = 10,height = 12)
      }else{
        ggsave(file.path(out.dir,out.name),ret,width = 10,height = 12)
      }
    }
  }
}

#' qtl.manhattan.plot
#'
#' This function creates a manhattan plot for the given methQTL result
#'
#' @param meth.qtl.result An object of type \code{\link{methQTLResult-class}} containing the methQTL
#' @param type Determines if either the CpG (default) or the SNP is to be visualized
#' @param stat Determines the statistic that is to be visualized. Can be either \code{P.value}, \code{Beta} or
#'              \code{p.val.adj.fdr}
#' @details A plot is shown that contains chromosome-wise interactions.
#' @author Michael Scherer
#' @export
qtl.manhattan.plot <- function(meth.qtl.result,type="CpG",stat="p.val.adj.fdr"){
  if(!type %in% c("CpG","SNP")){
    stop("Invalid value for type, needs to be 'CpG' or 'SNP'")
  }
  if(!stat %in% c("P.value","Beta","p.val.adj.fdr")){
    stop("Invalid value for stat, needs to be 'P.value', 'Beta', or'p.val.adj.fdr'")
  }
  if(!requireNamespace("qqman")){
    stop("Please install the 'qqman' package")
  }
  to.plot <- getResult(meth.qtl.result)
  colnames(to.plot)[c(5,ifelse(type=="CpG",6,7),ifelse(stat=="P.value",4,ifelse(stat=="Beta",3,9)))] <- c("CHR","BP","P")
  to.plot$CHR <- as.numeric(gsub("chr","",to.plot$CHR))
  to.plot$P <- abs(to.plot$P)
  qqman::manhattan(to.plot,
            main=paste("Manhattan plot for",type,", statistic",stat),
            suggestiveline = F, genomewideline = F,
            ylab=ifelse(stat=="Beta","absolute beta",ifelse(stat=="P.value","-log10(P-value)","-log10(adjusted P-value)")),
            logp=stat!="Beta")
}

#' qtl.venn.plot
#'
#' This function creates a venn plot from a list of methQTL results, showing the overlap between the interactions
#'
#' @param meth.qtl.result.list A named list with each entry being an object of type \code{\link{methQTLResult-class}}.
#'                       The names are used in the visualization.
#' @param type Determines if either the CpG (default) or the SNP is to be visualized
#' @param out.folder Required argument specifying the location to store the resulting plot
#' @param out.name Optional argument for the name of the plot on disk (ending needs to be .png)
#' @param ... Further argument passed to \code{\link[VennDiagram]{venn.diagram}}
#' @return The venn plot object
#' @export
#' @author Michael Scherer
qtl.venn.plot <- function(meth.qtl.result.list,type="CpG",out.folder,out.name=NULL,...){
  if(length(meth.qtl.result.list)>4){
    stop("Venn plot only supports up to 4 results, consider using qtl.upset.plot")
  }
  if(!type %in% c("CpG","SNP")){
    stop("Invalid value for type, needs to be 'CpG' or 'SNP'")
  }
  if(!requireNamespace("VennDiagram")){
    stop("Please install the 'VennDiagram' package")
  }
  res.all <- lapply(meth.qtl.result.list,function(res){
    getResult(res)[,type]
  })
  names(res.all) <- names(meth.qtl.result.list)
  venn.colors <- rnb.getOption("colors.category")
  if(is.null(out.name)){
    out.file <- file.path(out.folder,"vennplot.tiff")
  }else{
    out.file <- file.path(out.folder,out.name)
  }
  ret <- VennDiagram::venn.diagram(res.all,
                      filename = out.file,
                      category.names = names(res.all),
                      fontfamily="sans",
                      fill=venn.colors[1:length(res.all)],
                      cat.cex=0.5,
                      ...)

  return(ret)
}

#' qtl.upset.plot
#'
#' This function creates an UpSet plot from the given methQTL results
#'
#' @param meth.qtl.result.list A named list with each entry being an object of type \code{\link{methQTLResult-class}}.
#'                       The names are used in the visualization.
#' @param type Determines if either the CpG (default) or the SNP is to be visualized
#' @param ... Further argument passed to \code{\link[UpSetR]{upset}}
#' @details The plot is directly drawn and can be stored on disk using the known R graphic devices
#' @export
#' @author Michael Scherer
qtl.upset.plot <- function(meth.qtl.result.list,type="CpG",...){
  if(!type %in% c("CpG","SNP")){
    stop("Invalid value for type, needs to be 'CpG' or 'SNP'")
  }
  if(!requireNamespace("UpSetR")){
    stop("Please install the 'UpSetR' package")
  }
  res.all <- lapply(meth.qtl.result.list,function(res){
    getResult(res)[,type]
  })
  names(res.all) <- names(meth.qtl.result.list)
  UpSetR::upset(UpSetR::fromList(res.all),
               nsets = length(res.all),
               order.by = "freq",
               mainbar.y.label = paste("Number of overlapping",type),
               sets.x.label = "methQTL per class",
               text.scale = c(1,1,1.5,1.5,2,1.5),
               number.angles = 30)
}
