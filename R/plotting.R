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
    plot <- ggplot(to.plot,aes(x=SNPDosage,y=CpG))+geom_point()+geom_smooth(method="lm",se=F)+theme_bw()+
      my_theme+
      ylab(paste(cpg,"methylation"))+xlab(paste(snp,"dosage"))
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

#' qtl.plot.SNP.correlation.block
#'
#' This functions creates a multi-facet plot with a panel for each CpG in the correlation block that has
#' a methQTL interaction with the SNP of interest.
#'
#' @param meth.qtl.res An object of type \code{\link{methQTLResult}} containing the methQTL interactions
#'           and the CpG correlation blocks.
#' @param meth.qtl An object of type \code{\link{methQTLInput}} containing genotyping, methylation and
#'           annotation information.
#' @param snp The SNP identifier, for which the methQTL interaction is to be visualized.
#' @return The plot as an \code{ggplot} object, containing a facet with a scatterplot between the CpG
#'    methylation state and the SNP dosage. Discrete genotypes are currently not supported.
#' @export
#' @author Michael Scherer
qtl.plot.SNP.correlation.block <- function(meth.qtl.res,meth.qtl,snp=NULL){
  if(!inherits(meth.qtl.res,"methQTLResult")){
    stop("Invalid value for meth.qtl.res, needs to be methQTLResult")
  }
  if(!inherits(meth.qtl,"methQTLInput")){
    stop("Invalid value for meth.qtl, needs to be methQTLInput")
  }
  cor.blocks <- getCorrelationBlocks(meth.qtl.res)
  res <- getResult(meth.qtl.res,cor.blocks)
  if(is.null(snp)){
    snp <- as.character(res[order(abs(res$Beta),decreasing = T),]$SNP[1])
  }
  sel.row <- res[res$SNP%in%snp,,drop=F]
  to.plot <- c()
  meth.data <- getMethData(meth.qtl)
  anno.meth <- getAnno(meth.qtl)
  anno.geno <- getAnno(meth.qtl,"geno")
  snp.data <- as.numeric(getGeno(meth.qtl)[row.names(anno.geno)%in%snp,,drop=F])
  for(i in 1:length(unlist(sel.row$CorrelationBlock))){
    cpg <- unlist(sel.row$CorrelationBlock)[i]
    add.mat <- cbind(SNP=snp.data,CpG=as.numeric(meth.data[row.names(anno.meth)%in%cpg,,drop=F]),ID=rep(cpg,length(snp.data)))
    to.plot <- rbind(to.plot,add.mat)
  }
  to.plot <- as.data.frame(to.plot)
  to.plot$CpG <- as.numeric(as.character(to.plot$CpG))
  to.plot$SNP <- as.numeric(as.character(to.plot$SNP))
  to.plot$Representative <- rep("No",nrow(to.plot))
  to.plot$Representative[to.plot$ID%in%sel.row$CpG] <- "Yes"
  plot <- ggplot(to.plot,aes(x=SNP,y=CpG))+geom_point(aes(color=Representative))+
    geom_smooth(method="lm",aes(color=Representative))+facet_grid(ID~.)+
    my_theme+theme(legend.position="none")+scale_color_manual(values=c("black","firebrick4"))
  return(plot)
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
#' @param out.folder Required argument specifying the location to store the resulting plot
#' @param type Determines if either the SNP (default), the CpG, or the correlation block
#'     \code{'cor.block'} is to be visualized
#' @param out.name Optional argument for the name of the plot on disk (ending needs to be .png)
#' @param ... Further argument passed to \code{\link[VennDiagram]{venn.diagram}}
#' @details The plot can be stored on disk using \code{out.folder} and \code{out.name}
#' @export
#' @author Michael Scherer
qtl.venn.plot <- function(meth.qtl.result.list,out.folder,type="SNP",out.name=NULL,
                          meth.qtl.list=NULL,...){
  if(length(meth.qtl.result.list)>4){
    stop("Venn plot only supports up to 4 results, consider using qtl.upset.plot")
  }
  if(!requireNamespace("VennDiagram")){
    stop("Please install the 'VennDiagram' package")
  }
  res.all <- overlap.QTLs(meth.qtl.result.list,type=type)
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
}

#' qtl.upset.plot
#'
#' This function creates an UpSet plot from the given methQTL results
#'
#' @param meth.qtl.result.list A named list with each entry being an object of type \code{\link{methQTLResult-class}}.
#'                       The names are used in the visualization.
#' @param type Determines if either the SNP (default), the CpG, or the correlation block
#'     \code{'cor.block'} is to be visualized
#' @param ... Further argument passed to \code{\link[UpSetR]{upset}}
#' @details The plot is directly drawn and can be stored on disk using the known R graphic devices
#' @export
#' @author Michael Scherer
qtl.upset.plot <- function(meth.qtl.result.list,type="SNP",...){
  if(!requireNamespace("UpSetR")){
    stop("Please install the 'UpSetR' package")
  }
  res.all <- overlap.QTLs(meth.qtl.result.list,type=type)
  UpSetR::upset(UpSetR::fromList(res.all),
               nsets = length(res.all),
               order.by = "freq",
               mainbar.y.label = paste("Number of overlapping",type),
               sets.x.label = "methQTL per class",
               text.scale = c(1,1,1.5,1.5,2,1.5),
               number.angles = 30)
}

#' qtl.correlate.cor.block.stat
#'
#' This function correlates the size of the correlation block, a particular CpG is part of, to the statistic
#' that has been found for a methQTL this CpG is involved in.
#'
#' @param meth.qtl.res An object of type \code{\link{methQTLResult-class}} containing methQTLs and correlation blocks.
#' @param stat The statistic that is to be compared. Can be one of \code{'P.value'}, \code{'Beta'},
#'          \code{'Distance'}, or \code{'p.val.adj.fdr'}(default).
#' @param size.type The size of correlation block to be used. Can be either \code{'num.CpGs'} (default)
#'          or \code{'genomic'}.
#' @return A scatterplot and associated correlations as an objec to type \code{ggplot}
#' @author Michael Scherer
#' @export
qtl.correlate.cor.block.stat <- function(meth.qtl.res,stat="p.val.adj.fdr",size.type='num.CpGs'){
  if(!inherits(meth.qtl.res,"methQTLResult")){
    stop("Invalid value for meth.qtl.res, needs to be methQTLResult")
  }
  if(!(stat %in% c("P.value","Beta","Distance","p.val.adj.fdr"))){
    stop("Invalid value for stat, needs to be one of 'P.value', 'Beta', 'Distance', 'p.val.adj.fdr'")
  }
  if(!(size.type%in%c("num.CpGs","genomic"))){
    stop("Invalid value for size.type, needs to be 'num.CpGs' or 'genomic'")
  }
  logger.start("Retriving correlation blocks")
  cor.blocks <- getCorrelationBlocks(meth.qtl.res)
  logger.completed()
  logger.start("Retrieving methQTL information")
  res.frame <- getResult(meth.qtl.res,cor.blocks)
  if(size.type%in%"num.CpGs"){
    res.frame$CorBlockSize <- sapply(res.frame$CorrelationBlock,length)
  }else{
    anno.meth <- getAnno(meth.qtl.res)
    res.frame$CorBlockSize <- sapply(res.frame$CorrelationBlock,function(block){
      if(length(block)==1){
        return(1)
      }else{
        start <- min(anno.meth[block,"Start"])
        end <- max(anno.meth[block,"Start"])
        return(end-start)
      }
    })
  }
  logger.completed()
  if(stat %in% c("P.value","p.val.adj.fdr")){
    val <- -log10(res.frame[,stat])
  }else if(stat%in%"Beta"){
    val <- abs(res.frame[,stat])
  }else{
    val <- res.frame[,stat]
  }
  cori <- cor(res.frame$CorBlockSize,val)
  cori.pval <- cor.test(res.frame$CorBlockSize,val)$p.value
  res.frame$val <- val
  plot <- ggplot(res.frame,aes_string(x="CorBlockSize",y=val))+geom_point()+geom_smooth(method="lm")+my_theme+
    ggtitle(paste("Correlation btw correlation block size and",stat,round(cori,3),"p-value:",round(cori.pval,3)))+
    ylab(paste0(ifelse(stat%in%c("P.value","p.val.adj.fdr"),"-log10(",ifelse(stat%in%"Beta","abs(","")),stat,
                ifelse(stat%in%c("P.value","p.val.adj.fdr","Beta"),")","")))
  return(plot)
}

#' qtl.lola.plot
#'
#' This function plots the LOLA enrichment results using the \code{\link{lolaBarPlot}} routine.
#'
#' @param meth.qtl.res An object of type \code{\link{methQTLResult-class}} or a list of such objects
#' @param type The type of methQTL to be visualized. Can be either \code{'SNP'}, \code{'CpG'},
#'     or \code{'cor.block'}
#' @param lola.db The location of a LOLA DB already downloaded
#' @param assembly The assembly used (defaul \code{'hg19'})
#' @param pvalCut The p-value cutoff employed
#' @return The LOLA enrichment bar plot as a \code{ggplot} object
#' @author Michael Scherer
#' @export
qtl.lola.plot <- function(meth.qtl.res,type="SNP",lola.db=NULL,assembly="hg19",pvalCut=0.01){
  res <- qtl.lola.enrichment(meth.qtl.res,type=type,assembly=assembly,lola.db=lola.db)
  plot <- lolaBarPlot(lolaDb=res$lola.db,lolaRes=res$lola.res,pvalCut=pvalCut)+my_theme
  return(plot)
}
