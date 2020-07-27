##########################################################################################
# compute_methQTL.R
# created: 2020-03-26
# creator: Michael Scherer
# ---------------------------------------------------------------------------------------
# Methods for calling correlation blocks from DNA methylation data.
##########################################################################################

#' computeCorrelationBlocks
#'
#' This function computes CpG correlation blocks from correlations of CpGs across samples by Louvian
#' clustering.
#'
#' @param meth.data A \code{data.frame} containing the methylation data with CpGs in the rows and samples in the columns.
#' @param annotation The genomic annotations of the CpG positions.
#' @param cor.threshold The correlation threshold used to discard edges from the correlation-based network.
#' @param sd.gauss Standard deviation of the Gauss distribution used to weight the distance
#' @param absolute.cutoff Absolute distance cutoff after which no methQTL interaction is to be considered.
#' @param max.cpgs Maximum number of CpGs used in the computation (used to save memory). 40,000 is a reasonable
#'             default for machines with ~128GB of main memory. Should be smaller for smaller machines and larger
#'             for larger ones.
#' @param assembly The assembly used
#' @param chromosome The chromosome for which correlation block calling is to be performed
#' @param segmentation If performed, DNA methylation based segmentation into PMDs/nonPMDs as a \code{GRanges}
#'        object using the 'epicPMDdetect' package
#' @return A list representing the clustering of CpGs into correlation blocks. Each element is a cluster, which contains
#'      row indices of the DNA methylation matrix that correspond to this cluster.
#' @details This method performs clustering of the correlation matrix obtaind from the DNA methylation matrix. Correlations
#'      are computed for each pair of CpGs across all the samples. We then compute a similarity matrix from this correlation
#'      matrix and set correlations lower than the given threshold to 0. In the next step, we weight the correlations
#'      by the distance between the CpGs: smaller distances get higher weights according to Gaussian distribution with
#'      mean 0 and standard deviation as specified above. Furthermore, similarities of CpGs that are further than
#'      \code{absolute.distance.cutoff} away from one another are discarded.
#'
#'      We then compute the associated weighted, undirected graph from the similarity matrix and execute Louvain clustering
#'      on the graph. The resulting clusters of CpGs are returned.
#'
#' @author Michael Scherer
#' @export
#' @import igraph
#' @import bigstatsr
computeCorrelationBlocks <- function(meth.data,
                                       annotation,
                                       cor.threshold=qtlGetOption("cluster.cor.threshold"),
                                       sd.gauss=qtlGetOption("standard.deviation.gauss"),
                                       absolute.cutoff=qtlGetOption("absolute.distance.cutoff"),
                                       max.cpgs=qtlGetOption("max.cpgs"),
                                       assembly="hg19",
                                       chromosome="chr1",
				       segmentation=NULL){
  logger.start("Compute correlation blocks")
  if(nrow(annotation)>max.cpgs){
    logger.info(paste("Split workload, since facing",nrow(annotation),"CpGs (Maximum is",max.cpgs,")"))
    bin.split <- round(nrow(annotation)/2)
    return(c(computeCorrelationBlocks(meth.data=meth.data[1:bin.split,],
                                        annotation=annotation[1:bin.split,],
                                        cor.threshold = cor.threshold,
                                        sd.gauss = sd.gauss,
                                        absolute.cutoff = absolute.cutoff,
                                        max.cpgs = max.cpgs,
					assembly=assembly,
                                        chromosome=chromosome,
					segmentation=segmentation),
             lapply(computeCorrelationBlocks(meth.data=meth.data[(bin.split+1):nrow(annotation),],
                                        annotation=annotation[(bin.split+1):nrow(annotation),],
                                        cor.threshold = cor.threshold,
                                        sd.gauss = sd.gauss,
                                        absolute.cutoff = absolute.cutoff,
                                        max.cpgs = max.cpgs,
					assembly=assembly,
                                        chromosome=chromosome,
					segmentation=segmentation),
                    function(x) x+bin.split
                    )
             ))
  }
  logger.start("Compute correlation matrix")
  if(qtlGetOption("correlation.type")=="pearson"){
    cor.all <- big_cor(as_FBM(t(as.matrix(meth.data)),type="double"))
  }else{
    cor.all <- cor(t(as.matrix(meth.data)),qtlGetOption("correlation.type"))
  }
  rm(meth.data)
  logger.completed()
  cor.all <- cor.all[,,drop=F]
  if(qtlGetOption("hdf5dump")){
    cor.all <- writeHDF5Array(cor.all)
  }
  rep.vals <- cor.all<cor.threshold
  if(qtlGetOption("hdf5dump")){
    rep.vals <- writeHDF5Array(rep.vals)
  }
  cor.all[rep.vals] <- 0
  genomic.positions <- annotation$Start
  logger.start("Compute pairwise distances")
  gc()
  pairwise.distance <- abs(as.data.frame(lapply(genomic.positions,function(x)x-genomic.positions)))
  logger.completed()
  rep.vals <- pairwise.distance>absolute.cutoff
  if(qtlGetOption("hdf5dump")){
    rep.vals <- writeHDF5Array(rep.vals)
  }
  cor.all[rep.vals] <- 0
  gc()
  logger.start("Weight distances")
  if(qtlGetOption("hdf5dump")){
    weighted.distances <- matrix(nrow=nrow(cor.all),ncol=ncol(cor.all))
    weighted.distances <- writeHDF5Array(weighted.distances)
    chunk.size <- 10000
    i <- 1
    while(i < nrow(cor.all)){
      if((i + chunk.size)>nrow(cor.all)){
        do.work <- i:nrow(cor.all)
        weighted.distances[do.work,] <- as.matrix(cor.all[do.work,])*dnorm(as.matrix(pairwise.distance[do.work,]),0,sd.gauss)
        break
      }
      do.work <- i:(i+chunk.size)
      weighted.distances[do.work,] <- as.matrix(cor.all[do.work,])*dnorm(as.matrix(pairwise.distance[do.work,]),0,sd.gauss)
      i <- i+chunk.size+1
    }
  }else{
    weighted.distances <- cor.all*dnorm(as.matrix(pairwise.distance),0,sd.gauss)
  }
  logger.completed()
  weighted.distances <- weightSegmentation(weighted.distances,annotation,segmentation)
  weighted.distances <- weightFunctionalAnnotation(weighted.distances,annotation,chromosome=chromosome)
  colnames(weighted.distances) <- as.character(1:ncol(weighted.distances))
  rownames(weighted.distances) <- as.character(1:nrow(weighted.distances))
  rm(rep.vals)
  rm(cor.all)
  gc()
  logger.start("Compute graph")
  graph.ad <- graph.adjacency(as.matrix(weighted.distances),weighted = T,mode = "undirected",diag=F)
  logger.completed()
  logger.start("Compute clustering")
  clust <- cluster_louvain(graph.ad)
  rm(weighted.distances)
  gc()
  logger.completed()
  logger.completed()
  return(lapply(groups(clust),function(x)as.numeric(x)))
}

#' weightFunctionalAnnotation
#'
#' This functions weights the distance matrix according to predefined functional annotations accoding to the ENSEMBL
#' regulatory build regions.
#'
#' @param input.matrix The input distance/correlation matrix on which the clustering will be performed.
#' @param genomic.annotation The genomic annotation of the CpGs used as a \code{data.frame}.
#' @param assembly The assembly used for the annotation
#' @param chromosome The chromosome used for the annotation and the data
#' @return The modified distance matrix with CpGs in the same functional annotation category prioritized
#' @author Michael Scherer
#' @noRd
weightFunctionalAnnotation <- function(input.matrix,
                                         genomic.annotation,
                                         assembly="hg19",
                                         chromosome=chromosome){
  if(qtlGetOption("use.functional.annotation")){
    logger.start("Weighting functional annotation")
    genomic.annotation <- makeGRangesFromDataFrame(genomic.annotation)
    ensembl.types <- c("ctcf","distal","dnase","proximal","tfbs","tss")
    for(ensembl.type in ensembl.types){
      ensembl.type <- paste0("ensembleRegBuildBP",ensembl.type)
      if(!(ensembl.type%in%rnb.region.types())){
        rnb.load.annotation.from.db(types=ensembl.type,assembly=assembly)
      }
      anno.ensembl <- rnb.get.annotation(ensembl.type,assembly = assembly)
      op <- findOverlaps(genomic.annotation,anno.ensembl)
      input.matrix[queryHits(op),queryHits(op)] <- input.matrix[queryHits(op),queryHits(op)]*qtlGetOption("functional.annotation.weight")
    }
    logger.completed()
  }
  return(input.matrix)
}

#' weightSegmentation
#'
#' This functions weights the distance matrix according to DNA methylation based segmentation into
#' PMDs/nonPMDs using the 'epicPMDdetect' package
#'
#' @param input.matrix The input distance/correlation matrix on which the clustering will be performed.
#' @param genomic.annotation The genomic annotation of the CpGs used as a \code{data.frame}.
#' @param segmentation The segmentation into PMDs/nonPMDs as a \code{GRanges object}.
#' @return The modified distance matrix with CpGs in the same category PMD/nonPMD priortized
#' @author Michael Scherer
#' @noRd
weightSegmentation <- function(input.matrix,
				genomic.annotation,
                                segmentation){
  if(qtlGetOption("use.segmentation")){
      logger.start("Weighting segmentation")
      if(is.null(segmentation)){
        return(input.matrix)
      }
      genomic.annotation <- makeGRangesFromDataFrame(genomic.annotation)
      op <- findOverlaps(genomic.annotation,segmentation)
      input.matrix[queryHits(op),queryHits(op)] <- input.matrix[queryHits(op),queryHits(op)]*qtlGetOption("functional.annotation.weight")
      logger.completed()
  }
  return(input.matrix)
}
