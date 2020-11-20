##########################################################################################
# GWAPMap.R
# created: 2020-07-30
# creator: Michael Scherer
# ---------------------------------------------------------------------------------------
# methods to import GWAS map output and for working with the results
##########################################################################################

#' parseGWASMapOutput
#'
#' This method takes as input a path to a CSV file containing the GWASMap output
#'
#' @param file The path to the GWASMap CSV file
#' @param sep The table separator
#' @param dec The decimal point separator
#' @return A list object with each element being the first line of the GWASMap output file.
#'      Each of the distinct elements in the first column is used as a category and the associated
#'      table is parsed separately as a table.
#' @noRd
#' @author Michael Scherer
parseGWASMapOutput <- function(file,
				sep=",",
				dec="."){
	first.line <- unlist(read.table(file,sep=sep,nrows=1))
	first.line[is.na(first.line)] <- ""
	has.info <- sapply(first.line,nchar)
	parts <- sum(has.info>0)
	res <- list()
	i <- 1
	done <- rep(FALSE,length(first.line))
	while(!all(done)){
		done[1:i] <- TRUE
		n <- (min(c(which(has.info>0&!done),length(first.line)))-1)
		cols.to.read <- rep("NULL",length(first.line))
		cols.to.read[i:n] <- NA
		dat.table <- read.table(file,
				sep=sep,
				header=TRUE,
				skip=1,
				colClasses=cols.to.read,
				dec=dec)
		res[[first.line[i]]] <- dat.table
		i <- (n+1)
	}
	return(res)
}
