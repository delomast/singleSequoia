# function to make sequoia input for single parent assignment only
#takes baseline and mixture as dataframes with order of columns (column names ignored): IndividualName, marker1, marker2, ...
# markers are bialleleic SNPs with one call per column, AB format with "00" as unknown - potential genotypes of "AA", "AB", "BA", "BB", "00"
# assumes marker orders are the same
# returns list of life history data, genotype data

makeSequioaSingleParent <- function(baseline, mixture){
	# check for same column number
	if(ncol(baseline) != ncol(mixture)){
		cat("\nError: number of columns is not the same in the baseline and mixture inputs. Exiting.\n")
		return()
	}
	#check for same column names: throw a warning but attempt to continue
	colname_bool <- colnames(baseline) != colnames(mixture)
	if(sum(colname_bool) != 0){
		cat("\nWarning: column names of baseline and mixture are not identical. This function assumes")
		cat("that the first column of both is the sample name and the remaining columns are markers with the order of markers being identical")
		cat("in the two datasets.\n")
		cat("Columns with different names are", which(colname_bool), "\n")
		colnames(mixture) <- colnames(baseline)
	}

	#create sequoia input files
	#create life history dataframe
	lh_data <- data.frame(ID = c(baseline[,1], mixture[,1]),
					  Sex = rep(1, (nrow(baseline) + nrow(mixture))),
					  BY = c(rep(1, nrow(baseline)), rep(2, nrow(mixture))), stringsAsFactors = FALSE
				)
	#create genotype matrix
	geno_data <- rbind(baseline[,2:ncol(baseline)], mixture[,2:ncol(baseline)])
	for (i in 1:ncol(geno_data)){
		geno_data[geno_data[,i] == "00",i] <- -9
		geno_data[geno_data[,i] == "AA",i] <- 2
		geno_data[geno_data[,i] == "AB",i] <- 1
		geno_data[geno_data[,i] == "BA",i] <- 1
		geno_data[geno_data[,i] == "BB",i] <- 0
	}
	geno_data <- as.matrix(geno_data)
	# convert to numeric matrix
	geno_data <- apply(geno_data, 2, as.numeric)
	rownames(geno_data) <- lh_data[,1]
	return(list(lh_data, geno_data))
}
