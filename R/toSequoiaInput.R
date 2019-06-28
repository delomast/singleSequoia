#' Build Sequioa Input files for Tom's template script
#'
#' This interfaces with IDFGEN objects to build files expected by Tom's template Sequoia script
#'   for single parent assignment (actual assignments, not simulations)
#'
#' @param baselinePops Poplist of potential parents
#' @param mixturePops Poplist of offspring
#' @param markerList list of markers to use
#' @param prefix This is a string that will be used as the prefix for the two output files.
#'   Default is to have no prefix.
#' @return Writes files to the working directory that can be uploaded to the server and
#'   used to run single parentage with Sequoia using Tom's template script
#' @export

toSequoiaInput <- function(baselinePops, mixturePops, markerList, prefix = ""){
	scores_all <- matrix(nrow=0, ncol=length(markerList))
	colnames(scores_all) <- markerList
	#get parents
	parent_IDs <- c()
	for (pops in baselinePops){
		parent_IDs <- c(parent_IDs, inds(get(paste0(pops))))
		scores_all <- rbind(scores_all, scores(get(paste0(pops)))[,markerList])
	}
	#get offspring
	offspring_IDs <- c()
	for (pops in mixturePops){
		offspring_IDs <- c(offspring_IDs, inds(get(paste0(pops))))
		scores_all <- rbind(scores_all, scores(get(paste0(pops)))[,markerList])
	}
	#assign lifehistory fields
	life_history <- data.frame(ID = c(parent_IDs, offspring_IDs),
						  Sex = rep(1, length(c(parent_IDs, offspring_IDs))),
						  BY = c(rep(1, length(parent_IDs)), rep(2, length(offspring_IDs))),
						  	  stringsAsFactors = FALSE)
	#output metadata file
	write.table(life_history, file=paste0(prefix, "life_history.txt"), sep= "\t", row.names=FALSE, quote=FALSE)
	#build genotype file
	###define empty matrix with row names and column names
	genotypes <- matrix(nrow=length(life_history$ID), ncol=length(markerList), dimnames=list(rows = unlist(life_history$ID), cols = markerList))
	to_remove <- c()
	###build genotype table
	for (i in 1:ncol(genotypes)){
		marker <- colnames(genotypes)[i]
		Alleles <- unique(c(substr(scores_all[,marker], 1, 1), substr(scores_all[,marker], 2, 2)))
		Alleles <- subset(Alleles, Alleles != "0")
		if (length(Alleles) == 0){
			print(paste("Warning: all genotypes failed for locus", marker))
			genotypes[,i] <- -9
		}
		if (length(Alleles) == 1){
			print(paste("Warning: no variation in the population for locus", marker))
			genotypes[,i] <- 2
		}
		if (length(Alleles) > 2){
			print(paste0("Warning: more than two alleles found for locus ", marker, ". Removing this locus from the dataset."))
			to_remove <- c(to_remove, marker)
		}
		else{		#if two alleles present
			genos <- scores_all[rownames(genotypes),marker]
			genotypes[genos == paste0(Alleles[1], Alleles[1]),i] <- 2
			genotypes[genos == paste0(Alleles[1], Alleles[2]),i] <- 1
			genotypes[genos == paste0(Alleles[2], Alleles[1]),i] <- 1
			genotypes[genos == paste0(Alleles[2], Alleles[2]),i] <- 0
			genotypes[genos == "00",i] <- -9
			rm(genos)
		}
	genotypes <- genotypes[,!(colnames(genotypes) %in% to_remove)]
	}
	#output genotype file
	write.table(genotypes, file=paste(prefix, "genotypes.txt"), sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)
}
