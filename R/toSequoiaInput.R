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
		indivNames <- inds(get(paste0(pops)))
		parent_IDs <- c(parent_IDs, indivNames)
		scores_all <- rbind(scores_all, scores(get(paste0(pops)))[indivNames,markerList])
	}
	#check parent names for NAs
	if(sum(is.na(parent_IDs), parent_IDs == "") > 0){
		stop("Error: one or more individuals in your baseline does not have a name. This may be caused by an error in IDFGEN when a Population only contains one individual.")
	}

	#get offspring
	offspring_IDs <- c()
	for (pops in mixturePops){
		indivNames <- inds(get(paste0(pops)))
		offspring_IDs <- c(offspring_IDs, indivNames)
		scores_all <- rbind(scores_all, scores(get(paste0(pops)))[indivNames,markerList])
	}
	#check offspring names for NAs
	if(sum(is.na(offspring_IDs), offspring_IDs == "") > 0){
		stop("Error: one or more individuals in your mixture does not have a name. This may be caused by an error in IDFGEN when a Population only contains one individual.")
	}

	#assign lifehistory fields
	life_history <- data.frame(ID = c(parent_IDs, offspring_IDs),
						  Sex = rep(1, length(c(parent_IDs, offspring_IDs))),
						  BY = c(rep(1, length(parent_IDs)), rep(2, length(offspring_IDs))),
						  	  stringsAsFactors = FALSE)

	#check names for uniqueness
	tempTable <- table(life_history$ID)
	if(sum(tempTable > 1) > 0){ # if not unique
		warning("The names of individuals are not unique across Populations. Concatenating Population and individual names to form names for Sequoia.")
		# remake names
		parent_IDs <- c()
		for (pops in baselinePops){
			parent_IDs <- c(parent_IDs, paste0(pops, "_", inds(get(paste0(pops)))))
		}
		offspring_IDs <- c()
		for (pops in mixturePops){
			offspring_IDs <- c(offspring_IDs, paste0(pops, "_", inds(get(paste0(pops)))))
		}
		life_history$ID <- c(parent_IDs, offspring_IDs)
		#check that all are now unique
		tempTable <- table(life_history$ID)
		if(sum(tempTable > 1) > 0){
			stop("Error: After concatenation of Population name, names are still not unique.")
		}
		rownames(scores_all) <- life_history$ID #so that rownames are output correectly in genotypes.txt
	}

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
			genos <- scores_all[,marker]
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
	write.table(genotypes, file=paste0(prefix, "genotypes.txt"), sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)
}

#' Build Sequioa Input files for Tom's template script
#'
#' This interfaces with EFGLmh objects to build files expected by Tom's template Sequoia script
#'   for single parent assignment (actual assignments, not simulations)
#'
#' @param x EFGLdata object with potential parents and offspring
#' @param baselinePops names of pops with potential parents
#' @param mixturePops names of pops with offspring
#' @param markerList list of markers to use, if not provided the default is to use all
#' @param prefix This is a string that will be used as the prefix for the two output files.
#'   Default is to have no prefix.
#' @return Writes files to the working directory that can be uploaded to the server and
#'   used to run single parentage with Sequoia using Tom's template script
#' @import EFGLmh
#' @export

toSequoiaInput_EFGLmh <- function(x, baselinePops, mixturePops, markerList = NULL, prefix = ""){
	if(is.null(markerList)) markerList <- getLoci(x)
	x <- x %>% removeLoci(lociRemove = getLoci(x)[!getLoci(x) %in% markerList])
	#get parents
	parent_IDs <- getInds(x, pops = baselinePops)

	#get offspring
	offspring_IDs <- getInds(x, pops = mixturePops)

	scores_all <- x$genotypes[match(c(parent_IDs, offspring_IDs), x$genotypes$Ind),3:ncol(x$genotypes)]

	#assign lifehistory fields
	life_history <- data.frame(ID = c(parent_IDs, offspring_IDs),
						  Sex = rep(1, length(c(parent_IDs, offspring_IDs))),
						  BY = c(rep(1, length(parent_IDs)), rep(2, length(offspring_IDs))),
						  stringsAsFactors = FALSE)

	#check names for uniqueness
	if(sum(table(life_history$ID) > 1) > 0) stop("The names of individuals are not unique.")

	#output metadata file
	write.table(life_history, file=paste0(prefix, "life_history.txt"), sep= "\t", row.names=FALSE, quote=FALSE)

	#build genotype file
	###define empty matrix with row names and column names
	genotypes <- matrix(nrow=nrow(life_history), ncol=length(markerList),
					dimnames=list(rows = c(parent_IDs, offspring_IDs), cols = markerList))
	to_remove <- c()
	###build genotype table
	for (i in 1:ncol(genotypes)){
		marker <- colnames(genotypes)[i]
		a1 <- scores_all[[paste0(marker, ".A1")]]
		a2 <- scores_all[[paste0(marker, ".A2")]]
		Alleles <- unique(c(a1, a2))
		Alleles <- Alleles[!is.na(Alleles)]
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
			genotypes[,i] <- (a1 == Alleles[1]) + (a2 == Alleles[1]) # count of first allele
			genotypes[is.na(genotypes[,i]),i] <- -9 # missing
		}
	}
	genotypes <- genotypes[,!(colnames(genotypes) %in% to_remove)]
	#output genotype file
	write.table(genotypes, file=paste0(prefix, "genotypes.txt"), sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)
}

