#' Build Sequioa Input files for \code{sequoiaSim}
#'
#' This interfaces with IDFGEN objects to build the parent_data_file expected by \code{sequoiaSim}
#'
#' @param baselinePops Poplist of potential parents
#' @param markerList list of markers to use
#' @param prefix This is a string that will be used as the prefix for the output files.
#'   Default is to have no prefix.
#' @return Writes a file to the working directory that can be uploaded to the server and
#'   used to run \code{sequoiaSim}
#' @export

toSequoiaSimulation <- function(baselinePops, markerList, prefix = ""){
	print("Note: this function assumes you have already removed all duplicate and failed individuals from your baseline populations.")
	parents <- matrix(nrow=0, ncol=(3 + length(markerList)))	#initiate emtpy matrix to fill with output
	colnames(parents) <- c("Population", "Name", "Sex", markerList)
	for (pop in baselinePops){	#iterate through populations
		all_genos <- scores(get(paste0(pop)))	#get genotypes so don't have to keep calling the function
		#build list of sex for each individual
		sex_marker <- 0
		for (i in 1:ncol(all_genos)){
			if (sum(all_genos[,i] == "XX") >= 1){
				sex_marker <- i
				break
			}
		}
		if ("Gender" %in% colnames(metaData(get(paste0(pop))))){
			pheno_sex <- metaData(get(paste0(pop)))[rownames(all_genos),"Gender"]
			pheno_found <- TRUE
		} else {
			pheno_found <- FALSE
		}
		sex <- rep("U", nrow(all_genos))
		if (sex_marker != 0){
			sex[all_genos[,sex_marker] == "XX"] <- "F"
			sex[all_genos[,sex_marker] == "XY" | all_genos[,sex_marker] == "YX" | all_genos[,sex_marker] == "YY"] <- "M"
		}
		if (pheno_found){
			sex[sex == "U" & pheno_sex == "F"] <- "F"
			sex[sex == "U" & pheno_sex == "M"] <- "M"
		}
		if(sum(sex == "U") > 0){
			sex[sex == "U"] <- sample(c("M", "F"), sum(sex == "U"), replace = TRUE)
		}

		parents_temp <- cbind(pop, rownames(all_genos), sex, all_genos[,markerList])
		colnames(parents_temp) <- c("Population", "Name", "Sex", markerList)
		parents <- rbind(parents, parents_temp)
	}
	#check names for NAs
	if(sum(is.na(parents[,2]), parents[,2] == "") > 0){
		stop("Error: one or more individuals in your baseline does not have a name. This may be caused by an error in IDFGEN when a Population only contains one individual.")
	}
	#check names for uniqueness
	tempTable <- table(parents[,2])
	if(sum(tempTable > 1) > 0){ # if not unique
		warning("The names of individuals are not unique across populations. Concatenating Population and individual names to form names for Sequoia.")
		#concatenate pop and individual names
		parents[,2] <- paste0(parents[,1], "_", parents[,2])
		#double check for uniqueness
		tempTable <- table(parents[,2])
		if(sum(tempTable > 1) > 0){
			stop("Error: After concatenation of Population name, names are still not unique.")
		}
	}
	#output file to upload to server
	write.table(parents, paste0(prefix, "parent_data.txt"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE)
	cat("\nA file with parent data and genotypes for use with the simulation script has been written to", paste0(prefix,"parent_data.txt"), "\n")
}

#' Build Sequioa Input files for \code{sequoiaSim}
#'
#' This interfaces with EFGLmh objects to build the parent_data_file expected by \code{sequoiaSim}
#'
#' @param x EFGLdata object with potential parents and offspring
#' @param baselinePops names of pops with potential parents
#' @param markerList list of markers to use
#' @param prefix This is a string that will be used as the prefix for the output files.
#'   Default is to have no prefix.
#' @return Writes a file to the working directory that can be uploaded to the server and
#'   used to run \code{sequoiaSim}
#' @export

toSequoiaSimulation_EFGLmh <- function(x, baselinePops, markerList = NULL, prefix = ""){
	print("Note: this function assumes you have already removed all duplicate and failed individuals from your baseline populations.")
	if(is.null(markerList)) markerList <- getLoci(x)
	x <- x %>% removeLoci(lociRemove = getLoci(x)[!getLoci(x) %in% markerList])
	parents <- matrix(nrow=0, ncol=(3 + length(markerList)))	#initiate emtpy matrix to fill with output
	colnames(parents) <- c("Population", "Name", "Sex", markerList)
	for (pop in baselinePops){	#iterate through populations
		all_genos <- x$genotypes %>% filter(Pop == pop)
		#build list of sex for each individual
		sex_marker <- 0
		for (i in 3:ncol(all_genos)){
			if (sum(all_genos[[i]] == "X", na.rm = TRUE) >= 1){
				sex_marker <- i
				break
			}
		}
		if ("Gender" %in% colnames(x$metadata)){
			pheno_sex <- x$metadata$Gender[match(all_genos$Ind, x$metadata$Ind)]
			pheno_found <- TRUE
		} else {
			pheno_found <- FALSE
		}
		sex <- rep("U", nrow(all_genos))
		if (sex_marker != 0){
			sex[!is.na(all_genos[[sex_marker]]) & all_genos[[sex_marker]] == "X" & all_genos[[sex_marker + 1]] == "X"] <- "F"
			sex[!is.na(all_genos[[sex_marker]]) & (all_genos[[sex_marker]] == "Y" | all_genos[[sex_marker + 1]] == "Y")] <- "M"
		}
		if (pheno_found){
			sex[sex == "U" & !is.na(pheno_sex) & pheno_sex == "F"] <- "F"
			sex[sex == "U" & !is.na(pheno_sex) & pheno_sex == "M"] <- "M"
		}
		if(sum(sex == "U") > 0){
			sex[sex == "U"] <- sample(c("M", "F"), sum(sex == "U"), replace = TRUE)
		}

		parents_temp <- data.frame(Population = pop, Name = all_genos$Ind, Sex = sex)
		for(m in markerList){
			a1 <- all_genos[[paste0(m, ".A1")]]
			a2 <- all_genos[[paste0(m, ".A2")]]
			a1[is.na(a1)] <- "0"
			a2[is.na(a2)] <- "0"
			parents_temp <- cbind(parents_temp, data.frame(m = paste0(a1, a2)))
		}
		parents <- rbind(parents, parents_temp)
	}
	#check names for NAs
	if(sum(is.na(parents[,2]), parents[,2] == "") > 0){
		stop("Error: one or more individuals in your baseline does not have a name.")
	}
	#check names for uniqueness
	tempTable <- table(parents[,2])
	if(sum(tempTable > 1) > 0){ # if not unique
		stop("The names of individuals are not unique across populations.")
	}
	#output file to upload to server
	write.table(parents, paste0(prefix, "parent_data.txt"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE)
	cat("\nA file with parent data and genotypes for use with the simulation script has been written to", paste0(prefix,"parent_data.txt"), "\n")
}

