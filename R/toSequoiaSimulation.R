#' Build Sequioa Input files for \code{sequoiaSim}
#'
#' This interfaces with IDFGEN objects to build the parent_data_file expected by \code{sequoiaSim}
#'
#' @param baselinePops Poplist of potential parents
#' @param markerList list of markers to use
#' @return Writes a file to the working directory that can be uploaded to the server and
#'   used to run \code{sequoiaSim}
#' @export

toSequoiaSimulation <- function(baselinePops, markerList){
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
	#output file to upload to server
	write.table(parents, "parent_data.txt", sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE)
	print("A file with parent data and genotypes for use with the simulation script has been written to parent_data.txt")
}
