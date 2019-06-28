#' Estimate error expected with your baseline
#'
#' \code{sequoiaSim} runs simulations to estimate the error
#' expected when using sequoia for single parent assignments with your baseline
#'
#' This runs simulations given the parameters you input, see the separate writeup
#' for a full description of how the simulations are run.
#'
#' @param input_parameters This is a dataframe of input parameters to control the simulations. For example:
#'     parameters <- data.frame(Proportion_baseline_sampled = c(1, .9, .85, .75),
#'					Number_of_simulations_to_run =    c(100, 200, 200, 200),
#'					Number_of_offspring =             c(1000, 3000, 3000, 3000),
#'					Parent_expansion = c(TRUE, TRUE, TRUE, TRUE)
#'					)
#' @param LLR_min This is the minimum LLR you want to use to accept a single parent assignment
#' @param parent_data_file This is either a dataframe or a file path to a csv with a header and your baseline data
#'     in it. Each row is an individual, and columns are: PopulationName,IndividualName,Sex,SNP1,SNP2,SNP3,...
#'     Sex is either M or F, SNPs are one call per column and one character per allele (ie, AT, TT, AA)
#'     Individual names must be unique, even across populations. This file can be generated using
#'     the \code{toSequoiaSimulation} function.
#' @param min_genotyped This is the proportion of genotypes that must be non-missing (ie 0.9 for 90\% of genotypes) to
#'   include a simulated parent or offspring in the analysis
#' @param allele_error_rate This is the per allele genotyping error rate
#' @param Seq_MaxMismatch This is the MaxMismatch parameter to use for Sequoia. If not specified, the default is
#'   to use 5\% of the number of markers in your panel (rounded up to the nearest integer).
#' @param prefix This is the prefix to add to the output file name. Default is no prefix.
#' @return This function writes its output as a csv to the working directory afeter all simulations have finished.
#'   This function will print warnings() after each simulation, so if you have warnings in your R session from
#'   previous commands, these will be printed.
#' @importFrom stats rbinom rnbinom
#' @importFrom utils read.csv read.table write.table
#' @export

sequoiaSim <- function(input_parameters, LLR_min = 0.5, parent_data_file, min_genotyped = .9, allele_error_rate = .001,
				   Seq_MaxMismatch = NA, prefix = ""){

	#load parent data file
	if (is.character(parent_data_file)) {
		parent_data <- read.csv(parent_data_file, header = TRUE, colClasses = "character")
	} else {
		parent_data <- parent_data_file
	}

	###input error checking
	if(nrow(parent_data) < 1){
		stop("Error: empty parent data file")
	}
	for (i in 1:nrow(input_parameters)){
		if(input_parameters[i,1] < 0 || input_parameters[i,1] > 1){
			stop("Error: invalid Proportion_baseline_sampled value on row", i)
		}
		if(input_parameters[i,2] < 1){
			stop("Error: invalid Number_of_simulations_to_run value on row", i)
				}
		if(input_parameters[i,3] < 1){
			stop("Error: invalid Number_of_offspring value on row", i)
				}
		if(!(is.logical(input_parameters[i,4]))){
			stop("Error: invalid Parent_expansion value on row", i)
		}
	}
	#build list of simulations to run and proportion of parents to include
	sim_list <- c() #this is the proportion of baseline sampled in each simulation
	num_offspring <- c()
	parent_expansion <- c()
	for (i in 1:nrow(input_parameters)){
		sim_list <- c(sim_list, rep(input_parameters[i,1], input_parameters[i,2]))
		num_offspring <- c(num_offspring, rep(input_parameters[i,3], input_parameters[i,2]))
		parent_expansion <- c(parent_expansion, rep(input_parameters[i,4], input_parameters[i,2]))
	}
	#build blank output list
	output <- data.frame(Sim_number = 1:length(sim_list),
					 Proportion_baseline_sampled = sim_list,
					 Size_baseline = rep(0, length(sim_list)),
					 Number_true_parents_in_baseline = rep(0, length(sim_list)),
					 Number_offspring = rep(0, length(sim_list)),
					 Number_possible_assignments = rep(0, length(sim_list)),
					#####five possible groups to assess accuracy:
					###### Fish that could be assigned correctly (parent in baseline), and were assigned correctly
					###### Fish that could be assigned correctly (parent in baseline), and were assigned incorrectly
					###### Fish that could be assigned correctly (parent in baseline), and were NOT assigned
					###### Fish that could NOT be assigned correctly (parent NOT in baseline), and were assigned incorrectly
					###### Fish that could NOT be assigned correctly (parent NOT in baseline), and were NOT assigned correctly
					 Number_assignable_correct_assignments = rep(0, length(sim_list)),
					 Number_assignable_incorrect_assignments = rep(0, length(sim_list)),
					 Number_assignable_incorrect_UNassignments = rep(0, length(sim_list)),
					 Number_not_assignable_incorrect_assignments = rep(0, length(sim_list)),
					 Number_not_assignable_correct_UNassignments = rep(0, length(sim_list))
	)

	#change parent genotype data to A/B format with 0 still as missing genotype
	to_remove <- c()
	for (i in 4:ncol(parent_data)){
		alleles <- unique(c(substr(parent_data[,i],1,1), substr(parent_data[,i],2,2)))
		alleles <- subset(alleles, alleles != "0")
		if (length(alleles) == 0){
			print(paste0("Warning: all genotypes failed for locus ", colnames(parent_data)[i], ". Removing this locus from the dataset."))
			to_remove <- c(to_remove, colnames(parent_data)[i])
		}
		else if (length(alleles) == 1){
			print(paste("Warning: no variation in the parents for locus", colnames(parent_data)[i]))
			#all samples either has same geno or failed
			parent_data[parent_data[,i] != "00",i] <- "AA"
		}
		else if (length(alleles) > 2){
			print(paste0("Warning: more than two alleles found for locus ", colnames(parent_data)[i], ". Removing this locus from the dataset."))
			to_remove <- c(to_remove, colnames(parent_data)[i])
		}
		else{
			#change to A/B
			parent_data[parent_data[,i] == paste0(alleles[1], alleles[1]),i] <- "BB"
			parent_data[parent_data[,i] == paste0(alleles[1], alleles[2]),i] <- "AB"
			parent_data[parent_data[,i] == paste0(alleles[2], alleles[1]),i] <- "AB"
			parent_data[parent_data[,i] == paste0(alleles[2], alleles[2]),i] <- "AA"	#change to AA last b/c AA does occur in unchanged data
		}
	}
	parent_data <- parent_data[,!(colnames(parent_data) %in% to_remove)]	#remove loci with more than two alleles, or all fails
	parent_col_num <- ncol(parent_data)	#use this number repeatedly, so just calculate once

	#get failure to genotype rates, assume due to genotyping process and so is the same for all populations
	#will use this for simulating parents and offspring
	num_parent <- nrow(parent_data)	#save number of parents b/c use frequently
	# failure_rate_loci is a vector of genotyping failure rate (per individual) for each locus
	failure_rate_loci <- apply(parent_data[,4:parent_col_num], 2, function(x){sum(x == "00")/num_parent})

	#save parent data to reset after baseline expansion
	parent_data_orig <- parent_data

	#### beginning of for loop for multiple simulations
	for (sim_num in 1:length(sim_list)){
		print(paste("Starting simulation", sim_num, "of", length(sim_list)))
		#define parameters for this simulation
		mixture_size <- num_offspring[sim_num]
		baseline_expand <- parent_expansion[sim_num]
		prop_parent_sampled <- sim_list[sim_num]

		parent_data <- parent_data_orig
		pops <- unique(parent_data[,1])
		# simulate unsampled parents
		if(baseline_expand && prop_parent_sampled < 1){
			if(prop_parent_sampled < 0.01){	#avoid dividing by zero or overly large parent expansions
				prop_change <- 0.01
				print("Proportion_baseline_sampled was less than 0.01 and Parent_expansion was TRUE. Limiting parent population increase to 100x.")
			} else {
				prop_change <- prop_parent_sampled
			}
			add <- round(((nrow(parent_data)/prop_change) - nrow(parent_data)),0) #calculate number to add
			add_by_pop <- table(sample(pops, add, replace = TRUE, prob = table(parent_data[,1])[pops]))
			# call function to simulate new parents for each population
			#add in for each population, make prefixes different for each, select p_data to be only that population
			for (pop in pops){
				if(!(pop %in% names(add_by_pop)) || add_by_pop[pop] == 0){next}
				parent_data <- rbind(parent_data, sim_parents(p_data = parent_data[parent_data[,1] == pop,], num = add_by_pop[pop], prefix = paste0(pop, "_1"), min_geno = min_genotyped, fail_rates = failure_rate_loci))
			}
		}

		#fix genotyping errors so they aren't "inherited", which also increases genotyping error rate in offspring
		#within each population, keeps allele frequencies the same and assumes HWE when assigning genotypes to overwright failed loci
		new_parent_data <- parent_data # save as separate variable, because will use only for creating offspring genotypes
		for(pop in pops){	#loop through populations
			working_data <- new_parent_data[new_parent_data[,1] == pop,] # save to prevent mutiple calcs of accessing
			for (i in 4:parent_col_num){
				# determine if missing genotypes, and how many
				missing_num <- sum(working_data[,i] == "00")
				if(missing_num == 0) {next}
				#get parent allele frequencies
				parent_alleles <- c(substr(working_data[,i],1,1), substr(working_data[,i],2,2))
				A_freq <- sum(parent_alleles == "A")
				B_freq <- sum(parent_alleles == "B")
				both <- A_freq + B_freq
				if(both == 0){	#if all individuals in a population failed, use allele frequencies for entire parent data set
					parent_alleles <- c(substr(parent_data[,i],1,1), substr(parent_data[,i],2,2))
					A_freq <- sum(parent_alleles == "A")
					B_freq <- sum(parent_alleles == "B")
					both <- A_freq + B_freq
				}
				A_freq <- A_freq/both
				B_freq <- B_freq/both
				A_hom <- A_freq^2
				B_hom <- B_freq^2
				het <- 2*A_freq*B_freq
				working_data[working_data[,i] == "00",i] <- sample(c("AA", "AB", "BB"), missing_num, replace = TRUE, prob = c(A_hom, het, B_hom))
			}
			#will use new_parent_data for simulating offspring, but not for anything else
			new_parent_data[new_parent_data[,1] == pop,] <- working_data
		}
		rm(working_data)

		# simulate offspring

		#separate males and females for sampling
		dams <- parent_data[parent_data[,3] == "F",1:2]
		sires <- parent_data[parent_data[,3] == "M",1:2]

		#create offspring

		#choose females and number of offspring for each
		#parameters for negative binomial distribution (mu and dispersion parameter) of reproductive success
		#correspond to variance = 2.3 * mean, which was observed from SY2012 Snake River hatchery CHNK broodstock Females (mu=1.86, var=4.28)
		mu <- mixture_size/sum(parent_data[,3] == "F")	#mean is number of offspring chosen divided by number of females
		disp <- mu/1.3	#dispersion parameter to give Var = 2.3 * mean

		pairs <- data.frame(
			dam = dams[,2],
			sire = rep("empty", nrow(dams)),
			repro_succ = rnbinom(nrow(dams), size = disp, mu = mu), stringsAsFactors = FALSE
		)
		### remove zeroes
		pairs <- pairs[pairs[,3] > 0,]
		## choose sires from same populations
		for(pop in pops){
			bool_dams <- pairs[,1] %in% dams[dams[,1] == pop,2]
			pairs[bool_dams,2] <- sample(sires[sires[,1] == pop, 2], sum(bool_dams), replace = TRUE)
		}

		#sample offspring to get exactly the correct number
		mixture <-matrix(nrow = mixture_size, ncol = parent_col_num)	#matrix with name, true dam, true sire, genotypes
		colnames(mixture) <- c("Name", "true_dam", "true_sire", colnames(parent_data)[4:parent_col_num])
		mixture[,1] <- paste0("Offspring_", 1:mixture_size)
		mixture[,2:3] <- as.matrix(pairs[sample(1:nrow(pairs), mixture_size, replace = TRUE, prob = pairs[,3]),1:2])

		#save memory
		rm(dams)
		rm(sires)

		#create genotypes for offspring
		rownames(new_parent_data) <- new_parent_data[,2] # for fast accessing
		for (i in 1:mixture_size){	#iterate through offspring
			dam_genos <- new_parent_data[mixture[i,2],4:parent_col_num]
			sire_genos <- new_parent_data[mixture[i,3],4:parent_col_num]
			dam_pos <- sample(1:2, (parent_col_num - 3), replace = TRUE)	#choose allele one or two randomly for all loci
			sire_pos <- sample(1:2, (parent_col_num - 3), replace = TRUE) #choose allele one or two randomly for all loci
			mixture[i,4:parent_col_num] <- paste0(substr(dam_genos, dam_pos, dam_pos), substr(sire_genos, sire_pos, sire_pos))
		}

		#add in genotyping errors for offspring
		for(i in 4:parent_col_num){ # for each marker
			error <- rbinom(mixture_size, 2, allele_error_rate) #determine number of errors made by sampling binomial
			#make changes to observed genotypes
			for(j in which(error == 1)){
				if(mixture[j,i] == "AB"){
					mixture[j,i] <- sample(c("AA", "BB"), 1)
				} else {
					mixture[j,i] <- "AB"
				}
			}
			for(j in which(error == 2)){ # note that AB with two errors is AB
				if(mixture[j,i] == "AA"){
					mixture[j,i] <- "BB"
				} else {
					mixture[j,i] <- "AA"
				}
			}
		}

		## add in fails, preventing any below min_genotyped
		m_num <- parent_col_num-3
		fails <- sapply(failure_rate_loci, function(x){ # this produces matrix with 1 row per indiv, 1 col per marker with T if failed and F if not
			sample(c(TRUE,FALSE), mixture_size, replace = TRUE, prob = c(x, 1-x))
			})
		if(!is.matrix(fails)){fails <- t(as.matrix(fails))} # if only one sample, make a matrix
		fail_to_geno <- apply(fails, 1, sum)/m_num #determine genotyping success
		fail_to_geno <- fail_to_geno > (1-min_genotyped)
		# while not all successfully genotyped, re-sample fails, with a "ratcheting" mechanism
		while(sum(fail_to_geno) > 0){
			fails <- fails[!fail_to_geno,]
			fails2 <- sapply(failure_rate_loci, function(x){
				sample(c(TRUE,FALSE), sum(fail_to_geno), replace = TRUE, prob = c(x, 1-x))
			})
			if(!is.matrix(fails2)){fails2 <- t(as.matrix(fails2))}
			fails <- rbind(fails,fails2)
			fail_to_geno <- apply(fails, 1, sum)/m_num
			fail_to_geno <- fail_to_geno > (1-min_genotyped)
		}
		for(i in 1:m_num){
			mixture[fails[,i],(i+3)] <- "00"
		}

		#add in genotyping errors for parents
		num_parent <- nrow(parent_data)
		for(i in 4:parent_col_num){
			error <- rbinom(num_parent, 2, allele_error_rate)
			for(j in which(error == 1)){
				if(parent_data[j,i] == "AB"){
					parent_data[j,i] <- sample(c("AA", "BB"), 1)
				} else {
					parent_data[j,i] <- "AB"
				}
			}
			for(j in which(error == 2)){
				if(parent_data[j,i] == "AA"){
					parent_data[j,i] <- "BB"
				} else {
					parent_data[j,i] <- "AA"
				}
			}
		}

		#remove one parent from each cross from parent_data b/c only want to test single parent assignemnts
		true_parents <- unique(c(mixture[,2], mixture[,3]))	#list of true parents of offspring
		true_pairs <- pairs[pairs[,1] %in% true_parents,]	#only check dam b/c dams are not reused
		remove1 <- c()
		sires_to_save <- c()
		for (i in 1:nrow(true_pairs)){
			if (true_pairs[i,2] %in% remove1){	#only check sire b/c dams are not reused
				next
			}
			if (true_pairs[i,2] %in% sires_to_save) {
				remove1 <- c(remove1, true_pairs[i,1])
			}
			else{
				choice <- sample(1:2, 1)
				remove1 <- c(remove1, true_pairs[i,choice])
				if (choice == 1) {
					sires_to_save <- c(sires_to_save, true_pairs[i,2])
				}
			}
		}
		parent_data <- parent_data[!(parent_data[,2] %in% remove1),]

		#if parent expansion is TRUE
		#add in a simulated parent for each true parent that was removed; this keeps the baseline size the same as the input baseline
		if (baseline_expand){
			#calculate number of new parents needed
			add <- length(remove1)
			add_by_pop <- table(sample(pops, add, replace = TRUE, prob = table(parent_data[,1])[pops]))
			# call function to simulate new parents for each population
			#add in for each population, make prefixes different for each, select p_data to be only that population
			for (pop in pops){
				if(!(pop %in% names(add_by_pop)) || add_by_pop[pop] == 0){next}
				parent_data <- rbind(parent_data, sim_parents(p_data = parent_data[parent_data[,1] == pop,], num = add_by_pop[pop], prefix = paste0(pop, "_2"), min_geno = min_genotyped, fail_rates = failure_rate_loci))
			}
		}

		#remove parents from parent_data in accordance with Proportion_baseline_sampled
		if(prop_parent_sampled == 0){	#handle case of 0 differently
			print("Notice: Proportion_baseline_sampled is 0. Instead of randomly sampling individuals in the baseline to remove, all true parents will be removed from the baseline, and the rest of the individuals in the baseline will be included.")
			remove2 <- true_parents[!(true_parents %in% remove1)]
		} else {
			num_remove <- round((1 - prop_parent_sampled)*nrow(parent_data), 0)
			#removes parents at random, so actual proportion of true parents sampled will be variable, but mean of Proportion_baseline_sampled
			remove2 <- sample(parent_data[,2], num_remove, replace = FALSE)
		}
		parent_data <- parent_data[!(parent_data[,2] %in% remove2),]


		#make input and run sequoia
		seqInput <- makeSequioaSingleParent(parent_data[,c(2,4:ncol(parent_data))], mixture[,c(1,4:ncol(mixture))])
		if(is.na(Seq_MaxMismatch)){# allow user to specify, but this is default
			maxMis <- ceiling(0.05*ncol(seqInput[[2]]))
		} else {
			maxMis <- Seq_MaxMismatch
		}
		results <- sequoia::sequoia(GenoM=seqInput[[2]], LifeHistData=seqInput[[1]], MaxSibIter=0, MaxMismatch = maxMis, FindMaybeRel = FALSE)
		if (length(warnings()) > 0){
			print(warnings())
		}
		lh_data <- seqInput[[1]]
		rm(seqInput)
		#calculate summary statistics and save to output dataframe
		#calculation as variable and then assignment to matrix separately for error checking
		#calculate possible assignments, third term should always be zero, but left it in for error checking, or in case this script is repurposed
		pos_assign <- sum(mixture[,2] %in% lh_data[,1]) + sum(mixture[,3] %in% lh_data[,1]) - sum(mixture[,3] %in% lh_data[,1] & mixture[,2] %in% lh_data[,1])
		#calculate the number of correct and incorrect assignments made
		assigned <- results$PedigreePar 			#get assignments
		assigned <- assigned[!is.na(assigned$LLRdam),]	#filter unassigned
		assigned <- assigned[assigned$LLRdam >= LLR_min,]	#filter based on LLR
		####five possible groups to assess accuracy:
		###### Fish that could be assigned correctly (parent in baseline), and were assigned correctly
		###### Fish that could be assigned correctly (parent in baseline), and were assigned incorrectly
		###### Fish that could be assigned correctly (parent in baseline), and were NOT assigned
		###### Fish that could NOT be assigned correctly (parent NOT in baseline), and were assigned incorrectly
		###### Fish that could NOT be assigned correctly (parent NOT in baseline), and were NOT assigned correctly
		assignable_offspring <- mixture[mixture[,2] %in% lh_data[,1] | mixture[,3] %in% lh_data[,1], 1]
		unassignable_offspring <- mixture[!(mixture[,1] %in% assignable_offspring), 1]
		assignable_assigned_bool <- assigned[,1] %in% assignable_offspring	#used multiple times below

		assignable_correct <- sum(assigned[,2] == mixture[match(assigned[,1], mixture[,1]),2] | assigned[,2] == mixture[match(assigned[,1], mixture[,1]),3])
		assignable_incorrect <- sum(assigned[assignable_assigned_bool,2] != mixture[match(assigned[assignable_assigned_bool,1], mixture[,1]),2] & assigned[assignable_assigned_bool,2] != mixture[match(assigned[assignable_assigned_bool,1], mixture[,1]),3])
		assignable_unassigned <- sum(!(assignable_offspring %in% assigned[,1]))
		unassignable_assigned <- sum(unassignable_offspring %in% assigned[,1])
		unassignable_unassigned <- sum(!(unassignable_offspring %in% assigned[,1]))


		#add data to output data frame
		output[sim_num,3] <- sum(lh_data[,3] == 1)	#size of baseline fed to Sequoia
		output[sim_num,4] <- sum(true_parents %in% lh_data[,1])	#number of true parents in baseline
		output[sim_num,5] <- sum(lh_data[,3] == 2)	#number of offspring (mixture size)
		output[sim_num,6] <- pos_assign	#number of possible correct assignments
		output[sim_num,7] <- assignable_correct	# Fish that could be assigned correctly (parent in baseline), and were assigned correctly
		output[sim_num,8] <- assignable_incorrect # Fish that could be assigned correctly (parent in baseline), and were assigned incorrectly
		output[sim_num,9] <- assignable_unassigned	# Fish that could be assigned correctly (parent in baseline), and were NOT assigned
		output[sim_num,10] <- unassignable_assigned	# Fish that could NOT be assigned correctly (parent NOT in baseline), and were assigned incorrectly
		output[sim_num,11] <- unassignable_unassigned	# Fish that could NOT be assigned correctly (parent NOT in baseline), and were NOT assigned correctly
	}
	#output combined output file
	write.table(output, paste0(prefix, "Sequoia_simulation_summary.txt"), sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
	return("Simulations complete")
}
