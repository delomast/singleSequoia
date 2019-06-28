# simulation of parents function
# takes data for only one population
# assumes HWE

sim_parents <- function(p_data, num, prefix, min_geno, fail_rates){
	m_num <- ncol(p_data) - 3
	all_fail <- 0 #keep track of failed loci in pop for preventing min_geno violations
	## define blank output
	new_parents <- matrix(nrow = num, ncol=m_num+3)
	##metadata
	new_parents[,1] <- rep(p_data[1,1], num) #pop name
	new_parents[,2] <- paste0(prefix, "_NewParent_", 1:num) #individual name
	new_parents[,3] <- sample(p_data[,3], num, replace = TRUE) #sex
	##genotypes
	### calculate allele frequencies
	freqs <- apply(p_data[,4:(m_num+3)], 2, function(x){
		x <- c(substr(x,1,1), substr(x,2,2))
		both <- sum(x == "A" | x =="B")
		if(both == 0){return(c(0,0))}
		return(c(sum(x == "A")/both, sum(x == "B")/both))})
	for(i in 1:m_num){
		if(sum(freqs[,i]) == 0){
			## if all in pop failed, make all new ones fail
			new_parents[,(i+3)] <- "00"
			all_fail <- all_fail + 1
			next
		}
		#calculating HWE probs is much faster than sampling alleles individually
		Ahom <- freqs[1,i]^2
		Bhom <- freqs[2,i]^2
		het <- 2*freqs[1,i]*freqs[2,i]
		new_parents[,(i+3)] <- sample(c("AA", "AB", "BB"), num, replace = TRUE, prob = c(Ahom, het, Bhom))
	}
	## add in fails, preventing any below min_geno
	fails <- sapply(fail_rates, function(x){ # this produces matrix with 1 row per indiv, 1 col per marker with T if failed and F if not
		sample(c(TRUE,FALSE), num, replace = TRUE, prob = c(x, 1-x))
		})
	if(!is.matrix(fails)){fails <- t(as.matrix(fails))} # if only one sample, make a matrix
	fail_to_geno <- (apply(fails, 1, sum) + all_fail)/m_num #determine genotyping success
	fail_to_geno <- fail_to_geno > (1-min_geno)
	while(sum(fail_to_geno) > 0){ # while not all successfully genotyped, re-sample fails, with a "ratcheting" mechanism
		fails <- fails[!fail_to_geno,]
		fails2 <- sapply(fail_rates, function(x){
			sample(c(TRUE,FALSE), sum(fail_to_geno), replace = TRUE, prob = c(x, 1-x))
		})
		if(!is.matrix(fails2)){fails2 <- t(as.matrix(fails2))}
		fails <- rbind(fails,fails2)
		fail_to_geno <- (apply(fails, 1, sum) + all_fail)/m_num
		fail_to_geno <- fail_to_geno > (1-min_geno)
	}
	for(i in 1:m_num){
		new_parents[fails[,i],(3+i)] <- "00"
	}
	colnames(new_parents) <- colnames(p_data)
	return(new_parents)
}
