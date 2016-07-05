import os, sys, random, math, time, pickle

##########################################################################################
#	Name: findRhoMaximizingFitness
#
#	Assumptions: When calculating the fitness given a rho, the environmental optimum is 
#		     assumed to be zero. The reason for this is that this function is called
#		     at the population is being initialized and the environmental optimum 
#		     is 0 and the alpha and beta genes' sums are also 0. 
#
#	Purpose: This function exists so that, given an lDel gene, the user can find the rho
#	         value that maximizes the fitness of an individual. This is done by starting
#		 rho at .2. Rho is then multiplied by .99 until it is less than 10**-13, For 
#		 each one of these intermediate rho values, the fitness is calculated with the
#		 the given lDel gene. The highest fitness and the rho corresponding to it are
#		 stored.
#
#	Arguments: s - an integer corresponding to the fitness cost associated with 
#		       proofreading
#	           lDelGene - a list of 1's and 0's whose optimal rho value will be found
#		   pDel -  the probability of a loci in an lDel gene going from benign to 
#			   deleterious. This is only used to calculate the fitness.
#		   pMinusDel - the probability of a loci in lDel gene going from deleterious
#			       to benign. This is only used to calculate fitness.
#			       mulDel - the probability of a given loci in an lDel gene mutating. This is
#						passed to calculate fitness.
#		   alphaGene - the alpha gene that will be used to calculate the fitness
#	           betaGene - the beta gene used to calculate the fitness 
#
#	Returns: A float value in the interval (10**-13, 2] is returned. This float 
#		 corresponds to the rho value that would maximize fitness given the lDel gene
#		 passed by the user. 
#
##########################################################################################
def findRhoMaximizingFitness(s, lDelGene, pDel, pMinusDel, mulDel, alphaGene, betaGene):
	# Start r1ho at .2
	currentMaxFitness = 0
	rhoCorrespondingToMaxFit = 0
	lDelGeneLength = len(lDelGene)
	rho = .2
	
	# Multiply rho by .99  until rho < 10**-13
	while rho > 10**-13:
		# Find the fitness given the current rho and the assumption that the environmental
		# optimum is 0
		fitness = getFitness(s, lDelGene, pDel, pMinusDel, mulDel, rho, alphaGene, 
		                     betaGene, 0)
		                     
	 	# Keep track of the highest fitness and the rho corresponding to it
		if fitness > currentMaxFitness:
			currentMaxFitness = fitness
			rhoCorrespondingToMaxFit = rho
		
		# Increment rho
		rho *= .99
		
	# Return the rho corresponding to the maximum fitness found
	return(rhoCorrespondingToMaxFit)

##########################################################################################
#       Name: getEnvScore
#
#       Assumptions: It is assumed that the alpha sequences is at most as long as the lDel
#	             sequence and that the beta sequences is exactly as long as the alpha
#		     sequence. 
#
#       Purpose: The environmental score is calculated given alpha, lDel, and beta 
#		 sequences and a read through rate rho. Let's say an individual has k 
#		 cryptic sequences and l <= k alpha and l  beta loci. The environmental 
#		 readiness of an individual is then 
#		 sum(alpha[i] + (1 - lDel[i]) * tho * beta[i]) for i = 1, 2, ... , k. 
#
#       Arguments: Two lists, alphaGene and betaGene, of equal length containg floats. 
#		   Rho is a float between 0 and 1. lDelGene is a list of 1s and 0s.
#
#       Returns: A float corresponding in the formula in the purpose section of this
#		 block.
#
##########################################################################################
def getEnvScore(alphaGene, lDelGene, betaGene, rho):
	if len(alphaGene) != len(betaGene):
		print("Alpha and beta are different lengths.")
		print("Alpha: ", alphaGene)
		print("Beta: ", betaGene)
		exit(1)
	
	envScore = 0.0
	for i in range(len(alphaGene)):
		if lDelGene[i]:
			envScore += alphaGene[i]
		else:
			envScore += alphaGene[i] + rho * betaGene[i]
	return envScore

##########################################################################################
#	Name: getFitness
#
#	Assumptions: None
#
#	Purpose: getFitness exists so that one can find the fitness given a collection of 
#		 genetic information, the most important being an lDel gene, a rho value,
#		 an alpha gene and a beta gene.  There are four components of fitness. The
#		 first three components are calculated according the equations in the 2011
#		 Rajon-Masel paper. The last is calculated from the alpha and beta gene beta
#		 gene. A sum over the two lists will be calculated an used as an
#		 environmental readiness that will be compared to the optimum. The sum 
#		 includes alpha genes and possibly beta genes; whether or not these beta
#		 genes are included depends on a strip of lDels at the starrt  of the
#		 lDel gene.
#			
#	Arguments: s - an integer corresponding to the fitness cost associated with 
#		       proofreading
#		   lDelGene - a list of 1's and 0's whose optimal rho value will be found
#		   pDel -  the probability of a loci in an lDel gene going from benign to 
#			   deleterious.
#		   pMinusDel - the probability of a loci in lDel gene going from deleterious
#			       to benign.
#		   mulDel - the probability of a given loci in an lDel gene mutating.
#		   rho - the read through rate of a stop codon 
#		   alphaGene - a list containing ten floats representing 10 different 
#			       pre-stop codon trait values of an individuals
#		   betaGene - a list containing ten floats representing 10 different
#			      post-stop codon trait values
#		   envOptimum - the optimum environmental trait value for an individual
#
#	Returns: A float corresponding to the fitness of an individual is returned. It will
#		 will be in the interval [0, 1] and is the product of the four components of 
#		 fitness.
#
##########################################################################################
def getFitness(s, lDelGene, pDel, pMinusDel, mulDel, rho, alphaGene, betaGene, envOptimum):
	# Calculate the first three the fitness components
	lDelGeneLength = len(lDelGene)
	envScore = getEnvScore(alphaGene, lDelGene, betaGene, rho)

	tempDelFitness = max(0, 1 - s * ((rho *  lDelGene.count(1) / float(lDelGeneLength)) + \
	                     (1 - lDelGene.count(1) / float(lDelGeneLength)) * rho**2 * \
	                     pDel / (pDel + pMinusDel)))		
	permDelFitness = (1 - 23.0/9.0 * mulDel)**lDelGene.count(1)
	proofFitness =  1 / (1 - math.log(rho) * 10**-2.5)
	envFitness = math.exp(-((envScore - envOptimum)**2 / (2 * .5**2)))
	
	# Multiply all the components together and return the result 
	fitness = tempDelFitness * permDelFitness * proofFitness * envFitness
	return fitness

##########################################################################################
#	Name: initializePopulation
#
#	Assumptions: It is assumed that all of the individual in the population being 
#		     initialized should have the rho value that makes their fitness the 
#		     highest. It is also assumed the hat each individual has five pieces 
#		     of information the that define them: a rho gene, an lDel, a fitness, a 
#		     beta gene, and an alpha gene.
#
#	Purpose: This function exists to assign genes and a fitness value to each individual 
#		 in a population of a desired size. The genes of an individual include:
#                rho, lDel, alpha, beta. lDel is a list of 1's and 0's, rho is a float, and 
# 		 alpha and beta are lists of float values. An initial proportion for 1's
# 		 is given by the user, and each individual gets roughly that proportion of 
#		 1's in their lDel gene. Once their lDel gene is built, the rho value 
# 		 maximizing their fitness (see findRhoMaximizingFitness) is calculated. 
# 		 Then the individuals alpha and beta genes are initialized to all 0's. 
# 		 Finally, all of these are stored in a list and stored into the entry of 
#		 a larger list. Their index in the larger list is how they are referred to 
# 		 later on.
#
#	Arguments: population - the list where the genes and fitness of an individual are 
# 			        are stored
# 		   pOne - the proportion of 1's in the lDel gene for each individual in the 
#			  population
#		   s - an integer corresponding to the fitness cost associated with 
#		       proofreading
# 		   lDelGeneLength - the length of the lDel gene
#		   pNonDelToDel - the probability of a 1 changing to a 0 in lDel via mutation
#		   pDelToNonDel - the probability of a 0 changing to a 1 in lDel via mutation
#		   plDelLociMutation - the probability of a 1 or a 0 changing to a 0 or 1 
# 				       in lDel via mutation
# 		   alphaGeneLength - the number of floats in each alpha gene
#		   betaGeneLength - the number of floats in each beta gene
#
#	Returns: Nothing
#
##########################################################################################
def initializePopulation(populationSize, pOne, s, lDelGeneLength, 
                         pNonDelToDel, pDelToNonDel, plDelLociMutation, alphaGeneLength,
                         betaGeneLength):

	# Store the myPop as a list
	myPop = []                
	# Creates rho and lDel genes, in the form of two arrays, for each individual 
	# and calculates their fitness
	for i in range(populationSize):
		# Randomly lDel. Make approximately pOne of the lDels 1's and the rest 0's
		lDelGene = []
		for j in range(lDelGeneLength):
			if random.random() < pOne:
				lDelGene.append(1)
			else:
				lDelGene.append(0)
	
		# Initialize alpha and beta genes to zeroes
		alphaGene = [0.0 for i in range(alphaGeneLength)]
		betaGene = [0.0 for i in range(betaGeneLength)]
		
		# Find the rho value maximizing fitness given the lDel of the individual
		individualRho = findRhoMaximizingFitness(s, lDelGene, pNonDelToDel, pDelToNonDel, 
		                                         plDelLociMutation, alphaGene, betaGene)
		
		# Calculate considering the lDel gene and rho gene 
		fitness = getFitness(s, lDelGene, pNonDelToDel, pDelToNonDel, plDelLociMutation, 
							 individualRho, alphaGene, betaGene, 0)
	
		# An individual consists of a rho gene, an lDel gene, and a fitness
		# Each is an entry in an myPopSize size array
		myPop.append([individualRho, lDelGene, fitness, alphaGene, betaGene])
	
	# Return a reference to the newly made population	
	return myPop
