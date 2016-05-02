import os, sys, random, math, time, pickle
				
##########################################################################################
#	Name: pickDeadIndiv
#
#	Assumptions: It is assumed that list holding individuals, the ith entry in the list
#		     corresponds to the ith individual in the population. Furthermore, it is 
# 		     assumed the fitnessIndex entry in the ith index of the population array
# 		     is the fitness of the individual. Finally, it is assumed that the 
#		     fitness of an individual is less than one; if all indviduals in the 
#		     population have a fitness greater or equal to one, then an infinite 
#		     loop will occur.
#
#	Purpose: The simulation this function is used for is one that keeps a constant 
#		 population size. Therefore, to produce an offspring, an individual must
#		 "die" or be replaced. This function picks the dead indvidual. If selection
#		 is enabled, then a random index in the population representing an individual
#		 will be chosen. If selection is not on, then the first individual to have 
#		 a fitness less than a random number in the range [0, 1) will be picked as 
#		 the dead indvidual.
#
#	Arguments: selection - a boolean corresponding to whether or not selection will 
#			       dictate if the person who dies is random or someone whose
#			       fitness was not high enough
#		   population - a list containing the individuals of the population and all of
#				their genetic information
#		   fitnessIndex - the index in the an in individuals entry in the population 
#				  list holding the fitness of that individual. 
#				  I.e. population[indvidual][fitnessIndex] gives individual's
#				  fitness
#
#	Returns: An integer is returned. This integer is greater than zero and less than the
#		 number of individuals in the population. This integer corresponds to which
#		 individual has been chosen to die.
#
##########################################################################################
def pickDeadIndiv(selection, population, fitnessIndex):
	# Calculate how big the population is
	populationSize = len(population)
	
	# If selection is on, consider an individual's fitness
	if selection:
		# Pick a random individual until their fitness is lower than a random float 
		# between zero and one.
		while True:
			potientialDeadIndiv = int(random.random() * populationSize) 
			
			# If someone's fitness less than a random float, return their index in 
			# the population list
			if population[potientialDeadIndiv][fitnessIndex] < random.random():
				deadIndivIndex = potientialDeadIndiv
				break
		return deadIndivIndex
	# If selection is not on, then just pick a random person
	else: 
		return int(random.random() * populationSize) 
		
##########################################################################################
#	Name: replaceDeadWithOffspring
#
#	Assumptions: It is assumed that the likelihood of getting a rho gene from either
#		     in the case of recombination is .5; it is also assumed that sections
#		     of the lDel contributed from both parents can differ in size. In the case
#		     of no recombination, it is assumed that the person who died cannot be
#                    the parent of the offspring replacing them.
#
#	Purpose: This method exists so that a dead individual can be replaced with an 
#	         offspring. There are two cases: recombination is off or it is on. If it
#		 is off, then the offspring gets copies of each of the genes of another 
#		 individual in the population who is randomly chosen. If it is on, then 
#		 the offspring gets a combination of different genes from each parent and 
# 		 a fitness reflecting the difference in gene.
#			 
#
#	Arguments: deadIndex - the index in the population list of the individual who 
#			       died. This will be the index of the offspring.
#		   population - a list containing the individuals of the population and all of
#				their genetic information
#		   recombination - a boolean corresponding to whether or not recombination
#				   is allowed when developing offspring.
#
#	Returns: Nothing
#
##########################################################################################
def replaceDeadWithOffspring(deadIndex, recombination, population):
	populationSize = len(population)
	if recombination:
		# Pick two parents
		mateOne = deadIndex
		while mateOne == deadIndex:
			mateOne = int(random.random() * populationSize)
		mateTwo = deadIndex
		while mateTwo == deadIndex:
			mateTwo = int(random.random() * populationSize)	
		
		# Randomly assign one of the rho values from the parents to the individual
		if random.random() < .5:
			population[deadIndex][0] = population[mateOne][0]
		else:
			population[deadIndex][0] = population[mateTwo][0]
		
		# Pick a recombination site for every fifty loci in the not including index first
		# or last index for the lDel gene
		numberOflDelRecombSites = int(lDelGeneLength / 50)
		lDelRecombSites = []
		lDelRecombSites.append(0)
		for i in range (numberOflDelRecombSites):
			lDelRecombSites.append(random.randint(0,lDelGeneLength - 1))
		lDelRecombSites.append(lDelGeneLength)
		lDelRecombSites.sort()
		
		# Splice lDel from recombination site to recombination site, alternating between
		# parent one and parent two
		population[deadIndex][1] = []
		for i in range (0, len(lDelRecombSites) - 1):
			if i % 2 == 1:
				population[deadIndex][1] += population[mateTwo][1]\
				                            [lDelRecombSites[i]:lDelRecombSites[i + 1]] 
			else:
				population[deadIndex][1] +=  population[mateOne][1]\
				                             [lDelRecombSites[i]:lDelRecombSites[i + 1]]
				                      
		# Pick recombination sites for alpha
		alphaLength = len(population[deadIndex][3])
		alphaRecombSites = 2 
		alphaDividers = []
		alphaDividers.append(0)
		for i in range (alphaRecombSites):
			alphaDividers.append(random.randint(0, alphaLength - 1))
		alphaDividers.append(alphaLength)
		alphaDividers.sort()
		
		# Splice together sections of alphas from both parents
		population[deadIndex][3] = []
		for i in range (0, len(alphaDividers) - 1):
			if i % 2 == 1:
				population[deadIndex][3] += population[mateTwo][3]\
				                            [alphaDividers[i]:alphaDividers[i + 1]] 
			else:
				population[deadIndex][3] += population[mateOne][3]\
							    [alphaDividers[i]:alphaDividers[i + 1]]
											
		# Pick recombination sites for beta
		betaLength = len(population[deadIndex][4])
		betaRecombSites = 2
		betaDividers = []
		betaDividers.append(0)
		for i in range (betaRecombSites):
			betaDividers.append(random.randint(0, betaLength - 1))
		betaDividers.append(betaLength)
		betaDividers.sort()
		
		# Splice together beta genes
		population[deadIndex][4] = []
		for i in range (0, len(betaDividers) - 1):
			if i % 2 == 1:
				population[deadIndex][4] += population[mateTwo][4]\
				                            [betaDividers[i]:betaDividers[i + 1]] 
			else:
				population[deadIndex][4] +=  population[mateOne][4]\
							     [betaDividers[i]:betaDividers[i + 1]]
		
	else:
		# Pick an individual to be the parent who isn't the person who just died
		mateOne = deadIndex
		while mateOne == deadIndex:
			mateOne = int(random.random() * populationSize)
			
		# Assigne the give the offspring the rho of the parent
		population[deadIndex][0] = population[mateOne][0]
		
		# Copy each lDel loci from the parent to the lDel gene of the offspring.
		# (this is done manually to avoid multiple pointers to one gene)
		population[deadIndex][1] = []
		for i in range (len(population[mateOne][1])):
			population[deadIndex][1].append(population[mateOne][1][i])
		
		# Copy the alpha gene from the parent to the offspring
		population[deadIndex][3] = []
		for i in range(len(population[mateOne][3])):
			population[deadIndex][3].append(population[mateOne][3][i])
			
		# Copy the bea gene from the parent to the offspring
		population[deadIndex][4] = []
		for i in range(len(population[mateOne][4])):
			population[deadIndex][4].append(population[mateOne][4][i])

				    
##########################################################################################
#	Name: mutateIndividual
#
#	Assumptions: Beta mutations only occur when an lDel mutation has already occured.
#
#	Purpose: This method checks to see if there's a rho, lDel, alpha, or beta
#		 mutation according to the mutation rates giveen by the user. Let k
# 		 be the number of beta locus. If the index of an lDel mutation is less
#		 or equal to k, then the beta gene is also mutated in the beta gene. If
#		 the lDel mutation loci is greater than k, then only that lDel loci is
#		 changed.
#
#	Arguments: population - a reference to the array containing the genetic 
#				information of each inividual
#
#		   plDelmutation - the probability of there being at least one mutation 
#				   in the given lDel gene
#
#		   pRhoMutation - the porbability of a mutation of the read through 
#				  rho
#
#		   pDelToNonDel - the probability of a given lDel loci to go from 
#				  deleterious to benign
#
#		   mutantIndex - the index in the population array that holds the geneitc
#				 information on the individual to be me mutated
#
#		   pAlphaMutations - the mutation rate for a given alplha gene 
#
#	Returns: Nothing
#
##########################################################################################
def mutateIndividual(mutantIndex, population, pRhoMutation, plDelMutation, pAlphaMutation, pBetaMutation,
		     pDelToNonDel, pNonDelToDel):
	if random.random() < pRhoMutation:
		population[mutantIndex][0] *= 10**random.gauss(0, .2)

	if random.random() < plDelMutation:
		# Pick the locus to change
		betaLength = len(population[mutantIndex][4])
		lDelLength = len(population[mutantIndex][1])
		changeLoci = random.randint(0, betaLength + lDelLength)
		
		# The lDel locus has a corresponding beta locus
		if changeLoci < betaLength:
			mutationOccured = 0
			# Loci to be changed is a 1
			if population[mutantIndex][1][changeLoci]:
				if random.random() < pDelToNonDel:
					population[mutantIndex][1][changeLoci] = 0
					mutationOccured += 1
			# Loci to ba changed is a 0
			else:			
				if random.random() < pNonDelToDel:
					population[mutantIndex][1][changeLoci] = 1
					mutationOccured += 1
					
			# THIS DOES BETA MUTATION
			# Only mutat a beta if an lDel was changed
			#if mutationOccured:
				# Add a number drawn out of a normal distribution
				#population[mutantIndex][4][changeLoci] += \
                        	#random.gauss(-1 * population[mutantIndex][4][changeLoci] / 50.0, 
				#population[mutantIndex][4][changeLoci]/ \
				#len(population[mutantIndex][4]))
		else:
			# Change an lDel locus in the same way as above but mod
			# the loci by the length in case it's >= 500
			if population[mutantIndex][1][changeLoci % lDelLength]:
				if random.random() < pDelToNonDel:
					population[mutantIndex][1][changeLoci % lDelLength] = 0
			else:
				if random.random() < pNonDelToDel:
                              		population[mutantIndex][1][changeLoci % lDelLength] = 1
			
	# THIF PART DOESN ALPHA MUTATION
	#if random.random() < pAlphaMutation:
		# Pick a an alpha locus to change
		#alphaLength = len(population[mutantIndex][3])
		#changeLoci = random.randint(0, alphaLength - 1)
	
		# Add another number drawn from a normal distribution
		#population[mutantIndex][3][changeLoci] += \
		#random.gauss(-1 * population[mutantIndex][3][changeLoci] / 50.0,
		#population[mutantIndex][3][changeLoci]/ \
            	#len(population[mutantIndex][3]))
