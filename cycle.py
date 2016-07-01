import os, sys, random, math, time, pickle

###############################################################################
#       Name: recombine
#
#       Assumption: It is assumed that the lists passed to the function do not
#		    not contain objects. If the list does contain objects, the
#		    function works correctly, but each list involved will have
#	            a reference to the same object (shallow copy).
#
#       Purpose: A given number of recombination sites are randomly chosen.
#		 Between each two recombination sites, a corresponding piece
#		 of list from listOne or listTwo is added. Which list is
#		 chosen to contribute is random.  
#
#       Arguements: - sites: an integer corresponding to the number of 
#		      indices where the new list alternates
#
#                   - listOne: a list containing information that will be 
#                              sectioned to form a new list
#
#                   - listTwo: another list containing information that will
#                              be section off to form a new list
#
#       Returns: A new list that is a recombined version of the two lists
#		 given.
#
###############################################################################
def recombine(sites, listOne, listTwo):
        # Check that the two lists are the same size
        if len(listOne) != len(listTwo):
                print "Lists are different size"
                exit(1)
        
        # Choose the indices where the segments will alternate
        length = len(listOne)
	indices = [0, len(listOne)]
	for i in range(sites):
		indices.append(random.randint(0, len(listOne)))
	indices.sort()		
	
        # Build the new list. Take a section from listTwo first.
        product = []
        for i in range(len(indices) - 1):
                if random.random() <= .5: 
                        product += listOne[indices[i]:indices[i + 1]]
                else:
                        product += listTwo[indices[i]:indices[i + 1]]

        # Check that the new list is the size of the two original lists
        if len(listOne) == len(listTwo) and len(listTwo) == len(product):
                return product
        else:
                print "There's a bug in the recombination code."
		exit(1)

################################################################################
#       Name: myCopy 
#
#       Assumptions: It is assumed that lists passed to this function will
#                    only integers, lists, or boolean types. 
#
#       Purpose: myCopy makes a copy of a list by recursively appending each
#                item in the original list to a new list.
#
#       Arguements: origList is a the list to be copied
#
#       Returns: Nothing if the list contains a type that is not one of the
#                three stated in the assumptions. If the list contanins valid
#                types, then a copy of the original list is returned.
################################################################################
def myCopy(origList):
        newList = []
        for i in range(len(origList)):
    		# Append if the item in the list is an integer, string, 
     		# or boolean
		if type(origList[i]) is int or \
		   type(origList[i]) is str or \
		   type(origList[i]) is bool or \
		   type(origList[i]) is float:
			newList.append(origList[i])
                # Recursively copy the item if it is a list
                elif type(origList[i]) is list:
                        newList.append(myCopy(origList[i]))
                # Freak out if some other type of object is in the list
                else:
                       	print "Error: List contains invalid types"
                        exit(1)

        return newList

###############################################################################
#       Name: cooption 
#
#       Assumptions: The distribution that the new beta has the assumptions 
#		     made the Rajon Masel 2011 paper. The equations come 
#		     specifically from the the supplemental instruction.
#
#       Purpose: This functions checks to see whether a cooption event
#		 happens. A cooption occurs when a stop codon is lost and the
#		 the coding sequence then permanently includes the cryptic 
#		 sequence. This manifests itself as an alpha trait being the
#		 sum of the alpha and beta and a new beta being chosen 
#		 randomly.
#
#       Arguements: indiv - an integer corresponding to the index holding 
#			    the genetic information of the individual who
#			    may experience cooption
#		    pCooption - a float corresponding to the probability of 
#				a cooption happening
#		    alphaIndex - an integer corresponding to the index of 
#				 of the alpha gene in an individual's entry
#				 in the population array
#		    betaIndex - an integer corresponding to the index of 
#                               of the beta gene in an individual's entry
#                               in the population array
#
#       Returns: Nothing
# 
###############################################################################
def cooption(indiv, pCooption, alphaIndex, betaIndex):
	if random.random() < pCooption:
		population[indiv][alphaIndex] += population[indiv][betaIndex]
		
		population[indiv][betaIndex][index] = \
		random.gauss(0.0, 16.0 / (9.0/25.0))
				
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
	lDelGeneLength = len(population[0][1])
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
		
		# Recombine the parents and then give the offspring the result for 
		# the lDel, alpha, and beta gene of the osspring  
		population[deadIndex][1] = recombine(5, population[mateOne][1], 
							population[mateTwo][1])
	
		population[deadIndex][3] = recombine(3, population[mateOne][3], 
							population[mateTwo][3])
				
		population[deadIndex][4] = recombine(3, population[mateOne][4], 
							population[mateTwo][4])
		
	else:
		# Pick an individual to be the parent who isn't the person who just died
		mateOne = deadIndex
		while mateOne == deadIndex:
			mateOne = int(random.random() * populationSize)
		# Copy the the parent into the offsprings slot recursively	
		population[deadIndex] = myCopy(population[mateOne])	
		    
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
def mutateIndividual(mutantIndex, population, pRhoMutation, plDelMutation, pAlphaMutation, 
		     pBetaMutation, pDelToNonDel, pNonDelToDel):
	if random.random() < pRhoMutation:
		population[mutantIndex][0] *= 10**random.gauss(0, .2)

	if random.random() < plDelMutation:
		# Pick the locus to change
		betaLength = len(population[mutantIndex][4])
		lDelLength = len(population[mutantIndex][1])
		changeLoci = random.randint(0, lDelLength - 1) 
		
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
			if mutationOccured:
				# Add a number drawn out of a normal distribution
				population[mutantIndex][4][changeLoci] += \
                        	random.gauss(-1 * population[mutantIndex][4][changeLoci] / 50.0, 
				population[mutantIndex][4][changeLoci]/ \
				len(population[mutantIndex][4]))
			
		else:
			# Change an lDel locus in the same way as above
			if population[mutantIndex][1][changeLoci]:
				if random.random() < pDelToNonDel:
					population[mutantIndex][1][changeLoci] = 0
			else:
				if random.random() < pNonDelToDel:
                              		population[mutantIndex][1][changeLoci] = 1
			
	# THIF PART DOESN ALPHA MUTATION
	if random.random() < pAlphaMutation:		
		# Pick a an alpha locus to change
		alphaLength = len(population[mutantIndex][3])
		changeLoci = random.randint(0, alphaLength - 1)
			
		# Add another number drawn from a normal distribution
		delta = random.gauss(-1 * population[mutantIndex][3][changeLoci] / 50.0,
                10 / float(alphaLength))	
	
		population[mutantIndex][3][changeLoci] += delta	

##########################################################################################
#       Name: envShift
#
#       Assumptions: None
#
#       Purpose: This method takes the current environmental optimum, increments it with 
#		 a number taken from a noram distriubtion a mean of the negative of the 
#		 current optimum and a standard deviation of four. The mean is that so
#		 that the optimum doesn't stray too far zero. The standard deviation is
#		 that just because that's what whoever tweaked it thought was the best 
#		 paramater 
#
#       Arguments: A float called currentOpt that is the current environmental optimum            
#
#       Returns: A float that is the incremented previous optimum
#
##########################################################################################
def envShift(currentOpt):
	# 5 and 4 are paramters that are found to work the best. They're not from an assumption
	currentOpt += random.gauss((-1 * currentOpt) / float(5), 4)
	return currentOpt
