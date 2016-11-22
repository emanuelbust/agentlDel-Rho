import random

###############################################################################
#       Name: recombine
#
#       Assumption: None
#
#       Purpose: Recombine takes two lists and change the contents of a third 
#		 to corresponding chunks of the first two lists. This is done
#		 by copying over elements from either the first or second list
#		 to the third in a random fashion.
#
#       Arguments: recombined - the array that is made of chunks of the 
#				first two arrays#					
#
#		   listOne - the 1st list that contributes to the recombined
#			     list
#
#		   listTwo - the 2nd list that contributes to the recombined
#			     list
#
#       Returns: Nothing
#
###############################################################################
def recombine(sites, listOne, listTwo, dest, length):
	# Choose the indices where the segments will alternate
	indices = [0, length] + [int(random.random() * length) for i in range(sites)]
        indices.sort()
	
	# Build the new list. Take a section from listTwo first.
	for i in range(sites + 1):
		if random.random() < .5: 
			dest[indices[i]:indices[i + 1]] = listOne[indices[i]:indices[i + 1]]
		else:
			dest[indices[i]:indices[i + 1]] = listTwo[indices[i]:indices[i + 1]]

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
		if type(origList[i]) is int or type(origList[i]) is str or type(origList[i]) is bool or type(origList[i]) is float:
			newList.append(origList[i])
        # Recursively copy the item if it is a list
		elif type(origList[i]) is list:
			newList.append(myCopy(origList[i]))
        # Freak out if some other type of object is in the list
		else:
			print("Error: List contains invalid types")
			exit(1)

	return newList

##########################################################################################
#	Name: pickParents
#
#	Assumptions: It is assumed that the population has entries 0-populationSize - 1 
#
#	Purpose: pickParents exists for the sake of picking parents to produce an
#		 offspring. If selection is turned off, then two distinct individuals are
#		 chosen at random to be parents. If recombination is turned on, the first
#		 two distinct individuals with a fitness higher than a random number are
#		 chosen to the parents.
#
#	Arguments: selection - a boolean corresponding to whether or not selection is
#			       turned on
#
#		   population - a list containing the individuals of the population
#
#		   fitnessIndex - the index of an individuals fitness. 
#				  population[i][fitnessIndex] gives the fitness of 
#				  individual i
# 
#		   populationSize - the number of individuals contained the population.
#				    len(population) == populationSize 
#
#	Returns: A tuple containing two distinct integers in the interval
#		 [0, populationSize - 1. Each integer corresponds to the the entry of the
#		 population containg the genetic information of a parent. 
#
##########################################################################################
def pickParents(selection, population, fitnessIndex, populationSize):
	# If selection is on, discriminate against low fitness
	if selection:
		# Pick the first mom who has a fitness greater than a random number
		while True:
			mom = int(random.random() * populationSize)
			if population[mom][fitnessIndex] > random.random():
				break
		
		# Pick the first dad who has a fitness greater than a random number and isn't already the mom
		while True:
			dad = int(random.random() * populationSize)
			if population[dad][fitnessIndex] > random.random() and dad != mom:
				break

	# If slection is off, don't consider fitness
	else:
		# Pick a mom at random
		mom = int(random.random() * populationSize)
		
		# Pick the first dad who isn't already the mom
		dad = mom
		while dad == mom:
			dad = int(random.random() * populationSize)
	
	return (mom, dad)

##########################################################################################
#	Name: replaceDeadWithOffspring
#
#	Assumptions: It is assumed that the likelihood of getting a rho gene from either
#		     in the case of recombination is .5; it is also assumed that sections
#		     of the lDel contributed from both parents can differ in size. In the case
#		     of no recombination, it is assumed that the person who died cannot be
#                    the parent of the offspring replacing them.
#
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
def replaceDeadWithOffspring(deadIndex, recombination, population, populationSize, selection, 
			     fitnessIndex):
	if recombination:
		# Pick two parents
		parents = pickParents(selection, population, fitnessIndex, populationSize)
		mateOne = parents[0]
		mateTwo = parents[1]
		
		# Randomly assign one of the rho values from the parents to the individual
		population[deadIndex][0] = population[mateOne][0]
		if random.random() < .5:
			population[deadIndex][0] = population[mateTwo][0]
		
		# Recombine the parents and then give the offspring the result for 
		# the lDel, alpha, and beta gene of the osspring  
		recombine(5, population[mateOne][1], population[mateTwo][1], 
				  population[deadIndex][1], population[deadIndex][1].size)
	
		recombine(2, population[mateOne][3], population[mateTwo][3], 
				  population[deadIndex][3], population[deadIndex][3].size)
				
		recombine(2, population[mateOne][4], population[mateTwo][4], 
				  population[deadIndex][4], population[deadIndex][4].size)
		
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
		     pBetaMutation, pDelToNonDel, pNonDelToDel, pCooption):
	if random.random() < pRhoMutation:
		population[mutantIndex][0] *= 10**random.gauss(0, .2)

	if random.random() < plDelMutation:
		# Pick the locus to change
		betaLength = population[mutantIndex][4].size
		lDelLength = population[mutantIndex][1].size
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
					
			# Only mutat a beta if an lDel was changed
			if mutationOccured:
				# Add a number drawn out of a normal distribution
				population[mutantIndex][4][changeLoci] += \
				random.gauss(-1 * population[mutantIndex][4][changeLoci] / 50.0,
                		10 / float(betaLength))
				
				# Check for a cooption mutation
				if random.random() < pCooption: 
					population[mutantIndex][3][changeLoci] += population[mutantIndex][4][changeLoci]
					population[mutantIndex][4][changeLoci] = random.gauss(0.0, 16.0 / (9.0/25.0))
	
		else:
			# Change an lDel locus in the same way as above
			if population[mutantIndex][1][changeLoci]:
				if random.random() < pDelToNonDel:
					population[mutantIndex][1][changeLoci] = 0
			else:
				if random.random() < pNonDelToDel:
                              		population[mutantIndex][1][changeLoci] = 1
	
	# THIF PART DOES ALPHA MUTATION
	if random.random() < pAlphaMutation:		
		# Pick a an alpha locus to change
		alphaLength = population[mutantIndex][3].size
		changeLoci = random.randint(0, alphaLength - 1)
			
		# Add another number drawn from a normal distribution
		delta = random.gauss(-1 * population[mutantIndex][3][changeLoci] / 50.0,
                10 / float(alphaLength))	
	
		population[mutantIndex][3][changeLoci] += delta	
