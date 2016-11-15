import numpy, setup

##########################################################################################
#	Name: recordStatistics
#
#	Assumptions: This assumes that the population list's entries will each have 5
#		     sub entries. The 5 sub entries are as following: the rho value of 
#		     the individual as a float, the lDel gene of the individual as a list of
#		     ones and zeroes, the fitness of the individual, the alpha gene of the 
#		     individual, and the beta gene of the individual. It is assumed that they
#	             appear in the order of the last sentences
#
#	Purpose: This function exists for the purpose finding the mean fitness, number of 
#	 	 ones in an lDel gene, and mean rho of an individual in a population. 
#		 Furthermore, the variance of each one of the aforementioned means is 
#		 calculated. Once these statistics are calculated, they are separated by
#		 tabs and written as a line and are written to the outFile object. 
#
#	Arguments: - outFile is file object representing the file where statistics will be 
#		     written to.
#		   - population is a list representing the population that the statistics
#		     will be about
#
#	Returns: Nothing
#
##########################################################################################
def recordStatistics(outFile, population, envOptimum):
	# The statistics to be recorded are initialized to 0
	populationSize = len(population)
	meanFitness = 0.0
	meanFitnessVariance = 0.0
	meanRho = 0.0
	meanlDelLoci = 0.0
	meanRhoVariance = 0.0
	meanlDelLociVariance = 0.0
	meanAlpha = 0.0
	meanAlphaVariance = 0.0
	meanBeta = 0.0
	meanBetaVariance = 0.0
	meanEnvScore = 0.0
	meanEnvScoreVariance = 0.0
	thirdMoment = 0.0
		
	# Iterate through the population once to find the mean
	for i in range (populationSize):
		meanRho += population[i][0] / float(populationSize)
		meanlDelLoci += numpy.sum(population[i][1]) / float(populationSize)
		meanFitness += population[i][2] / float(populationSize)
		meanAlpha += numpy.sum(population[i][3]) / float(populationSize)
		meanBeta += numpy.sum(population[i][4]) / float(populationSize)
		meanEnvScore += setup.getEnvScore(population[i][3], population[i][1],
					    population[i][4], population[i][0]) / \
				float(populationSize)
		
	# Iterate another time to find the variance using the means
	for i in range (populationSize):
		meanRhoVariance += (population[i][0] - meanRho)**2 / float(populationSize)
		meanlDelLociVariance += (numpy.sum(population[i][1]) - meanlDelLoci)**2 / \
		                        float(populationSize)
		meanFitnessVariance += (population[i][2] - meanFitness)**2 / float(populationSize)
		meanAlphaVariance += (numpy.sum(population[i][3]) - meanAlpha)**2 / \
		                          float(populationSize)
		meanBetaVariance += (numpy.sum(population[i][4]) - meanBeta)**2 / \
		                          float(populationSize)
		meanEnvScoreVariance += (setup.getEnvScore(population[i][3], population[i][1],
                                                     population[i][4], population[i][0]) - \
					 meanEnvScore) / float(populationSize)
	
	# Use the variance and mean to find the third moment
	for i in range(populationSize):
		thirdMoment += (population[i][2] - meanFitness)**(3) / meanFitnessVariance**(1.5)

	# Write the numbers into a text using the join method. 
	numbers = [meanFitness, meanFitnessVariance, meanlDelLoci, meanlDelLociVariance, 
	           meanRho, meanRhoVariance, meanAlpha, meanAlphaVariance, meanBeta, 
	           meanBetaVariance, envOptimum, meanEnvScore, meanEnvScoreVariance,
		   thirdMoment]
	data = [str(number) for number in numbers]
	outFile.write("\t".join(data) + "\n")
	outFile.flush()
	
##########################################################################################
#	Name: outputlDelCount
#
#	Assumptions: It is assumed that the individuals in a population are stored in a 
#		     list. Furthermore, it is assumed that each entry in the population
#		     list is another list containing genes of a given individual. 
#	             Finally, it is assumed that an lDel gene is a list and each 1 in the
#		     lDel gene is a loci that is deleterious.
#
#	Purpose: The purpose of this method is to output the lDel count to a file. This
#		 file will then be used to make a density plot in R. To accomplish this,
#		 each individual's lDel gene will be counted for lDels and then written
#		 to the desired output file on a new line. To be able to distinguish 
#		 between lDels of on group and another in the output file, a blank line
#		 will be used as a divider in the file, hence the new line write at the
#		 end of the function.
#
#	Arguments: - population is a list containing the lDels of each idividual
#		   - lDelIndex is the index of the sub list where the lDel of an indivudal
#		     can be accessed. I.e. population[individual][lDelIndex] gives access
#		     to individual's lDel gene
#		   - outfile is the file object representing the text file that the lDel
#		     measurements will be written to.
#
#	Returns: Nothing
#
##########################################################################################
def outputlDelCount(population, lDelIndex, outfile):
	for i in range (len(population)):
		outfile.write(str(numpy.sum(population[i][lDelIndex])) + " ")
	# The line below ensures that each group of ldel outputs is 
	# seperated by a new line
	outfile.write("\n")
