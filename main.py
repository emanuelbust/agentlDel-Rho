import os, sys, random, math, time, pickle
import setup, cycle, output


# Output the time the program started
sys.stdout.write("[" + time.asctime() + "]: ")
sys.stdout.write("Program started...\n")
sys.stdout.flush()


# All of the other paramaters. Most of these are taken from the 2011 paper
selection = True
environment = True
pNonDelToDel = .4
pDelToNonDel = .1 
plDelLociMutation = 10**-8
delLociFitnessCost = 7
pRhoMutation = 10**-5
nucleotidesPerlDel = 30
alphaGeneLength = 10
betaGeneLength = alphaGeneLength
pAlphaMutation = .001 
pBetaMutation = .001
pCooption = (1 + 7.0/9 + 7.0/9) * plDelLociMutation
envOptChangePerGeneration = 2000 
envOpt = 0.0

continuation = None
lDelGeneLength = -1
# Population is to be intiailized
if len(sys.argv) == 7:
	# All of the parameters that are set via the command line
	recombination = int(sys.argv[1])
	os.chdir(sys.argv[2])
	populationSize = int(sys.argv[3])
	lDelGeneLength = int(sys.argv[4])
	generations = int(sys.argv[5])
	pOneLoci = float(sys.argv[6])
	pZeroLoci = 1 - pOneLoci

	# Call for the population to be initialized
	population = setup.initializePopulation(populationSize, pOneLoci, 
		     delLociFitnessCost, lDelGeneLength, pNonDelToDel, 
		     pDelToNonDel, plDelLociMutation, alphaGeneLength, 
		     betaGeneLength)

	print("Recomb", recombination, " | Directory", sys.argv[2], 
	" | PopSize", populationSize, " | # of lDels", lDelGeneLength, 
	" | Gens", generations, " | pOne", pOneLoci, " |")

	continuation = False
	
# A population is to be imported		
elif len(sys.argv) == 5:
	# Read in the population from a binary file
	with open(sys.argv[4], "rb") as inFile:
		population = pickle.load(inFile)

	recombination = int(sys.argv[1])
	os.chdir(sys.argv[2])
	populationSize = len(population)
	lDelGeneLength = len(population[0][1])
	generations = int(sys.argv[3])	
		
	print("Recomb", recombination, " | Directory", sys.argv[2], 
	" | PopSize", populationSize, " | # of lDels", lDelGeneLength, 
	" | Gens", generations, " |")

	continuation = True

# Incorrect arguements are given
else:
	print("You must have 6 for a new run or 4 to continue a run.")
	print("You have %d." % (len(sys.argv) - 1))
	exit(1)
	
# How man times statstics are to be recorded to text files
NUMBER_OF_RESULT_WRITES = 1000
NUMBER_OF_LDEL_COUNT_WRITES = 5
plDelMutation = plDelLociMutation * nucleotidesPerlDel * lDelGeneLength

# Output the time the population is done being initialized 
sys.stdout.write("[" + time.asctime() + "]: ")
sys.stdout.write("Population initialized...\n")
sys.stdout.flush() 
# Open and write the columns for the results text file
results = open("results.txt", "w")
columns = ["MeanFitness", "MeanFitnessVariance", "MeanlDels", "MeanlDelsVariance",
           "MeanRho", "MeanRhoVariance", "MeanAlpha", "MeanAlphaVariance", 
           "MeanBeta", "MeanBetaVariance", "EnvOpt", "meanEnvScore", 
	   "meanEnvScoreVariance"]
results.write("\t".join(columns) + "\n")
lDelOutName = "lDels.txt"
lDelOut = open(lDelOutName, "w")


replacements = populationSize * generations
replacementNumber = 0
# Pick an individual to die, produce an offspring, mutates the offspring and report
# statstics on the population. One generation is N deaths and replacements
while replacementNumber < replacements:
	
	# The two conditions ensure that shift happen as often as they should 
	# and don't happen immediately
	if ((float(replacementNumber) / populationSize) % envOptChangePerGeneration) == 0 \
	    and (float(replacementNumber) / populationSize) != 0 and environment:
		envOpt = cycle.envShift(envOpt)
	

	deadIndex = cycle.pickDeadIndiv(selection, population, 2)	
	
	cycle.replaceDeadWithOffspring(deadIndex, recombination, population)
	
	cycle.mutateIndividual(deadIndex, population, pRhoMutation, plDelMutation, pAlphaMutation, pBetaMutation,
			 pDelToNonDel, pNonDelToDel)

	cycle.cooption(deadIndex, pCooption, 3, 4, population)
	
	population[deadIndex][2] = setup.getFitness(delLociFitnessCost, population[deadIndex][1], 
	                                      pNonDelToDel, pDelToNonDel, plDelLociMutation, 
	                                      population[deadIndex][0], 
	                                      population[deadIndex][3], 
	                                      population[deadIndex][4], envOpt)
		
	# Records numbers decribing the population NUMBER_RESULT_WRITES times
	if replacementNumber % (populationSize * generations / NUMBER_OF_RESULT_WRITES) == 0:
		output.recordStatistics(results, population, envOpt)
		
				############
				# BUG HERE #
				############
	# Report the progress of the script at each percent
	if replacementNumber % (populationSize * generations / 100) == 0:
		sys.stdout.write("[" + time.asctime() + "]: ")
		sys.stdout.write(str(int(replacementNumber / \
		                 (populationSize * generations) * 100)) + \
		                 "% complete\n")
		sys.stdout.flush()
		
	# Records numbers decribing the population NUMBER_OF_LDEL_COUNT_WRITES times	
	if replacementNumber % (populationSize * generations / \
	   NUMBER_OF_LDEL_COUNT_WRITES) == 0:
		output.outputlDelCount(population, 1, lDelOut)

	replacementNumber += 1

# Save the population as a binary file
popFile = str(time.time()) + "_n" + str(populationSize) + "_ldels" + str(lDelGeneLength) + ".pop"
with open(popFile, "wb") as storage:
        pickle.dump(population, storage)

# Write the parameters to the end of the text file and close the file
if not continuation:
	parameters = "[populationSize = " + str(populationSize) + ", Generations = " + \
		     str(generations) + ", pOneLoci = " + str(pOneLoci) + \
		     ", Recombination = " + str(recombination) + ", Ldel length = " + \
		     str(lDelGeneLength) + "]" 

else:
 
	parameters = "[populationSize = " + str(populationSize) + ", Generations = " + \
		     str(generations) + ", Recombination = " + str(recombination) + \
		     ", Ldel length = " + str(lDelGeneLength) + "]" 

results.write(parameters + "\n")
results.write("Population file: " + popFile + "\n")
results.close()
lDelOut.close()

# Output the time that the program finishes
sys.stdout.write("[" + time.asctime() + "]: ")
sys.stdout.write("Done\n")
sys.stdout.write("Population file: " + popFile + "\n")
sys.stdout.flush()

# Make all of the rows in the lDel output into columns
os.system("python ~emanuelb/agentSim/row2col.py " + os.getcwd()+ "/"  + lDelOutName)
