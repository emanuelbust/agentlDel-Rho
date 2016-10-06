import os, sys, random, math, time, pickle, numpy
import setup, cycle, output, terminate


# Output the time the program started
sys.stdout.write("[" + time.asctime() + "]: ")
sys.stdout.write("Program started...\n")
sys.stdout.flush()


# Most of these are taken from the 2011 paper
selection = True
environment = True
pNonDelToDel = .4
pDelToNonDel = .1 
plDelLociMutation = 10**-8
delLociFitnessCost = 20
pRhoMutation = 10**-5
nucleotidesPerlDel = 30
alphaGeneLength = 10
betaGeneLength = alphaGeneLength
pAlphaMutation = .0006 
pBetaMutation = .0006
pCooption = (1.0 + 7.0/9 + 7.0/9) * plDelLociMutation
envOptChangePerGeneration = 2000 
envOpt = 0.0

continuation = None
lDelGeneLength = -1
coopt = 0
replicateNumber = -1
# Population is to be intiailized
if len(sys.argv) == 9:
	# All of the parameters that are set via the command line
	recombination = int(sys.argv[1])
	os.chdir(sys.argv[2])
	populationSize = int(sys.argv[3])
	lDelGeneLength = int(sys.argv[4])
	generations = int(sys.argv[5])
	pOneLoci = float(sys.argv[6])
	pZeroLoci = 1 - pOneLoci
	coopt = int(sys.argv[7])
	replicateNumber = int(sys.argv[8])
	
	# Call for the population to be initialized
	population = setup.initializePopulation(populationSize, pOneLoci, 
		     delLociFitnessCost, lDelGeneLength, pNonDelToDel, 
		     pDelToNonDel, plDelLociMutation, alphaGeneLength, 
		     betaGeneLength)

	print("Recombinaton: ", recombination, " | Directory: ", sys.argv[2], 
	" | Population Size: ", populationSize, " | Number of lDels: ", lDelGeneLength, 
	" | Generations: ", generations, " | Initial lDel Frequency: ", pOneLoci, 
	" | Replication: ", replicateNumber)

	continuation = False
	
# A population is to be imported		
elif len(sys.argv) == 7:
	# Read in the population from a binary file
	with open(sys.argv[4], "rb") as inFile:
		population = pickle.load(inFile)

	recombination = int(sys.argv[1])
	os.chdir(sys.argv[2])
	populationSize = len(population)
	lDelGeneLength = len(population[0][1])
	generations = int(sys.argv[3])	
	coopt = int(sys.argv[5])
	replicateNumber = int(sys.argv[6])
	
	print("Recombination: ", recombination, " | Directory: ", sys.argv[2], 
	" | Population Size: ", populationSize, " | Number of lDels", 
	lDelGeneLength, " | Generations: ", generations, "| Replication: ",
	replicationNumber)

	continuation = True

# Incorrect arguements are given
else:
	print("You must have 8 for a new run or 7 to continue a run.")
	print("You have %d." % (len(sys.argv) - 1))
	exit(1)
	
# How man times statstics are to be recorded to text files
NUMBER_OF_PROGRESS_WRITES = 100
NUMBER_OF_RESULT_WRITES = generations
NUMBER_OF_POPULATION_WRITES = 100
NUMBER_OF_LDEL_COUNT_WRITES = 5
plDelMutation = plDelLociMutation * nucleotidesPerlDel * lDelGeneLength

# Output the time the population is done being initialized 
sys.stdout.write("[" + time.asctime() + "]: ")
sys.stdout.write("Population initialized...\n")
sys.stdout.flush()

# Create a stamp to distinguish this run 
stamp = str(time.time()) + "_rec=" + str(recombination) + \
"_n=" + str(populationSize) + "_lTot=" + str(lDelGeneLength) + \
"_initLdelFreq=" + str(pOneLoci) + "_cooption=" + str(coopt) + \
"_replicateNumber=" + str(replicateNumber) 

# Open and write the columns for the results text file
results = open(stamp + ".results" , "w")
columns = ["MeanFitness", "MeanFitnessVariance", "MeanlDels", "MeanlDelsVariance",
           "MeanRho", "MeanRhoVariance", "MeanAlpha", "MeanAlphaVariance", 
           "MeanBeta", "MeanBetaVariance", "EnvOpt", "meanEnvScore", 
	   "meanEnvScoreVariance"]
results.write("\t".join(columns) + "\n")
lDelOutName = stamp + ".out"
lDelOut = open(lDelOutName, "w")
lastPopFile = -1

# Replace n individuals each generation
replacements = populationSize * generations
replacementNumber = 1
while replacementNumber <= replacements:
	# Shift the environment after a given number of replacements
	if ((float(replacementNumber) / populationSize) % envOptChangePerGeneration) == 0 \
	    and (float(replacementNumber) / populationSize) != 0 and environment:
		envOpt = cycle.envShift(envOpt)
	
	# Kill someone
	deadIndex = cycle.pickDeadIndiv(selection, population, 2)	
	
	# Make a baby
	cycle.replaceDeadWithOffspring(deadIndex, recombination, population)
	
	# Mutate the baby
	cycle.mutateIndividual(deadIndex, population, pRhoMutation, plDelMutation, pAlphaMutation, pBetaMutation,
			 pDelToNonDel, pNonDelToDel)
	if coopt:
		cycle.cooption(deadIndex, pCooption, 3, 4, population)
	
	# Calculate the baby's fitness
	population[deadIndex][2] = setup.getFitness(delLociFitnessCost, population[deadIndex][1], 
	                                      pNonDelToDel, pDelToNonDel, plDelLociMutation, 
	                                      population[deadIndex][0], 
	                                      population[deadIndex][3], 
	                                      population[deadIndex][4], envOpt)
		
	# Records numbers decribing the population
	if replacementNumber % (populationSize * generations / NUMBER_OF_RESULT_WRITES) == 0:
		output.recordStatistics(results, population, envOpt)
		
	# Report progress
	if replacementNumber % (populationSize * generations / NUMBER_OF_PROGRESS_WRITES) == 0:
		sys.stdout.write("[" + time.asctime() + "]: ")
		sys.stdout.write(str(int((float(replacementNumber) / \
		                 (populationSize * generations) * 100))) + \
		                 "% complete\n")
		sys.stdout.flush()
	
	# Write the population to a file
	if replacementNumber % (populationSize * generations / NUMBER_OF_POPULATION_WRITES) == 0:
		# Save the population as a binary file
		if lastPopFile != -1:
			os.system("rm " + lastPopFile)
		popFile = stamp + ".pop"
		with open(popFile, "wb") as storage:
			pickle.dump(population, storage)
		lastPopFile = popFile
		
	# Count everyone's lDels	
	if replacementNumber % (populationSize * generations / NUMBER_OF_LDEL_COUNT_WRITES) == 0:
		output.outputlDelCount(population, 1, lDelOut)
	'''
	if replacementNumber > 2 * populationSize:
		if terminate.isDoneLax(.01, stamp + ".results"):
			results.write("LAX DONE\n")
		if terminate.isDoneStrict(.01, stamp + ".results"):
			results.write("STRICT DONE\n")	
	'''	

	replacementNumber += 1

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
#os.system("python ~emanuelb/agentSim/row2col.py " + os.getcwd()+ "/"  + lDelOutName)
