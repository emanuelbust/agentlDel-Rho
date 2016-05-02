import os, sys, random, math, time, pickle
import setup, cycle, output


# Output the time the program started
sys.stdout.write("[" + time.asctime() + "]: ")
sys.stdout.write("Program started...\n")
sys.stdout.flush()

# All of the parameters that are set via the command line
recombination = int(sys.argv[1])
os.chdir(sys.argv[2])
populationSize = int(sys.argv[3])
lDelGeneLength = int(sys.argv[4])
generations = int(sys.argv[5])
pOneLoci = float(sys.argv[6])
pZeroLoci = 1 - pOneLoci

print("Recomb", recombination, "Directory", sys.argv[2], "PopSize", populationSize, "# of lDels", lDelGeneLength, "Gens", generations,"pOne", pOneLoci)

# All of the other paramaters. Most of these are taken from the 2011 paper
selection = True
pNonDelToDel = .4
pDelToNonDel = .1 
# per base mutation rate
plDelLociMutation = 10**-8
delLociFitnessCost = 7
# multiply below times 1000
pRhoMutation = 10**-5
nucleotidesPerlDel = 30
plDelMutation = plDelLociMutation * nucleotidesPerlDel * lDelGeneLength
alphaGeneLength = 10
betaGeneLength = alphaGeneLength
pAlphaMutation = plDelLociMutation * nucleotidesPerlDel * alphaGeneLength
pBetaMutation = plDelLociMutation * nucleotidesPerlDel * betaGeneLength
envOptChangeRate = 15000
envOpt = 0

# How man times statstics are to be recorded to text files
NUMBER_OF_RESULT_WRITES = 100
NUMBER_OF_LDEL_COUNT_WRITES = 5

# Call for the population to be initialized
population = setup.initializePopulation(populationSize, pOneLoci, delLociFitnessCost, 
				  lDelGeneLength, pNonDelToDel, pDelToNonDel, 
				  plDelLociMutation, alphaGeneLength, betaGeneLength)

# Output the time the population is done being initialized 
sys.stdout.write("[" + time.asctime() + "]: ")
sys.stdout.write("Population initialized...\n")
sys.stdout.flush() 
# Open and write the columns for the results text file
results = open("results.txt", "w")
columns = ["MeanFitness", "MeanFitnessVariance", "MeanlDels", "MeanlDelsVariance",
           "MeanRho", "MeanRhoVariance", "MeanAlpha", "MeanAlphaVariance", 
           "MeanBeta", "MeanBetaVariance", "EnvOpt"]
results.write("\t".join(columns) + "\n")
lDelOut = open("lDelOutput.txt", "w")
	
# Pick an individual to die, produce an offspring, mutates the offspring and report
# statstics on the population. One generation is N deaths and replacements
for replacementNumber in range (populationSize * generations):
	#if replacementNumber % envOptChangeRate == 0 and replacementNumber != 0:
		# 5 and 4 are paramters that are found to work the best. They're not from an assumption
		#envOpt += random.gauss((-1 * envOpt) / float(5), 4)

	deadIndex = cycle.pickDeadIndiv(selection, population, 2)	
	
	cycle.replaceDeadWithOffspring(deadIndex, recombination, population)
	
	cycle.mutateIndividual(deadIndex, population, pRhoMutation, plDelMutation, pAlphaMutation, pBetaMutation,
			 pDelToNonDel, pNonDelToDel)

	population[deadIndex][2] = setup.getFitness(delLociFitnessCost, population[deadIndex][1], 
	                                      pNonDelToDel, pDelToNonDel, plDelLociMutation, 
	                                      population[deadIndex][0], 
	                                      population[deadIndex][3], 
	                                      population[deadIndex][4], envOpt)
		
	# Records stats 100 times every run and report the progress of the script
	if replacementNumber % (populationSize * generations / NUMBER_OF_RESULT_WRITES) == 0:
		output.recordStatistics(results, population)

		sys.stdout.write("[" + time.asctime() + "]: ")
		sys.stdout.write(str(int(replacementNumber / \
		                 (populationSize * generations) * NUMBER_OF_RESULT_WRITES)) + \
		                 "% complete\n")
		sys.stdout.flush()
		
	if replacementNumber % (populationSize * generations / \
	   NUMBER_OF_LDEL_COUNT_WRITES) == 0:
		output.outputlDelCount(population, 1, lDelOut)

# Write the parameters to the end of the text file and close the file
parameters = "[populationSize = " +  str(populationSize) +  ", Generations = " + \
             str(generations) +  ", pOneLoci = " +  str(pOneLoci) + \
             ", Recombination = "  + str(recombination) + ", Ldel length = " + \
             str(lDelGeneLength) +"]" 
results.write(parameters)
results.write("\n")
results.close()
lDelOut.close();

# Pickle the population into a binary file
with open("frozenPopulation.bin", "wb") as storage:
	pickle.dump(population, storage)

# Output the time that the program finishes
sys.stdout.write("[" + time.asctime() + "]: ")
sys.stdout.write("Program done\n")
sys.stdout.flush()
