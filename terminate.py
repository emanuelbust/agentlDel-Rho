import csv
import math

##########################################################################################
#       Name: getlDels
#
#	Purpose: Gets a  list of lDels from a text file
#
#       Assumptions: It is assumed that each column in the text file is tab delineated. 
#		     Furthermore, it is assumed that the lDel sequence is the 3rd column 
#		     in the text file.
#
#       Arguments: inputFile - a string that of the file containing the lDel sequence         
#
#       Returns: A list of floats corresponding the the lDels in the text file.
#
##########################################################################################
def getlDels(inputFile):
	# Get rid of the header and convert the numbers to floats
	return [float(i) for i in getCol(inputFile, "\t", 2)[1:]]
	
##########################################################################################
#       Name: getCol
#
#	Purpose: Creates a list of an arbitrary column from a text file given a delimiter
#
#       Assumptions: The first column is the 0th column. 
#
#       Arguments: inputFile - a string that of the file containing the lDel sequence
#		   col - the number of the column that is desired
#		   delim - a string that is the delimiter of columns in the input file
#
#       Returns: A list of floats corresponding the entires of desired column.
#
##########################################################################################
def getCol(inputFile, delim, col):
	# Read the results text file as a matrix including column names
	with open(inputFile, "r") as inFile:
		reader = csv.reader(inFile, delimiter=delim)
		
		# Subset the the matrix corresponding to the text file by taking the col-th entry 
		# of each row
		desiredCol = []
		for row in reader:
			if len(row) > 5:
				desiredCol.append(row[col])
			
		# Return the column
		return desiredCol

##########################################################################################
#       Name: averageDiff
#
#	Purpose: To calculate the average distance between the last n adjacent entries in 
# 		 a list.
#
#       Assumptions: None
#
#       Arguments: values - a list of floats corresponding the the values whose average
#		 	    average difference is to be calculated
#		   n - the number of values whose average difference will be calculated
#
#       Returns: A float corresponding to the average difference of the last in values in
#		 the given list.
#
##########################################################################################
def averageDiff(values):
	# Check to see there are enough values to add
	if len(values) < 2:
		print("Not enough values.")
		exit(1)
	
	# Find the average difference of two adjacent values in the list
	sum = 0
	terms = 0
	position = len(values) - 1
	while position != 0:
		sum += values[position] - values[position - 1]
		terms += 1
		position -= 1
	
	return float(sum) / float(terms)

##########################################################################################
#       Name: isDoneLax
#
#	Purpose: Compares the average distance between the last n adjacent entries in a 
#		 a list to a theshold value. If the average distance is less than the 
#		 threshold, the test for convergence passes. If not, the test fails. This
#		 test is meant to be a lax version of isDoneStrict. 
#
#       Assumptions: It is assumed that the third row of the text file is the sequence
#		     is being tested for convergence. 
#
#       Arguments: threshold - a float serving as the the upper bound for lists that 
#			       have converged
#		   inFile - the text file containing the sequences that is being tested
#			    for convergence
#		   n - the number of terms in the sequence that will be tested for 
#		       convergence
#
#       Returns: 1 if the sequence converges and 0 if the sequence does not.
#
##########################################################################################
def isDoneLax(threshold, inFile):
	# Check for a valid threshold
	if threshold < 0:
		print "Threshold values must be positive."
		exit(1)

	# Compare the sum to the threshold value
	lDels = getlDels(inFile)
	if abs(averageDiff(lDels)) < threshold:
		return True
	else:
		return False

##########################################################################################
#       Name: isDoneStrict
#
#	Purpose: Finds the greatest distance between any given two points in a sequence 
#		 and compares it to a threshold value. If the us less than a given 
#		 theshold value, then the sequence convergence. Otherwise, the sequence
#		 does not converge. This is meant to be a more strict test than isDoneLax.
#
#       Assumptions: It is assumed that the third row of the text file is the sequence
#		     being tested for convergence. 
#
#       Arguments: threshold - a float serving as the the upper bound for lists that 
#			       have converged
#		   inFile - the text file containing the sequences that is being tested
#			    for convergence
#		   n - the number of terms in the sequence that will be tested for 
#		       convergence        
#
#       Returns: 1 if the sequence converges and 0 if the sequence does not.
#
##########################################################################################	
def isDoneStrict(threshold, inFile):
	# Check for a valid threshold
	if threshold < 0:
		print "Threshold values must be positive."
		exit(1)

	# Find the range of the last n lDels
	lDels = getlDels(inFile)
	maxlDel = max(lDels)
	minlDel = min(lDels)
	range = abs(maxlDel - minlDel)
	
	# Compare the range to the threshold value 
	if abs(range / float(len(lDels) - 1)) < threshold:
		return True
	else:
		return False
