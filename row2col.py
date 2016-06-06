##########################################################################################
#       Name: row2col.py
#
#       Assumptions: It is assumed that the first command line arguement given to this 
#		     file will be an absolute path to a file. It is also assumed that 
#		     there are the number of rows and columsn are equal. 
#
#       Purpose: This script takes rows in a given file and writes them as columns in a
#		 new file. Essentially is takes a square matrix from one file and writes
#		 its transpose to a new file with a .out extension.
#
#       Arguments: The first arguement is an absolute file path
#
#       Output: A file with the same name as the inputed file---except with a .out 
#		extension.
#
##########################################################################################
import sys, os

with open(sys.argv[1], "r") as inFile:
	with open(sys.argv[1][:-3] + "out", "w") as outFile: 
		rows = []
		for row in inFile:
			entries = row.split(' ')
			rows.append(entries)
		
		numRows = len(rows[0])
		numCols = len(rows)	
		
		for i in range(numRows):
			row = []
			for j in range(numCols):
				row.append(rows[j][i])
			outFile.write(("\t".join(row)).replace("\n", "") + "\n")
