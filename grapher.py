import matplotlib.pyplot as plt
import sys

##########################################################################################
#	Name: parse
#
#	Assumptions: None
#
#	Purpose: parse takes a text file and separates all of the entries in the line. 
#			 An entry is a string of text in between two delimiters, the beginning of a 
#			 line and the delimiter, or the delimiter and the end of a line. Entries on 
#			 the same line are put in a list, and the entries from each line is put in a 
#			 larger list in the order that the lines appear.
#
#	Arguments: filePath - a string that is the path of the file to be read
#
#			   delimiter - 	the string that serves as a delimiter for entries		
#
#	Returns: a list whose elements correspond to the entries in each line of the file.
#			 I.e. lines[i][j] gives the jth entry in the ith line of the text file.
#
##########################################################################################
def parse(filePath, delimiter):
	# Read the the input file and the delimiter
	inFile = open(filePath, "r")
	
	# Split the columns
	lines = []
	for line in inFile:
		lines.append(line.split(delimiter))

	# Close the file and return
	inFile.close()
	return lines

##########################################################################################
#	Name: trimMatrix
#
#	Assumptions: trimMatrix assumes that the matrix has at least one row. If the matrix
#				 matrix does not have at least one row, this function will try to access
#				 a list entry that doesn't exist.
#
#	Purpose: trimMatrix makes takes a matrix and makes it rectangular. The dimension of 
#			 the resulting matrix has the width of the first row of the original matrix
#			 and a height of the number of consecutive rows that have the same width of
#			 the header row (including the header row).
#
#	Arguments: matrix - the matrix (2D array) to be trimmed
#
#	Returns: A matrix containing that is a subset of the original matrix. The returned 
#			 matrix contains the first stretch of rows of the original sequence that have
#			 the same width. This matrix is returned in a dictionary along with its 
#			 dimensions.
#
##########################################################################################
def trimMatrix(matrix):
	# Make the width of the square matrix the width of the header
	desiredWidth = len(matrix[0])
	
	# Only keep rows that make the matrix a square
	rectangle = []
	for row in matrix:
		if len(row) == desiredWidth:
			rectangle.append(row)
		else:
			break
			
	# Return the matrix along with its dimensions
	return {"matrix": rectangle, "rows": len(rectangle), "columns": desiredWidth}

##########################################################################################
#	Name: getColumn
#
#	Assumptions: getColumn assumes that the matrix being subsetted is rectangular/
#
#	Purpose: getColumn returns a specified column.
#
#	Arguments: matrix - the matrix to be subsetted
#			
#			   col - the integer corresponding to which column will be returned
#
#	Returns: A list that has all of the entries in the specified column of the given
#			 matrix
#
##########################################################################################
def getColumn(matrix, col):
	desiredCol = []
	for row in matrix:
		desiredCol.append(row[col])
		
	return desiredCol

##########################################################################################
#	Name: graphColumn
#
#	Assumptions: graphColumn assumes that the first entry in any column is a label and the
#				 rest of the entries are numeric data.
#
#	Purpose: graphColumn gets a column from a matrix, graphs it, and then saves the file.
#			 This is done by extracting a column from a matrix and using its first entry
#			 as a title and using the rest of the entries as data. The plot is saved as
#			 the title name with a png extension.
#
#	Arguments: matrix - the matrix of data containing the column to be plotted
#
#			   col - the integer corresponding to which column will be plotted
#
#	Returns: Nothing
#
##########################################################################################
def graphColumn(matrix, col):
	# Get the data and label for the graph
	column = getColumn(matrix, col)
	title = column[0].replace("\n", "")
	
	# Cast data to floats and remove and new line characters
	column = [float(num.replace("\n", "")) for num in column[1:]]
	
	# Plot the figure and save it
	plt.clf()
	plt.plot(column)
	plt.ylabel(title)
 	plt.savefig(title + ".png")
 	plt.clf()

