import getopt
import math
import sys
import argparse

from TargetGrid import *

"""
TargetSearch main function
- Uses TargetGrid to search 2 proteins.
- Takes 2 PDB files, optional distance variable, optional output file NameError
- All atoms in PDB files are ran with TargetGrid.py
"""

##read both input PDB files
def parsePDB(PDB1,PDB2):
	f=open(PDB1,"r")
	lines1=f.readlines()
	f.close()
	
	##X,Y,Z coordinates found in lines[7],lines[8],lines[9] of PDB file
	for	x	in range(len(lines1)):
		#if lines1[x][0]	== 'A':
		lines1[x]	=	lines1[x].split()
		
	f=open(PDB2)
	lines2=f.readlines()
	f.close()
	
	for	x	in range(len(lines2)):
		#if lines2[x][0]	== 'A':
		lines2[x]	=	lines2[x].split()
	
	##coordinates are added to a list and then appended to make a list of points
	points1 = []
	l_index1 = [] 
	for x	in lines1:
		coord = []
		if x[0]	== "ATOM":
			coord.append(float(x[6]))
			coord.append(float(x[7]))
			coord.append(float(x[8]))
			points1.append(coord)
			l_index1.append(x)

	points2 = []
	l_index2 = [] 
	for x	in lines2:
		coord = []
		if x[0]	== "ATOM":
			coord.append(float(x[6]))
			coord.append(float(x[7]))
			coord.append(float(x[8]))
			points2.append(coord)
			l_index2.append(x)

	return points1, points2, l_index1, l_index2


##write into PDB format
def writePDB(input,output,l_index,outputFileName):
	##remove .pdb from name
	reformatInput = input.replace('.pdb','')
	
	##if no output file name given as command line parameter, set output name to "[input]_TargetSearch.pdb"
	if outputFileName == None:
		outputFileName = str(reformatInput) + "_TargetSearch.pdb"
	r = open(input,"r")
	lines = r.readlines()
	r.close()
	f=open(outputFileName,"w+")
	for i in l_index:
		f.write(lines[i])
#	for i in range(len(output)):
#		for j in lines:
#			if j[0] == 'ATOM':
#				if str(output[i]) == j[1]:
#					j[0] = j[0].ljust(6)
#					j[1] = j[1].rjust(5) + str('  ')
#					j[2] = j[2].ljust(4)
#					j[3] = j[3].ljust(3)
#					j[4] = j[4].rjust(2)
#					j[5] = j[5].rjust(4) + str('    ')
#					j[6] = str('%8.3f' % (float(j[6]))).rjust(8)
#					j[7] = str('%8.3f' % (float(j[7]))).rjust(8)
#					j[8] = str('%8.3f' % (float(j[8]))).rjust(8)
#					j[9] = j[9].rjust(6)
#					j[10]= j[9].rjust(18)
#					f.write(''.join(map(str,j))+"\n")
	f.close
	print("File written to " + outputFileName)

##compare distances
def compareDist(m1,	m2):
	d	=	0
	sim1 = []
	sim2 = []
	for	i	in range(len(m1)):
		for	j	in range(len(m2)):
			d1 = (m1[i][6]-m2[j][6])**2
			d2 = (m1[i][7]-m2[j][7])**2
			d3 = (m1[i][8]-m2[j][8])**2
			d	=	math.sqrt(d1+d2+d3)
			if d < 5:
				sim1.append(m1[i])
				sim2.append(m2[j])
	return sim1, sim2


##Parse arguments from command line
def parseArg():
	#parse output to take two inputs -i
	parser = argparse.ArgumentParser(description = 'Identify all points between two proteins that are within a certain distance of each other.')
	parser.add_argument('-i', nargs = 2, metavar = 'InputPDB', help = 'Input PDB file to be compared.')
	parser.add_argument('-d', nargs = '?', metavar = 'distance', type = int, help = 'Resolution for distance checking.')
	parser.add_argument('-o', nargs = '?', metavar = 'OutputPDB', help = 'Output PDB file name' )

	#parse list of points from inputs
	args = parser.parse_args()
	args = vars(args)
	i1,i2 = args['i'][0], args['i'][1]
	d = args['d']
	o = args['o']
	return i1,i2,d,o


##Removes duplicates and also parses out any atoms with # = -1 or 0
def removeDupe(output):
		output = list(dict.fromkeys(output))
		for i in output:
			if i == -1:
				output.remove(-1)
			if i == 0:
				output.remove(0)
		return output

def main():
	
	input1, input2,dist,outputFileName = parseArg()
	if dist == None:
		dist = 5
	
	##list of points from first input PDB or second input PDB
	##lines1, lines2 are original parsed lines from PDB files used for writing back to output files
	points1,points2,l_index1,l_index2 = parsePDB(input1,input2)
	
	##querySphere for every point in points1 on the target grid for points 2 -> vice versa
	output1 = []
	t2 = TargetGrid(points2)
	for i in range(len(points1)):
		x = points1[i][0]
		y = points1[i][1]
		z = points1[i][2]
		t2.querySphere(x,y,z,dist,output1)

	print(t2)
	output2 = []
	t1 = TargetGrid(points1)
	for i in range(len(points2)):
		x = points2[i][0]
		y = points2[i][1]
		z = points2[i][2]
		t1.querySphere(x,y,z,dist,output2)
	
	##remove duplicate points in output1 and output2
	output1 = removeDupe(output1)
	output2 = removeDupe(output2)
	
	#print(output1)
	#print(output2)
	
	#print(lines1)
	#print(lines2)
	
	##takes input names, output lists of points and lines from input file and writes PDB
	writePDB(input2,output1,l_index2, outputFileName)
	writePDB(input1,output2,l_index1, outputFileName)

if __name__	== "__main__":
	main()

