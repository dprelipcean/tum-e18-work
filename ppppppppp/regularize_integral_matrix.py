#!/usr/bin/python
# regularize_integral_matrix.py
# Created: 2018-07-20 15:35:27.287781
# Author: Fabian Krinner
import os, sys
import numpy as np
import numpy.linalg as la
from getBranchFileEnding import getBranchFileEnding, getIntegralFileEnding

def isHermitian(matrix, numLim = 1.e-15):
	dim = len(matrix)
	nhCount = 0
	maxDiff = 0.
	for i in range(dim):
		for j in range(i+1):
			diff = abs(matrix[i,j] - np.conj(matrix[j,i]))
			if diff > numLim:
				nhCount += 1
				maxDixx = max(maxDiff, diff)
	if nhCount == 0:
		return True
	print nhCount,"non-hermitian pairs of entries (maxDiff =", maxDiff,")"
	return False

def toComplex(s):
	chunks = s[1:~0].split(',')
	return float(chunks[0]) + 1.j*float(chunks[1])

def toString(val):
	return '('+str(val.real)+','+str(val.imag)+')'

def parseMatrixFile(inFileName):
	matrix = []
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			lineVals = []
			chunks = line.split()
			for c in chunks:
				lineVals.append(toComplex(c))
			matrix.append(lineVals)
	return np.asarray(matrix, dtype = complex)

def writeMatrixFile(outFileName, matrix):
	with open(outFileName, 'w') as outFile:
		for i in range(len(matrix)):
			for j in range(len(matrix)):
				outFile.write(toString(matrix[i,j])+' ')
			outFile.write('\n')

def regulatrizeMatrix(matrix, silent = False):
	iterCount = 0
	while True:
		val, vec = la.eig(matrix)
		indices = []
		for i,v in enumerate(val):
			if v.real < 0.:
				if not silent:
					print v," !!!"
				indices.append(i)
		if len(indices) == 0:
			break
		for i in indices:
			for j in range(len(val)):
				for k in range(len(val)):
					matrix[j,k] -= 2*val[i]*vec[j,i]*vec[k,i]
		iterCount += 1
		if not silent:
			print "interation",iterCount
		for i in range(len(matrix)):
			matrix[i,i] = matrix[i,i].real
	return matrix

def regularize(inFileNames, outFileName):
	matrix = parseMatrixFile(inFileNames[0])
	for i in range(1,len(inFileNames)):
		matrix += parseMatrixFile(inFileNames[i])
	writeMatrixFile(outFileName, regulatrizeMatrix(matrix))

def main():
	bfe = '.' + getIntegralFileEnding()
	if len(sys.argv) == 1:
		folder = "./build/integralFiles/"
		for fn in os.listdir(folder):
			if 'regular' in fn:
				continue
			print fn,"<-------------------"
			if fn.endswith(bfe):
				inFileName = folder + fn
				print inFileName
#				mm = parseMatrixFile(inFileName)
				outFileName = inFileName.replace(bfe,"_regular"+bfe)
				regularize([inFileName], outFileName)
	else:
		inFileName = sys.argv[1]
		outFileName = inFileName.replace(bfe,"_regular"+bfe)
		regularize([inFileName], outFileName)

if __name__ == "__main__":
	main()
