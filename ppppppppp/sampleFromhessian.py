#!/usr/bin/python
# sampleFromhessian.py
# Created: 2018-12-06 13:33:51.387882
# Author: Fabian Krinner
import os, sys
import numpy as np
import numpy.linalg as la


def loadFile(inFileName):
	retVal = []
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			retVal.append([float(v) for v in line.split()])
	return np.array(retVal)

def checkHessian(matrix, numLim = 1.e-5):
	if not matrix.shape[0] == matrix.shape[1]:
		print("Matrix is not quadratic.")
		return False
	dim = len(matrix)
	for i in range(dim):
		for j in range(dim):
			if abs(matrix[i,j] - matrix[j,i]) > numLim:
				print("Matrix is not symmetric.")
				return False
	return True

def main():
	inFileName  = sys.argv[1]
	outFileName = inFileName.replace('~~','~~~~')

	nSample = 10000
	if len(sys.argv) > 2:
		nSample = int(sys.argv[2])

	hessian    = loadFile(inFileName)
	if not checkHessian(hessian):
		return False
	coma = la.inv(hessian)
	dim  = len(coma)

	mean = np.zeros(len(hessian))

	with open(outFileName,'w') as outFile:
		for i in range(nSample):
			smpl = np.random.multivariate_normal(mean, coma)
			for v in smpl:
				outFile.write(str(v) + ' ')
			outFile.write('\n')
	os.remove(inFileName)
	return True

if __name__ == "__main__":
	if main():
		sys.exit(0)
	else:
		sys.exit(1)
