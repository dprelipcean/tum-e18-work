#!/usr/bin/python
# check_strange_solution.py
# Created: 2019-03-07 15:15:00.689384
# Author: Fabian Krinner
import os, sys
import numpy as np
import numpy.linalg as la
from integral import toComplex
import ROOT

def is_pos_def(x):
    return np.all(la.eigvals(x) > 0)

def readFile(inFileName):
	retVal = []
	with open(inFileName, 'r') as inFile:
		for line in inFile:
			vals = [toComplex(v) for v in line.split()]
			retVal.append(vals)
	return np.array(retVal, dtype = complex)

def solveForLast(matrix):
	dim = len(matrix)-1
	if not matrix[dim,dim] == 1.:
		raise ValueError("Matrix not properly normalized.")
	A   = matrix[np.ix_(range(dim), range(dim))]
	B   = matrix[np.ix_(range(dim), [dim])]

	X   = -np.dot( la.inv(A), B) ## Check this

	return X

	delta = 1.e-7
	evl_center = np.dot(np.conjugate(X).T, np.dot(A, X)) + np.dot(np.conjugate(B).T,X) + np.dot(np.conjugate(X).T,B)
	for i in range(len(X)):
		X[i] += delta
		evl = np.dot(np.conjugate(X).T, np.dot(A, X)) + np.dot(np.conjugate(B).T,X) + np.dot(np.conjugate(X).T,B)
		print i,(evl - evl_center)/delta
		X[i] += delta*(1.j-1.)
		evl = np.dot(np.conjugate(X).T, np.dot(A, X)) + np.dot(np.conjugate(B).T,X) + np.dot(np.conjugate(X).T,B)
		print i,(evl - evl_center)/delta
		X[i] -= 1.j*delta


def getSubMatrix(fullMatrix, use):
	NNN = 10
	if not len(use) == NNN:
		raise ValueError("Need exactly " + str(NNN) + " flags.")
	numbers    = [28, 44, 32, 1,1,1,1,1,1,1] # freeS freeP freeD fixedS fixedP fixedD
	count      = 0
	indices    = []
	for w in range(NNN):
		if use[w]:
			for i in range(numbers[w]):
				indices.append(count + i)
		count += numbers[w]
	return fullMatrix[np.ix_(indices, indices)]

def getMinEvEv(matrix):
	val,vec = la.eig(matrix)
	minI    = None
	minEv   = float("inf")
	for i,v in enumerate(val):
		if v.imag > 1.e-10:
			raise ValueError("Complex valued eigenvalue encountered: " + str(v) + ".")
		if v.real < minEv:
			minI  = i
			minEv = v.real
	return val[minI].real, vec[:,minI]

def norm(matrix):
	dim = matrix.shape[0]
	for i in range(dim):
		for j in range(i):
			if i == j:
				continue
			fak = (matrix[i,i]*matrix[j,j])**.5
			if fak == 0.:
				matrix[i,j] = 0.
				matrix[j,i] = 0.
			else:
				matrix[i,j] /= fak
				matrix[j,i] /= fak
	for i in range(dim):
		if not  matrix[i,i] == 0.:
			matrix[i,i] = 1.

def plotVec(vec, stop = True):
	hist_re = ROOT.TH1D("h_re","h_re", len(vec), 0.,1.)
	hist_im = ROOT.TH1D("h_im","h_im", len(vec), 0.,1.)
	for i,v in enumerate(vec):
		hist_re.SetBinContent(i+1, v.real)
		hist_im.SetBinContent(i+1, v.imag)
	if stop:
		hist_re.SetMaximum(max(hist_re.GetMaximum(), hist_im.GetMaximum()))
		hist_re.SetMinimum(min(hist_re.GetMinimum(), hist_im.GetMinimum()))
		hist_im.SetLineColor(2)
		hist_re.Draw()
		hist_im.Draw("SAME")	

		raw_input("<enter_to_proceed>")

	return hist_re, hist_im

def mergeIndices(matrix, coefficients, indices):
	dim = len(matrix)
	mer = len(coefficients)
	if not mer == len(indices):
		raise ValueError("Number of indices to merge does not match")
	dimNew = dim - mer + 1
	newMatrix = np.zeros((dimNew, dimNew))
	count_i = 0
	for i in range(dim):
		if not i in indices:
			count_j = 0
			for j in range(dim):
				if not j in indices:
					newMatrix[count_i,count_j] = matrix[i,j]
					count_j += 1
			count_i +=1
	count_j = 0
	for j in range(dim):
		if not j in indices:
			for i,I in enumerate(indices):
				newMatrix[count_j,dim-mer] += matrix[j,I]*coefficients[i]
				newMatrix[dim-mer,count_j] += matrix[I,j]*np.conjugate(coefficients[i])
			count_j += 1	
	for i,I in enumerate(indices):
		newMatrix[dim-mer,dim-mer] += abs(coefficients[i])**2 * matrix[I,I]
	return newMatrix


def main():
	CLEO = [1.01108-0.196533j,0.675822+0.750576j,-1.95532-0.240083j,1.+0.j,0.675822+0.750576j,0.397326+0.211262j,0.117557-0.161803j]


	inFileName = "/home/iwsatlas1/fkrinner/Own_Software/ppppppppp/build/integralFiles/ps_integral_CP_MC_freed_and_fixedK-pi+pi+_L-012_0_1_2.CP_MC"
	matrix     = readFile(inFileName)

	norm(matrix)

	X = solveForLast(matrix)
	plotVec(X)
	return

	matrix_P_892 = getSubMatrix(matrix, [False, True,True, False, False, False, True, False, False, True])
	dim = len(matrix)
	newMatrix    = mergeIndices(matrix, CLEO, [dim-1-v for v in range(7)])
	norm(newMatrix)
	val, vec     =  getMinEvEv(newMatrix)
	print val
	plotVec(vec)

if __name__ == "__main__":
	main()

