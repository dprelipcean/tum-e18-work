#!/usr/bin/python
# splinalysis.py
# Created: 2019-11-06 16:51:40.494705
# Author: Fabian Krinner
import os, sys
import numpy as np
import numpy.linalg as la
import ROOT

from regularize_integral_matrix import parseMatrixFile
from splineFunctions            import makeBinnedSplineFunctions

GLOBAL_NFIXED = 7

cpls = np.array(
       [ 1.01108  - 0.196533j,
         0.675822 + 0.750576j,
        -1.95532  - 0.240083j,
         1.       + 0.j,
         0.675822 + 0.750576j,
         0.397326 + 0.211262j,
         0.117557 - 0.161803j ]
       , dtype = complex)

binning_S = [ 0.401002 , 0.501002 , 0.601002 , 0.701002 , 0.801002 , 0.901002 , 
              1.001    , 1.101    , 1.201    , 1.301    , 1.401    , 1.501    ,
              1.601    , 1.701    , 1.801    , 1.901    , 2.001    , 2.101    ,
              2.201    , 2.301    , 2.401    , 2.501    , 2.601    , 2.701    ,
              2.801    , 2.99293 ]

binning_P = [ 0.401002 , 0.501002 , 0.601002 , 0.701002 , 0.726002 , 0.751002 ,
              0.776002 , 0.801002 , 0.826002 , 0.851002 , 0.876002 , 0.901002 ,
              0.926002 , 0.951002 , 0.976002 , 1.001    , 1.101    , 1.201    ,
              1.301    , 1.401    , 1.501    , 1.601    , 1.701    , 1.801    ,
              1.901    , 2.001    , 2.101    , 2.201    , 2.301    , 2.401    ,
              2.501    , 2.601    , 2.701    , 2.801    , 2.99293 ]

binning_D = [ 0.401002 , 0.501002 , 0.601002 , 0.701002 , 0.801002 , 0.901002 ,
              1.001    , 1.101    , 1.201    , 1.301    , 1.401    , 1.501    ,
              1.601    , 1.701    , 1.801    , 1.851    , 1.901    , 1.951    ,
              2.001    , 2.051    , 2.101    , 2.151    , 2.201    , 2.251    ,
              2.301    , 2.351    , 2.451    , 2.551    , 2.651    , 2.751    ,
              2.851    , 2.99293 ]

GLOBAL_NWAVES = [len(binning_S) - 1, len(binning_P) - 1, len(binning_D) - 1]

def cutFixedFromMatrix(matrix):
	return matrix[GLOBAL_NFIXED:,GLOBAL_NFIXED:]

def getNormedMatrix(matrix):
	nm = np.outer(np.diagonal(matrix)**.5, np.diagonal(matrix)**.5)
	nm[nm != 0.] = 1/nm[nm!=0.]
	return matrix * nm

def getNwaves(L, degree):
	return GLOBAL_NWAVES[L] + degree + 1

def getBordersWithDegree(degree):
	borders = [0]
	for L in range(3):
		borders.append(borders[~0] + getNwaves(L,degree))
	return borders

def getBordersWithDegreeAndFixed(degree):
	bwd = getBordersWithDegree(degree)
	return [b + GLOBAL_NFIXED for b in bwd]

def loadResultFile(inFileName):
	vals = []
	with open(inFileName, 'r') as inFile:
		for line in inFile:
			vals += [float(v) for v in line.split()]
	vals = np.array(vals)
	return vals[::2] + 1.j*vals[1::2]

def realizeMatrix(complexMatrix):
	dim = complexMatrix.shape[0]
	if complexMatrix.shape[1] != dim:
		raise ValueError("Only valid for square matrices")
	realMatrix = np.zeros((2*dim, 2*dim))
	realMatrix[ ::2, ::2] = complexMatrix.real
	realMatrix[ ::2,1::2] =-complexMatrix.imag
	realMatrix[1::2, ::2] = complexMatrix.imag
	realMatrix[1::2,1::2] = complexMatrix.real
	return realMatrix

def getIntegralFileName(degree):
	return "/home/iwsatlas1/fkrinner/Own_Software/ppppppppp/build/integralFiles/ps_integral_CP_MC_withFixed_degree3_K-pi+pi+_L-012_0_1_2.CP_MC_EQD_SL"
#	return "/home/iwsatlas1/fkrinner/Own_Software/ppppppppp/build/integralFiles/ps_integral_CP_MC_withFixed_degree" + str(degree) + "_K-pi+pi+_L-012_0_1_2.CP_MC_SLINETEGRAL"

def getSplineTransformFileName(degree):
	return "localSplineTransform_deg"+str(degree)+"__sEq.dat"
#	return "localSplineTransform_deg"+str(degree)+"__.dat" # THe "__" in the end saves it from being

def loadComplexIntegralMatrix(degree, withFixed):
	if not degree == 3:
		raise ValueError("Only degree in tegral exsists at the moment")

	inMatrix = parseMatrixFile(getIntegralFileName(degree))
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#	indices  = range(GLOBAL_NFIXED)                                         #
#	count    = GLOBAL_NFIXED                                                # 
#	for i in range(len(binning_S) + len(binning_P) + len(binning_D) -3):    # 
#		for d in range(6):                                              # 
#			if d < degree+2:                                        # 
#				indices.append(count)                           # 
#			count += 1                                              # 
#	inMatrix = inMatrix[np.ix_(indices,indices)]                            # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
	if withFixed:
		return inMatrix
	return cutFixedFromMatrix(inMatrix)

def loadSplineTransform(degree, withFixed):
	mat      = [ ]
	fileName = getSplineTransformFileName(degree)
	print "Loading spline transform from '" + fileName + "'"
	with open(fileName, 'r') as inFile:
		for line in inFile:
			line = [float(v) for v in line.split()]
#			if len(mat) == 0:
#				for i in range(GLOBAL_NFIXED):
#					unitLine = [0.]* len(line)
#					unitLine[i] = 1.
#					mat.append(unitLine)
			mat.append(line)
	retVal = np.array(mat).T
	if withFixed:
		return retVal
	return cuitFixedFromMatrix(retVal)

def makeSmoothSplineEval(coeffs, splineTransform, functions, masses, degree, filterL = None):
	fax = np.array([1.] * len(coeffs))
	if filterL is not None:
		bordersWithDegree = getBordersWithDegree(degree)		
		if   filterL == 0:
			for i in range(bordersWithDegree[1] , bordersWithDegree[3]):
				fax[i] = 0.
		elif filterL == 1:
			for i in range(bordersWithDegree[1]):
				fax[i] = 0.
			for i in range(bordersWithDegree[2], bordersWithDegree[3]):
				fax[i] = 0.
		elif filterL == 2:
			for i in range(bordersWithDegree[2]):
				fax[i] = 0.
		else:
			raise ValueError("Unknwon wave " + str(i_isob))
	retVal = np.array([[f(m) for m in masses] for f in functions])
	retVal = np.dot(splineTransform.T,retVal)
	retVal = np.dot(coeffs*fax, retVal)
	return retVal

def makeZeroModes(matrix, minEV):
	val, vec = la.eig(matrix)
	dim      = len(matrix)
	retVal   = [ ]
	for i,v in enumerate(val):
		if v.real < 0.:
			print "WARNING: Eigenvalue smaller than zero:",v
		elif v.real < minEV and v.real != 0.:
			print v.real
			retVal.append(vec[:,i].real)
	return retVal

def prepareForInversion(matrix, nApprox, nFixed = GLOBAL_NFIXED):
	A = matrix[2*nFixed:, 2*nFixed:]
	B = -matrix[2*nApprox,2*nFixed:] - matrix[2*nFixed:,2*nApprox]
	C = matrix[2*nApprox, 2*nApprox]
	return A,B,C

def getSingleLbestSolution(i_isob, realMatrix, degree, nFixed = GLOBAL_NFIXED):
	LLL     = [0, 0, 0, 1, 1, 1, 2]
	borders = getBordersWithDegreeAndFixed(degree)
	indices = range(GLOBAL_NFIXED) + range(borders[LLL[i_isob]], borders[LLL[i_isob] + 1])
	Indices = [ ]
	for i in indices:
		Indices.append(2*i  )
		Indices.append(2*i+1)
	Indices = np.array(Indices)
	mti = realMatrix[np.ix_(Indices, Indices)]
	iNDICES = [ ]
	for i in range(borders[LLL[i_isob]]-GLOBAL_NFIXED, borders[LLL[i_isob] + 1]-GLOBAL_NFIXED):
		iNDICES.append(2*i  )
		iNDICES.append(2*i+1)
	iNDICES = np.array(iNDICES)
	A,B,C = prepareForInversion(mti,i_isob, nFixed = nFixed)

	solution = -np.dot(la.pinv(A + A.T), B)

	retVal = np.zeros(2*borders[~0] - 2*GLOBAL_NFIXED)
	retVal[iNDICES] = solution
	return retVal

def resolveZeroMode(i_isob, solution, zeroModes, realMatrix, degree, nFixed = GLOBAL_NFIXED):
	singleLsolution = getSingleLbestSolution(i_isob, realMatrix, degree, nFixed = nFixed)
	delta           = solution - singleLsolution
	sps             = [np.dot(zm, delta) for zm in zeroModes]
	print np.sum(delta**2),
	for i in range(len(sps)):
		delta -= sps[i] * zeroModes[i]
	print np.sum(delta**2), "this is the improvement"
	return singleLsolution + delta

def getBestSolutionToFixedShape(i_isobs, complexMatrix, zeroModes, degree, nFixed = GLOBAL_NFIXED):
	dim           = len(complexMatrix)
	realZeroModes = [ ]
	for zm in zeroModes:
		rzm       = np.zeros(2*len(zm))
		izm       = np.zeros(2*len(zm))
		rzm[ ::2] = zm.real
		izm[1::2] = zm.real
		realZeroModes.append(rzm)
		realZeroModes.append(izm)

	realMatrix = realizeMatrix(complexMatrix)
	retVal     = [ ]
	for i_isob in i_isobs:
		A, B, C      = prepareForInversion(realMatrix, i_isob, nFixed = nFixed)
		bestSolution = -np.dot(la.pinv(A + A.T), B)
		bestSolution = resolveZeroMode(i_isob, bestSolution, realZeroModes, realMatrix, degree)
		retVal.append(bestSolution)
	return retVal

def quickplot(arrays, outFileName):
	hists = [ ]
	c1    = ROOT.TCanvas()
	for a, array in enumerate(arrays):
		hists.append( ROOT.TH1D("h"+str(a),"h"+str(a),len(array), 0.,1.))
		for i,v in enumerate(array):
			hists[~0].SetBinContent(i+1, v)
		if a == 0:
			hists[~0].Draw()
			maxx = hists[0].GetMaximum()
		else:
			hists[~0].SetLineColor(a+1)
			hists[~0].SetMaximum(maxx)
			hists[~0].Draw("SAME")
	c1.Print(outFileName)

def abs2(v):
	return v.real**2 + v.imag**2

def main():
	degree = 3

	smax   = 2.9929346001000003

	fullComplexBinIntegralMatrix = loadComplexIntegralMatrix(degree, True)
	fullSplineTransform          = loadSplineTransform(degree, True)

	fullComplexSplineIntegralMatrix  = np.dot(fullSplineTransform.T,np.dot(fullComplexBinIntegralMatrix, fullSplineTransform))
	fullSplineIntegralNorms          = np.array([1./fullComplexSplineIntegralMatrix[i,i].real**.5 for i in range(len(fullComplexSplineIntegralMatrix))])
	
	normedFullComplexSplineIntegralMatrix = getNormedMatrix(fullComplexSplineIntegralMatrix)

	splineZeroModesFromNormedMatrix  = makeZeroModes(cutFixedFromMatrix(normedFullComplexSplineIntegralMatrix), 1.e-6)
	print len(splineZeroModesFromNormedMatrix),"zero-modes found."

	bestApproxes = getBestSolutionToFixedShape(range(GLOBAL_NFIXED), normedFullComplexSplineIntegralMatrix, splineZeroModesFromNormedMatrix, degree)

#	for i,bapp in enumerate(bestApproxes):
#		quickplot([abs2(bapp[::2] + 1.j*bapp[1::2])], "bapp_" + str(i) + ".pdf")
#	return 

	splineFcns  = [ ]
	splineFcns += makeBinnedSplineFunctions(degree, binning_S)
	splineFcns += makeBinnedSplineFunctions(degree, binning_P)
	splineFcns += makeBinnedSplineFunctions(degree, binning_D)

	startVals = np.zeros(len(bestApproxes[0])/2, dtype = complex)
	for i in range(GLOBAL_NFIXED):
#		if not i in [0,3]:
#			continue
		startVals += (bestApproxes[i][::2] + 1.j*bestApproxes[i][1::2]) * cpls[i]  #/ fullSplineIntegralNorms[i] no normalizations needed, since the bestApproxes are normed already

	startValueResult = "/home/iwsatlas1/fkrinner/Own_Software/ppppppppp/build/MC_fits/real_params/012_1573233168_-13099840.240761_deg3_dataLstringSPD_0_1_2.CP_MC"
	result = loadResultFile(startValueResult)

	print len(result), len(startVals)

	sss   = np.linspace(binning_S[0],smax, 1000)			
	evall = [makeSmoothSplineEval(result * fullSplineIntegralNorms[GLOBAL_NFIXED:] , cutFixedFromMatrix(fullSplineTransform), splineFcns, sss, degree, filterL = i) for i in range(3)]
	with open("sv3magebei" , 'w') as outFile:
#		for i,s in enumerate(sss):
#			outFile.write(str(i) + " " + str(s))
#			for evl in evall:
#				outFile.write(" " + str(evl[i].real) + " " + str(evl[i].imag))
#			outFile.write('\n')

		for i,v in enumerate(startVals):
			norm = fullSplineIntegralNorms[GLOBAL_NFIXED+i]
			outFile.write(str(v.real * norm) + " " +str(v.imag * norm) + '\n')
	return 
#	evall = [makeSmoothSplineEval((coeffs[::2] + 1.j*coeffs[1::2])*fullSplineIntegralNorms[GLOBAL_NFIXED:], cutFixedFromMatrix(fullSplineTransform), splineFcns, sss, degree, filterL = 0) for coeffs in bestApproxes]
#	with open("repr", 'w') as outFile:
#		for i,s in enumerate(sss):
#			outFile.write(str(i) + " " + str(s))
#			for evl in evall:
#				outFile.write(" " + str(evl[i].real) + " " + str(evl[i].imag))
#			outFile.write("\n")
#	evall = [makeSmoothSplineEval((coeffs[::2] + 1.j*coeffs[1::2])*fullSplineIntegralNorms[GLOBAL_NFIXED:], cutFixedFromMatrix(fullSplineTransform), splineFcns, sss, degree, filterL = 1) for coeffs in bestApproxes]
#	with open("Prepr", 'w') as outFile:
#		for i,s in enumerate(sss):
#			outFile.write(str(i) + " " + str(s))
#			for evl in evall:
#				outFile.write(" " + str(evl[i].real) + " " + str(evl[i].imag))
#			outFile.write("\n")
#	evall = [makeSmoothSplineEval((coeffs[::2] + 1.j*coeffs[1::2])*fullSplineIntegralNorms[GLOBAL_NFIXED:], cutFixedFromMatrix(fullSplineTransform), splineFcns, sss, degree, filterL = 2) for coeffs in bestApproxes]
#	with open("Drepr", 'w') as outFile:
#		for i,s in enumerate(sss):
#			outFile.write(str(i) + " " + str(s))
#			for evl in evall:
#				outFile.write(" " + str(evl[i].real) + " " + str(evl[i].imag))
#			outFile.write("\n")


if __name__ == "__main__":
	main()
