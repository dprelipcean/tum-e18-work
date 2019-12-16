#!/usr/bin/python
# integralAnalysis.py
# Created: 2019-08-30 10:33:37.091734
# Author: Fabian Krinner
import os, sys
import numpy as np
import numpy.linalg as la
from scipy.linalg import block_diag
from regularize_integral_matrix import parseMatrixFile#, isHermitian, regulatrizeMatrix
from splineFunctions import makeBinnedSplineFunctions
import ROOT

GLOBAL_NFIXED = 7

cpls = [ 1.01108  - 0.196533j,
         0.675822 + 0.750576j,
        -1.95532  - 0.240083j,
         1.       + 0.j,
         0.675822 + 0.750576j,
         0.397326 + 0.211262j,
         0.117557 - 0.161803j ]

binning_S = [ 0.401002, 0.453262, 0.508721, 0.567381, 0.629241, 0.694301, 0.76256,
              0.83402,  0.90868,  0.98654,  1.0676,   1.15186,  1.23932,  1.32998,
              1.42384,  1.5209,   1.62116,  1.72462,  1.83128,  1.94114,  2.0542,
              2.17046,  2.28992,  2.41258,  2.53844,  2.6675,   2.79976,  2.93522,
              3.07388 ]

binning_P = [ 0.401002, 0.453262, 0.508721, 0.567381, 0.629241, 0.645206, 0.661371,
              0.677736, 0.694301, 0.711066, 0.72803,  0.745195, 0.76256,  0.780125,
              0.79789,  0.815855, 0.83402,  0.852385, 0.87095,  0.889715, 0.90868,
              0.927845, 0.94721,  0.966775, 0.98654,  1.0065,   1.08836,  1.17342,
              1.26168,  1.35314,  1.4478,   1.54566,  1.64672,  1.75098,  1.85844,
              1.9691,   2.08296,  2.20002,  2.32028,  2.44374,  2.5704,   2.70026,
              2.83332,  2.96958,  3.10904 ]

binning_D = [ 0.401002, 0.453262, 0.508721, 0.567381, 0.629241, 0.694301, 0.76256,
              0.83402,  0.90868,  0.98654,  1.0676,   1.15186,  1.23932,  1.32998,
              1.42384,  1.5209,   1.62116,  1.72462,  1.83128,  1.88581,  1.94114,
              1.99727,  2.0542,   2.11193,  2.17046,  2.22979,  2.28992,  2.35085,
              2.47511,  2.60257,  2.73323,  2.86709,  3.00415 ]

def loadSplineTransform(deg):
	fileName = "localSplineTransform_deg"+str(deg)+".dat"
	mat      = [ ]
	with open(fileName, 'r') as inFile:
		for line in inFile:
			line = [0.]*GLOBAL_NFIXED + [float(v) for v in line.split()]
			if len(mat) == 0:
				for i in range(GLOBAL_NFIXED):
					unitLine = [0.]* len(line)
					unitLine[i] = 1.
					mat.append(unitLine)
			mat.append(line)
	retVal = np.array(mat).T
	return retVal

def makeSplineTransform(nWaves, to_DOS):
	N = 0
	for n in nWaves:
		N += n
	nSpline = 0
	for toDO in to_DOS:
		if toDO =='lin':
			nSpline += 1
	retVal   = np.zeros((2*N, N+nSpline))
	countIn  = 0
	countOut = 0
	for i,n in enumerate(nWaves):
		if to_DOS[i] == 'lin':
			for _ in range(n):
				retVal[countOut, countIn  ]   =  .5
				retVal[countOut, countIn+1]   =  .5
				retVal[countOut+N, countIn]   =-1.
				retVal[countOut+N, countIn+1] = 1.
				countOut += 1
				countIn  += 1
			countIn += 1
		elif to_DOS[i] == 'const':
			for _ in range(n):
				retVal[countOut,countIn] = 1.
				countOut += 1
				countIn  += 1
		elif to_DOS[i] == 'zero':
			countOut += n
			countIn  += n
		else:
			raise ValueError("Unknown toDo: "+to_DOS[i]+".")
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

def norm(matrix):
	nm = np.outer(np.diagonal(matrix)**.5, np.diagonal(matrix)**.5)
	nm[nm != 0.] = 1/nm[nm!=0.]
	matrix *= nm

def prepareForInversion(matrix, nApprox, nFixed = GLOBAL_NFIXED):
	A = matrix[2*nFixed:, 2*nFixed:]
	B = -matrix[2*nApprox,2*nFixed:] - matrix[2*nFixed:,2*nApprox]
	C = matrix[2*nApprox, 2*nApprox]
	return A,B,C

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

def zeroModeResolve(vector, zeroModes, indices):
	dim = vector.shape[0]/2
	fullZM = np.zeros((2*dim, 2*len(zeroModes)))
	for z,zm in enumerate(zeroModes):
		relZM = np.zeros(dim)
		relZM[indices] = zm[indices]
		fullZM[ ::2,2*z  ] = relZM
		fullZM[1::2,2*z+1] = relZM
	relevantData = np.zeros(2*dim)
	relevantData[2*indices  ] = vector[2*indices  ]
	relevantData[2*indices+1] = vector[2*indices+1]

	A = np.dot(fullZM.T, fullZM)
	B = -2*np.dot(fullZM.T, relevantData)
	C = np.dot(relevantData, relevantData)

	cpls = -np.dot(la.pinv(A),B)

	retVal = np.zeros(2*dim)
	for z,zm in enumerate(zeroModes):
		retVal[ ::2] -= cpls[2*z  ]*zm
		retVal[1::2] -= cpls[2*z+1]*zm
	retVal += vector
	return retVal

sectorBorders          = [0,28,72,104]
sectorBordersWithFixed = [b + GLOBAL_NFIXED for b in sectorBorders]

def getWaveIndices(L, toDOS):
	borders = [0,28,72,104]
	for i in range(3):
		if toDOS[i] == 'lin':
			for j in range(i+1,4):
				borders[j] += 1
	return range(borders[L],borders[L+1])

def getAntiWaveIndices(L, toDOS):
	retVal = [ ]
	for l in range(3):
		if l == L:
			continue
		retVal += getWaveIndices(l,toDOS)
	return np.array(retVal, dtype = int)

def getSingleLbestSolution(i_isob, realMatrix, borders):
	LLL     = [0, 0, 0, 1, 1, 1, 2]
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
	A,B,C = prepareForInversion(mti,i_isob)

	solution = -np.dot(la.pinv(A + A.T), B)

	retVal = np.zeros(2*borders[~0] - 2*GLOBAL_NFIXED)
	retVal[iNDICES] = solution
	return retVal

def resolveZeroMode(i_isob, solution, zeroModes, realMatrix, borders):
	singleLsolution = getSingleLbestSolution(i_isob, realMatrix, borders)
	delta           = solution - singleLsolution
	sps             = [np.dot(zm, delta) for zm in zeroModes]
	print np.sum(delta**2),
	for i in range(len(sps)):
		delta -= sps[i] * zeroModes[i]
	print np.sum(delta**2), "this is the improvement"
	return singleLsolution + delta

def getFunctionCoefficients(result, basisFunctions):
	dim              = len(result)
	realResult       = np.zeros(2*dim)
	realResult[ ::2] = result.real
	realResult[1::2] = result.imag

	nf = len(basisFunctions)
	funcAmpls = np.zeros((2*nf,2*dim))
	for i in range(nf):
		if len(basisFunctions[i]) != dim:
			raise ValueError("Wrong size of basis functions: " + str(len(basisFunctions[i])) + " instead of " + str(dim) + ".")
		funcAmpls[2*i  , ::2] = basisFunctions[i].real
		funcAmpls[2*i  ,1::2] =-basisFunctions[i].imag
		funcAmpls[2*i+1, ::2] = basisFunctions[i].imag
		funcAmpls[2*i+1,1::2] = basisFunctions[i].real

	A =    np.dot(funcAmpls,  funcAmpls.T)
	B = -2*np.dot(funcAmpls,  realResult)
	C =    np.dot(realResult, realResult)
#	for i in range(2*nf):
#		for j in range(2*nf):
#			if A[i,j] == 0.:
#				print "0",
#			else:
#				print "1",
#		print
#	raise Exception

	retVal = -np.dot(la.inv(A+A.T),B)
	return retVal

def grame(fcns):
	nd = len(fcns[0])
	nf = len(fcns)
	allvals = np.zeros((nf, 2*nd))
	for i,f in enumerate(fcns):
		allvals[i, ::2] = f.real
		allvals[i,1::2] = f.imag
	ube = np.dot(allvals,allvals.T)
	norm(ube)

	val,vec = la.eig(ube)

	print val

def abs2(z):
	return z.real**2 + z.imag**2

def EVAL(coeffs, splineTransform, functions, masses, L = None, bordersWithDegree = None):
	fax = np.array([1.] * len(coeffs))
	if L is not None:
		if bordersWithDegree is None:
			raise RuntimeError("L is not None, but no borders given.")
		if   L == 0:
			for i in range(bordersWithDegree[1] , bordersWithDegree[3]):
				fax[i] = 0.
		elif L == 1:
			for i in range(bordersWithDegree[1]):
				fax[i] = 0.
			for i in range(bordersWithDegree[2], bordersWithDegree[3]):
				fax[i] = 0.
		elif L == 2:
			for i in range(bordersWithDegree[2]):
				fax[i] = 0.
		else:
			raise ValueError("Unknwon wave " + str(i_isob))

	retVal = np.array([[f(m) for m in masses] for f in functions])
	retVal = np.dot(splineTransform.T,retVal)
	retVal = np.dot(coeffs*fax, retVal)
	return retVal

def main():
	i_isob     = int(sys.argv[1])
	from time   import sleep
	from random import random
	while True:
		print "Alles Mist!",
		if random() > .9:
			print
			sleep(0.1)
	

	L          = {0:0,1:0,2:0,3:1,4:1,5:1,6:2}[i_isob]
#	inFileName = "/home/iwsatlas1/fkrinner/Own_Software/ppppppppp/studyIntegralLin.ps"
	inFileName = "/home/iwsatlas1/fkrinner/Own_Software/ppppppppp/build/integralFiles/ps_integral_CP_MC_withFixed_degree4_K-pi+pi+_L-012_0_1_2.CP_MC_SLINETEGRAL"
#	inFileName = "/home/iwsatlas1/fkrinner/Own_Software/ppppppppp/build/integralFiles/ps_integral_CP_MCK-pi+pi+_L-012_0_1_2.CP_MC_SLINETEGRAL"

	complexMatrix = parseMatrixFile(inFileName)

	JPC_to_dos = sys.argv[1:4]
	cntLin     = 0
	degree     = 4
	for td in JPC_to_dos:
		if td == 'lin':
			cntLin += 1

	bordersWithDegreeAndFixed = [sectorBordersWithFixed[i] + i*(degree+1) for i in range(len(sectorBordersWithFixed))]
	bordersWithDegree         = [sectorBorders[i]          + i*(degree+1) for i in range(len(sectorBorders))         ]

	splineFcns  = [ ]
	splineFcns += makeBinnedSplineFunctions(degree, binning_S)
	splineFcns += makeBinnedSplineFunctions(degree, binning_P)
	splineFcns += makeBinnedSplineFunctions(degree, binning_D)

	nWaves          = [ 28, 44, 32 ]
	splineTransform = loadSplineTransform(degree)
	complexMatrix   = np.dot(splineTransform.T,np.dot(complexMatrix, splineTransform))
	unnormedCM      = np.copy(complexMatrix)
	normsFixed      = np.array([1./complexMatrix[i,i].real**.5 for i in range(GLOBAL_NFIXED)])
	norms           = np.array([1./complexMatrix[i,i].real**.5 for i in range(GLOBAL_NFIXED,len(complexMatrix))])

	norm(complexMatrix)

	zms     = makeZeroModes(complexMatrix[GLOBAL_NFIXED:,GLOBAL_NFIXED:], 1.e-6)
	realZMS = [ ]

	for zm in zms:
		rzm       = np.zeros(2*len(zm))
		izm       = np.zeros(2*len(zm))
		rzm[ ::2] = zm.real
		izm[1::2] = zm.real
		realZMS.append(rzm)
		realZMS.append(izm)
	print len(zms),'zero modes'

	dim          = len(complexMatrix)
	realMatrix   = realizeMatrix(complexMatrix)
	A, B, C      = prepareForInversion(realMatrix, i_isob)
	bestSolution = -np.dot(la.pinv(A + A.T), B)
	bestSolution = resolveZeroMode(i_isob, bestSolution, realZMS, realMatrix, bordersWithDegreeAndFixed)
	startValues  = None

	for i in range(7):
		sol      = getSingleLbestSolution(i, realMatrix, bordersWithDegreeAndFixed)
		complSol = cpls[i] / normsFixed[i] * (sol[ ::2] + 1.j*sol[1::2])
		if i == 0:
			startValues  = complSol
		else:
			startValues -= complSol

	startValues *= norms

	with open("StartValuesSplineFit_deg"+str(degree), 'w') as outFile:
		for f in startValues:
			outFile.write(str(f.real) + ' ' + str(f.imag) + '\n')

	bestResultFileName = "build/MC_fits/real_params/012_1571556784_-13099820.428558_deg4_dataLstringSPD_0_1_2.CP_MC"
#	bestResultFileName = "build/MC_fits/real_params/012_1571907812_-13099814.461731_deg4_dataLstringSPD_0_1_2.CP_MC"
	vals = [ ]
	with open(bestResultFileName, 'r') as inFile:
		for line in inFile:
			vals += [ float(v) for v in line.split() ]
	vals = np.array(vals)
	vals = (vals[ ::2] + 1.j*vals[1::2])  * norms

#	vals = np.conjugate( vals )

	allFcns = [ startValues ]
	for zm in zms:
		zm *= norms
		allFcns.append(zm)

#	sss   = np.linspace(.5, 2.8, 1000)
#	alll  = [ EVAL(vals, splineTransform[GLOBAL_NFIXED:,GLOBAL_NFIXED:], splineFcns, sss, L = int(sys.argv[1]), bordersWithDegree = bordersWithDegree) ]
#	alll.append(EVAL(norms, splineTransform[GLOBAL_NFIXED:,GLOBAL_NFIXED:], splineFcns, sss, L = int(sys.argv[1]), bordersWithDegree = bordersWithDegree))
#	for fcn in allFcns:
#		alll.append( EVAL(fcn, splineTransform[GLOBAL_NFIXED:,GLOBAL_NFIXED:], splineFcns, sss, L = int(sys.argv[1]), bordersWithDegree = bordersWithDegree) )

#	with open("theComparison_" + sys.argv[1], 'w') as outFile:
#		for i,s in enumerate(sss):
#			outFile.write(str(i) + " " + str(s))
#			for vls in alll:
#				outFile.write(" " + str(vls[i].real) + " " + str(vls[i].imag))
#			outFile.write("\n")
#	return

#	cffsss = getFunctionCoefficients(alll[0], alll[1:])
#	brv    = None
#	vals   = alll[0]
#	for i in range(1,len(alll)):
#		f = alll[i]
#		if i == 1:
##			brv = np.zeros(len(f), dtype = complex)
#			brv = (cffsss[2*(i-1)  ]+1.j*cffsss[2*(i-1)+1]) * f
#		else:
#			pass
#			vals -=  (cffsss[2*(i-1)  ]+1.j*cffsss[2*(i-1)+1]) * f


	vals        /= norms**2
#	startValues /= norms

	fakk = (np.sum(abs2(vals))/np.sum(abs2(startValues)))**.5
	startValues *= fakk

	I   = abs2(vals.real)
	II  = abs2(-startValues.real)
	III = vals.imag
	IV  = -startValues.imag

	quickplot([I],'awaiasistdas'+sys.argv[1]+'.pdf')

	return

	makeThePlotFrom = vals

	if   i_isob == 0:
		for i in range(bordersWithDegree[1], bordersWithDegree[3]):
			makeThePlotFrom[i] = 0.
	elif i_isob == 1:
		for i in range(bordersWithDegree[1]):
			makeThePlotFrom[i] = 0.
		for i in range(bordersWithDegree[2], bordersWithDegree[3]):
			makeThePlotFrom[i] = 0.
	elif i_isob == 2:
		for i in range(bordersWithDegree[2]):
			makeThePlotFrom[i] = 0.
	else:
		raise ValueError("Unknwon wave " + str(i_isob))

	sss   = np.linspace(binning_D[0],binning_D[~0], 1000)
	evals = EVAL(makeThePlotFrom, splineTransform[GLOBAL_NFIXED:,GLOBAL_NFIXED:], splineFcns, sss)
	with open("thePlot"+sys.argv[1], 'w') as outFile:
		for i,v in enumerate(evals):
			outFile.write(str(i) + ' ' + str(sss[i]) + ' ' + str(v.real) + ' ' + str(v.imag) + '\n')
	return 0

if __name__ == "__main__":
	main()
