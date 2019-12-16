#!/usr/bin/python
# analyzeBELLE_version2.py
# Created: 2018-10-08 10:02:46.680053
# Author: Fabian Krinner
import os, sys
import numpy as np
import numpy.linalg as la
freeNwaves  = [ 28, 43, 32, 28, 43, 32, 40, 34, 33]
fixedNwaves = [  1,  2,  1,  1,  3,  1,  1,  3,  1]

import freeDdecayBinnings as fDdb

from regularize_integral_matrix import parseMatrixFile, isHermitian, regulatrizeMatrix
import ROOT
from rootfabi import root_open
import datetime
from getBranchFileEnding import getBranchFileEnding, getIntegralFileEnding

from loadConstants import loadConstants

INF = float('inf')

def getNbin(mass):
	return int((mass-.5)/.04)

class integralClass:
	def __init__(self, inFileName):
		if inFileName.split(os.sep)[~0].startswith("ac_"):
			self.acFileName = inFileName
			self.psFileName = inFileName.replace("ac_","ps_")
		elif inFileName.split(os.sep)[~0].startswith("ps_"):
			self.acFileName = inFileName.replace("ps_","ac_")
			self.psFileName = inFileName
		else:
			raise ValueError("'" + inFileName + "' is not a valid matrix file name (does not start with 'ac_' or 'ps_'")
		self.acMatrix = parseMatrixFile(self.acFileName)
		self.psMatrix = parseMatrixFile(self.psFileName)
		if not isHermitian(self.acMatrix):
			raise ValueError("ps matrix is not hermitian")
		if not isHermitian(self.psMatrix):
			raise ValueError("ac matrix is not hermitian")

		self.dim      = len(self.acMatrix)
		self.freeMap  = getFreeMap(self.psFileName)
		self.fixMap   = None
		norm          = []
		for i in xrange(self.dim):
			nn = self.psMatrix[i,i]
			if nn == 0.:
				nn = 0.
			else:
				nn = 1./nn.real**.5
			norm.append(nn)
		self.norm   = np.asarray(norm)
		self.normed = False
		self.eigen  = False
		self.cut    = False
		self.numLim = 1.e-14

	def cutFreed(self):
		self.norm = np.asarray(cutFree(self.norm, self.freeMap, fixMap = self.fixMap))
		acMatrix  = []
		psMatrix  = []
		for i in xrange(self.dim):
			psMatrix.append(cutFree(self.psMatrix[i], self.freeMap, fixMap = self.fixMap))
			acMatrix.append(cutFree(self.acMatrix[i], self.freeMap, fixMap = self.fixMap))
		self.acMatrix = np.asarray(cutFree(acMatrix, self.freeMap, fixMap = self.fixMap), dtype = complex)
		self.psMatrix = np.asarray(cutFree(psMatrix, self.freeMap, fixMap = self.fixMap), dtype = complex)
		self.dim = len(self.norm)
		self.cut = True

	def cutFirstSector(self):
		nWaves = 0
		for i in xrange(9):
			if self.freeMap[i]:
				nWaves += freeNwaves[i]
			else:
				nWaves += fixedNwaves[i]
		self.norm = self.norm[:nWaves]
		self.acMatrix = self.acMatrix[:nWaves,:nWaves]
		self.psMatrix = self.psMatrix[:nWaves,:nWaves]
		self.dim = nWaves

	def makeZMs(self, maxZMev):
		val, vec = la.eig(self.psMatrix)
		self.allEVs = np.real(val)
		zms = []
		evs = []
		for i,v in enumerate(val):
			if abs(v.imag) > self.numLim:
				raise ValueError("Non vanishing imaginary part of eigenvalue: " + str(v))
			if v.real < 0.:
				raise ValueError("Negative eigenvalue encountered: "+ str(v))
			if v.real < maxZMev:
				zms.append(np.real(vec[:,i]))
				evs.append(v)
		self.zms       = np.asarray(zms)
		self.evs       = np.asarray(evs)
		self.zmsNormed = self.normed
		self.zmsCut    = self.cut
		self.eigen     = True

	def normalize(self):
		if self.normed:
			print "integralClass.norm(): WARNING: Already normed... no nothing"
			return
		normmatrix = np.diag(self.norm)
		self.acMatrix = regulatrizeMatrix(np.dot(normmatrix, np.dot(self.acMatrix, normmatrix)), silent = True)
		self.psMatrix = regulatrizeMatrix(np.dot(normmatrix, np.dot(self.psMatrix, normmatrix)), silent = True)

def makeFreeMap(string):
	if not len(string) == 9:
		raise ValueError("Could not makeFreeMap('" + string + "')")
	retVal = []
	for c in string:
		if c == '0':
			retVal.append(False)
		elif c == '1':
			retVal.append(True)
		else:
			raise ValueError("Could not makeFreeMap('"+string+"')")
	return retVal

def getFreeString(inFileName):
	lastPart = inFileName.split(os.sep)[~0]
	if lastPart.startswith("BELLE_fit_"):
		string = lastPart[10:19]
	elif "_integral_model_" in lastPart:
		string = lastPart[18:27]
	elif "BELLE_startValues" in lastPart:
		string = lastPart[17:26]
	return string

def getFreeMap(inFileName):
	string = getFreeString(inFileName)
	try:
		retVal = makeFreeMap(string)
	except ValueError:
		raise ValueError("Could not getFreeMap('" + inFileName + "')")
	return retVal

def getIntegralFileNames(inFileName):
	bfe = '.' + getIntegralFileEnding()
	string = getFreeString(inFileName)
	ps_fileName = "/home/iwsatlas1/fkrinner/Own_Software/ppppppppp/build/integralFiles/ps_integral_model_" + string + "_regular" + bfe
	ac_fileName = "/home/iwsatlas1/fkrinner/Own_Software/ppppppppp/build/integralFiles/ac_integral_model_" + string + "_regular" + bfe
	return ps_fileName, ac_fileName

def getHessianFileName(inFileName):
	if not "BELLE_startValues" in inFileName:
		theFileName = inFileName
	else:
		theFileName = getResultFileName(inFileName)
	return theFileName.replace("BELLE_fit_results/BELLE_fit_","BELLE_fit_results/hessians/BELLE_hessian_")

def getCOMAfileName(inFileName):
	if os.path.isfile(inFileName + ".COMA"):
		return inFileName + ".COMA"
	if not "BELLE_startValues" in inFileName:
		theFileName = inFileName
	else:
		theFileName = getResultFileName(inFileName)
	return theFileName.replace("BELLE_fit_results/BELLE_fit_","BELLE_fit_results/hessians/BELLE_COMA_")

def getResultFileName(inFileName):
	theEnd = inFileName.split("_")[~0]
	for fn in os.listdir("./build/BELLE_fit_results"):
		if fn.endswith(theEnd):
			return "./build/BELLE_fit_results/" + fn
	raise IOError("No hessian file for: '"+fn+"'")

def loadMatrix(inFileName, freeMap = None, fixMap = None):
	matrix = []
	with open(inFileName,'r') as inFile:
		for line in inFile.readlines():
			matrix.append([float(v) for v in line.split()])

	if freeMap is not None:
		hss     = cutFree(matrix,freeMap, mult = 2, fixMap = fixMap)
		matrix = []
		for line in hss:
			matrix.append(cutFree(line, freeMap, mult = 2, fixMap = fixMap))
	return np.asarray(matrix)

def getComaHist(resultFileName, fixMap, mass, conj = False):
	"""
	Creates a histogram for the covariance matrix
	"""
	freeMap         = getFreeMap(resultFileName)
	hessianFileName = getHessianFileName(resultFileName)
	COMAfileName    = getCOMAfileName(resultFileName)
	nBin            = getNbin(mass)
	if os.path.isfile(COMAfileName):
		print "Using COMA file"
		coma = loadMatrix(COMAfileName)
		coma = cutFree(coma, freeMap, mult = 2, fixMap = fixMap)
		COMA    = []
		for line in coma:
			COMA.append(cutFree(line,freeMap, mult = 2, fixMap = fixMap))

	elif os.path.isfile(hessianFileName):
		hessian = loadMatrix(hessianFileName)
		coma    = la.pinv(hessian)
		coma    = cutFree(coma,freeMap, mult = 2, fixMap = fixMap)
		COMA    = []
		for line in coma:
			COMA.append(cutFree(line,freeMap, mult = 2, fixMap = fixMap))
	else:
		print "No hessian file for '" + resultFileName + "' -> using unit matrix"
		dim = 0
		for i,f in enumerate(freeMap):
			if f:
				dim += freeNwaves[i]
		dim *= 2 # re and im
		coma = np.zeros(dim)
		COMA = [[0.]*dim for _ in xrange(dim)]
		for i in xrange(dim):
			COMA[i][i] = 1.

	hist = ROOT.TH2D('COMA_0_'+str(nBin),'COMA_0_'+str(nBin),len(COMA), 0.,1., len(COMA), 0.,1.)
	for i in xrange(len(coma)):
		for j in xrange(len(coma)):
			fakk = 1.
			if conj:
				fakk = (-1)**(i+j)
			hist.SetBinContent(i+1, j+1, fakk* COMA[i][j])
	return hist

def getLL(inFileName):
	lastPart = inFileName.split(os.sep)[~0]
	if not lastPart.startswith("BELLE_fit_"):
		raise ValueError("Could not getLL('"+inFileName+"')")
	return float(lastPart.split("_")[~1])

def parseResultFile(inFileName, conj = False):
	results = []
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			s = line.strip().replace(",","+").replace(")","j)").replace("+-","-")
			results.append(complex(s))
			if conj:
				results[~0] = results[~0].conjugate()
	if len(results) == 0:
		raise IOError("No amplitudes found in '" + inFileName + "'" )
	return results

def cutFree(arr, freeMap, fixMap, mult = 1):
	count  = 0
	retVal = []
	for i in xrange(9):
		if fixMap is not None:
			if fixMap[i]:
				count += freeNwaves[i]
				continue
		if freeMap[i]:
			for b in xrange(freeNwaves[i]):
				for mm in xrange(mult):
					retVal.append(arr[mult*count+mm])
				count += 1
		else:
			count += fixedNwaves[i]
	return retVal

def makeFullVector(cutVector, freeMap, double = True, addBG = True):
	fullVector = []
	count = 0
	for i in xrange(9):
		if freeMap[i]:
			for _ in xrange(freeNwaves[i]):
				fullVector.append(cutVector[count])
				count += 1
		else:
			for _ in xrange(fixedNwaves[i]):
				fullVector.append(0.)
	if double:
		nn = len(fullVector)
		for i in xrange(nn):
			fullVector.append(0.)
	if addBG:
		fullVector.append(0.)
	return np.asarray(fullVector, dtype = cutVector.dtype)

def getStartFileName(bestFileName):
	folder = "./build/BELLE_startValues/"
	fn     = bestFileName.split(os.sep)[~0]
	chunks = fn.split("_")
	svFN   = "BELLE_startValues"+chunks[2]+"_"+chunks[~0]
	return folder+svFN

def getInfoFileName(resultFileName):
	fileName = resultFileName.split(os.sep)[~0]
	chunks   = fileName.split("_")
	chunks.pop(3)
	fileName = "_".join(chunks)
	return "/home/iwsatlas1/fkrinner/Own_Software/ppppppppp/build/BELLE_fit_results/fitInfos/" + fileName

def getIsCopy(fn):
	infoFileName = getInfoFileName(fn)
	with open(infoFileName, 'r') as inFile:
		for line in inFile.readlines():
			if "-copy" in line:
				return True
	return False

def getBestFileName(freeMap,
	            resultFolder      = "/home/iwsatlas1/fkrinner/Own_Software/ppppppppp/build/BELLE_fit_results",
	            zeroMap           = None,
	            hasToBeAfterFile  = None,
	            hasToBeBeforeFile = None,
	            hasInfo           = [],
	            hasNotInfo        = []):
	freeString = ''
	for f in freeMap:
		if f:
			freeString += '1'
		else:
			freeString += '0'
	minLL  = INF
	bestFn = None

	afterDate = None
	if hasToBeAfterFile is not None:
		afterDate = os.path.getctime(hasToBeAfterFile)

	beforeDate = None
	if hasToBeBeforeFile is not None:
		beforeDate = os.path.getctime(hasToBeBeforeFile)

	allLL = []
	bfe = '.' + getBranchFileEnding()
	for fn in os.listdir(resultFolder):
		if not fn.endswith(bfe):
			continue

		if not freeString in fn:
			continue
		if afterDate is not None or beforeDate is not None:
			date = os.path.getctime(resultFolder + os.sep + fn)
			if date < afterDate and afterDate is not None:
#				print "File '" +fn+ "' is too early"
				continue
			if not date < beforeDate and beforeDate is not None:
#				print "File '" +fn+ "' is too late"
				continue

		infoFileName = getInfoFileName(fn)
		if os.path.isfile(infoFileName):
			print "Found '" + infoFileName + "'"
		else:
			continue
		with open(infoFileName, 'r') as inFile:
			fitInfo = [str(line.strip()) for line in inFile.readlines()]
		allInfo = True
		for info in fitInfo:
			if info in hasNotInfo:
				allInfo = False
				break
		for info in hasInfo:
			if not info in fitInfo:
				allInfo = False
				break
		if not allInfo:
			continue
		if zeroMap is not None:
			try:
				if not hasZeroWaves(resultFolder + os.sep + fn, zeroMap, freeMap):
					continue
			except IOError:
				continue
		ll = getLL(resultFolder + os.sep + fn)
		allLL.append(ll)
		if ll < minLL:
			minLL  = ll
			bestFn = resultFolder + os.sep + fn
	if bestFn is None:
		raise IOError("No result file found for freeString '" + freeString  + "'")
	return bestFn, allLL

def hasZeroWaves(fn, zeroMap, freeMap = None):
	if freeMap is None:
		freeMap = getFreeMap(fn)
	ampls = parseResultFile(fn)
	count = 0
	for w in xrange(9):
		if freeMap[w]:
			N = freeNwaves[w]
		else:
			N = fixedNwaves[w]
		for n in xrange(N):
			if not ampls[count] == 0.+0.j:
				if zeroMap[w]:
					return False
			count += 1
	return True

def getBinnings(freeMap, fixMap):
	retVal = []
	for i in xrange(9):
		if fixMap is not None:
			if fixMap[i]:
				continue
		if freeMap[i]:
			retVal.append(fDdb.allBinnings[i])
	return retVal

def getSectorNames(freeMap = [True]*9, fixMap = None):
	fullNames = ["KpiSright","KpiPright", "KpiDright", "KpiSwrong","KpiPwrong","KpiDwrong", "piPiS","piPiP","piPiD"]
	retVal    = []
	for i in xrange(9):
		if fixMap is not None:
			if fixMap[i]:
				continue
		if freeMap[i]:
			retVal.append(fullNames[i])
	return retVal

def makeAllWaveHists(sectorNames, prodAmps, binnings, norms, mass):
	if not len(sectorNames) == len(binnings):
		raise ValueError("Number of waves does not match: " + str(len(sectorNames)) + " vs. " + str(len(binnings)) + ".")
	nBins = 0
	for b in binnings:
		nBins += len(b) - 1
	if not nBins == len(prodAmps):
		raise ValueError("Total number of bins does not match: " + str(nBins) + " vs. " + str(len(prodAmps)) + ".")
	if not len(prodAmps) == len(norms):
		raise ValueError("Number of normalizations doe not match: " + str(len(prodAmps)) + " vs. " + str(len(norms)) + ".")
	hists = []
	cumulIndex = 0
	for i,sn in enumerate(sectorNames):
		b = binnings[i]
		hists += getSingleWaveHists(sn, b, cumulIndex, prodAmps, norms,mass)
		cumulIndex += len(b)-1
	return hists

def getSingleWaveHists(sectorName, binning, startIndex, prodAmps, norms, mass):
	npBinning = np.asarray(binning, dtype = np.float64)
	I    = ROOT.TH2D(sectorName + "_0_intens",sectorName + "_0_intens", 50, .5, 2.5, len(binning)-1, npBinning)
	re   = ROOT.TH2D(sectorName + "_0_real",sectorName + "_0_real", 50, .5, 2.5, len(binning)-1, npBinning)
	im   = ROOT.TH2D(sectorName + "_0_imag",sectorName + "_0_imag", 50, .5, 2.5, len(binning)-1, npBinning)
	nrm  = ROOT.TH2D(sectorName + "_0_norm",sectorName + "_0_norm", 50, .5, 2.5, len(binning)-1, npBinning)
	indx = ROOT.TH2D(sectorName + "_0_index",sectorName + "_0_index", 50, .5, 2.5, len(binning)-1, npBinning)
	nBin = getNbin(mass)
	for i in xrange(len(binning)-1):
		I.SetBinContent(   nBin+1, i+1, abs(prodAmps[startIndex + i])**2)
		re.SetBinContent(  nBin+1, i+1, prodAmps[startIndex + i].real)
		im.SetBinContent(  nBin+1, i+1, prodAmps[startIndex + i].imag)
		nrm.SetBinContent( nBin+1, i+1, norms[startIndex + i])
		indx.SetBinContent(nBin+1, i+1, startIndex + i)
	return I,re,im,nrm,indx

def makeZeroModeHists(sectorNames, zeroMode, eigenvalue,  binnings, nZero, mass, minAmount = 0.01):
	tot = np.dot(zeroMode,zeroMode)**.5
	parts = [0.]*len(binnings)
	count = 0
	for b,B in enumerate(binnings):
		for i in xrange(len(B)-1):
			parts[b] += zeroMode[count]**2
			count += 1
		parts[b] **= .5
		parts[b]  /= tot
	name  = ""
	nBins = 0
	for i,S in enumerate(sectorNames):
		if parts[i] >= minAmount:
			if not name == "":
				name += "<|>"
			name  += S
			nBins += len(binnings[i]) - 1
	nBin     = getNbin(mass)
	hist     = ROOT.TH2D("zero"+str(nZero)+"_0", name, 50, .5, 2.5, nBins, 0., 1.)
	eigen    = ROOT.TH1D("eigen"+str(nZero)+"_0", name, 50, .5, 2.5)
	count    = 0
	totCount = 0
	eigen.SetBinContent(nBin+1, eigenvalue)
	for i,B in enumerate(binnings):
		if parts[i] >= minAmount:
			for b in xrange(len(B)-1):
				hist.SetBinContent(nBin+1, count+1, zeroMode[totCount])
				count    += 1
				totCount += 1
		else:
			totCount += len(B) - 1
	scale = 0.
	for b in xrange(nBins):
		scale += hist.GetBinContent(nBin + 1, b+1)**2
	hist.Scale(1./scale**.5)
	return hist, eigen

def parseCmdLine(args, additionalOptions = []):
	argv = []
	for a in args:
		append = True
		for opt in additionalOptions:
			if a.startswith(opt):
				append = False
				break
		if append:
			argv.append(a)

	if len(argv) == 2 and len(argv[1]) > 1:
		freeString = argv[1]
		if not len(freeString) == 9:
			raise RuntimeError("Invalid length of freeString: '" + freeString + "'")
		freeMap    = []
		for c in freeString:
			if c == '0':
				freeMap.append(False)
			elif c == '1':
				freeMap.append(True)
			else:
				raise RuntimeError("Invalid cahracter in freeString '" + c + "'")
	elif len(argv) > 1:
		freeMap = [False]*9
		for i in xrange(1,len(argv)):
			n = int(argv[i])
			if n < 0 or n > 8:
				raise RuntimeError("Invalid wave number '" + argv[i] + "'")
			if freeMap[n]:
				print "parseCmdLine(...): WARNING: wave '" + argv[i] + "' freed twice"
			freeMap[n] = True
		freeString = ''
		for f in freeMap:
			if f:
				freeString += '1'
			else:
				freeString +='0'
	else:
		print "No arguments given. Exit."
		sys.exit(0)
	return freeMap, freeString

def main():
	bfe = '.' + getBranchFileEnding()

	conj     = False
	if "-conj" in sys.argv:
		conj = True

	makeZM   = False
	cutFreed = True
	loadConstants()
	print mD0, mKs, mPi
	freeMap, freeString = parseCmdLine(sys.argv, additionalOptions = ["-conj"])

	fixMap = [False, False,False, False, False, False, False, False, False]

#	bestFn  = getStartFileName(getBestFileName(freeMap))

#	firstFileWithoutPrior = "./build/BELLE_fit_results/BELLE_fit_111111111_-16956062.253189_1541618747.dat"
#	lastFileWithoutPrior  = "./build/BELLE_fit_results/BELLE_fit_111111111_-20638751.029967_1541673946.dat"

#	dataMarkerFile = "./build/BELLE_fit_results/BELLE_fit_000000000_-16927048.975653_1541596896.dat"

	bestFn, allLL = getBestFileName(freeMap, zeroMap = fixMap)
#	bestFn, allLL = getBestFileName(freeMap, zeroMap = zeroMap, hasToBeAfterFile = dataMarkerFile)

	isCopy = getIsCopy(bestFn)
	if isCopy:
		fixMap[0] = True # Remove [K,pi_right]_S
		fixMap[1] = True # Remove [K,pi_right]_P
		fixMap[2] = True # Remove [K,pi_right]_D


	allLL.sort()
#	for i in xrange(10):
#		print "bestLikelihoods",allLL[i]
	print bestFn

	print "-----------------------------"
	print "best result file: '"+bestFn+"'"
	print "info file is: '"+getInfoFileName(bestFn)+"'"
	print "created:",datetime.datetime.fromtimestamp(os.path.getctime(bestFn))
	print "-----------------------------"

	ps_fileName, ac_fileName = getIntegralFileNames(bestFn)
	integral                 = integralClass(ps_fileName)
	integral.fixMap          = fixMap

	if cutFreed:
		integral.cutFreed()
#	integral.cutFirstSector()

	prodAmps = parseResultFile(bestFn, conj = conj)
	prodAmps = cutFree(prodAmps, freeMap,fixMap = fixMap)
	binnings = getBinnings(freeMap, fixMap = fixMap)

	norms    = [integral.psMatrix[i,i].real for i in xrange(len(prodAmps))]

	integral.normalize()

	sectorNames = getSectorNames(freeMap, fixMap = fixMap)
	binning     = getBinnings(freeMap, fixMap = fixMap)

	hists = makeAllWaveHists(sectorNames, prodAmps, binnings, norms, mD0)

	integral.makeZMs(0.01)
	maxEV = np.max(integral.allEVs)

	if makeZM and cutFreed:
		for i,zeroMode in enumerate(integral.zms):
			zeroFileName = "./build/zeroModeFiles/"+freeString+"_"+str(i) + bfe
			fullZM = makeFullVector(zeroMode, freeMap, True, True)
			with open(zeroFileName, 'w') as outFile:
				for v in fullZM:
					outFile.write("("+str(v.real)+','+str(v.imag)+') ')
		return
	elif makeZM:
 # # # Use this after cutFirstSector
		print len(integral.zms),"::::::::::"
		for i,zeroMode in enumerate(integral.zms):
			zeroFileName = "./build/zeroModeFiles/"+freeString+"_"+str(i) + bfe
			with open(zeroFileName, 'w') as outFile:
				for v in zeroMode:
					outFile.write("("+str(v.real)+','+str(v.imag)+') ')
				for v in zeroMode:
					outFile.write("(0.,0.) ")
				outFile.write("(0.,0.) ")
		return

	titles = []
	for i, zeroMode in enumerate(integral.zms):
		zmh, eigh = makeZeroModeHists(sectorNames, zeroMode, integral.evs[i], binnings, i, mD0, minAmount = 0.01)
		title = zmh.GetTitle()
		while title in titles:
			title += "<|sec++>"
		titles.append(title)
		zmh.SetTitle(title)
		eigh.SetTitle(title)
		hists.append(zmh)
		hists.append(eigh)

	hists.append(getComaHist(bestFn, fixMap = fixMap, mass = mD0, conj = conj))

	with root_open("DdecayResults_"+freeString+"_"+getBranchFileEnding()+".root", "RECREATE"):
		for h in hists:
			h.Write()
	return

	c1 = ROOT.TCanvas()
	hist = ROOT.TH1D("h","h",len(integral.zms[0]),0.,1.)
	for i in xrange(len(integral.zms)):
		print len(integral.zms[i])
		for j in xrange(len(integral.zms[i])):
			hist.SetBinContent(j+1, integral.zms[i,j])
		hist.Draw()
		c1.Update()
		raw_input()

if __name__ == "__main__":
	main()
