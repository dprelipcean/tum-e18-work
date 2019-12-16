#!/usr/bin/python
# changeFitResults.py
# Created: 2018-12-04 14:31:37.624448
# Author: Fabian Krinner
import os, sys
import ROOT
from rootfabi import root_open
import numpy as np

from analyzeBELLE_version2 import parseResultFile

waves  = ["Kright", "Kwrong", "pipi"]
spins  = ['S','P','D']

def parseFile(inFileName):
	reals = []
	imags = []
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			vals = [float(v) for v in line.split()]
			reals.append(vals[0])
			imags.append(vals[2])
	return np.asarray(reals, dtype = np.float64), np.asarray(imags, dtype = np.float64)

def loadPlotPoints(corr = False):
	retVals   = []
	folder = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/zeroModeFitter_2D/BELLE_ampl_files/"
	typ    = "_data"
	if corr:
		typ = "_corr"

	for wave in waves:
		for L in spins:
			fileName = folder + wave + L + typ + '.argand'
			retVals.append(parseFile(fileName))
	return retVals

def parseRootMacro(inFileName):
	reals = []
	imags = []
	with open(inFileName, 'r') as inFile:
		inData     = False
		realsDone  = False
		for line in inFile.readlines():
			if "Double_t" in line and " = {" in line:
				inData = True
				continue
			if inData:
				val = float(line.replace(",","").replace("};",""))
				if realsDone:
					imags.append(val)
				else:
					reals.append(val)
				if "}" in line:
					if realsDone:
						return np.asarray(reals, dtype = np.float64), np.asarray(imags, dtype = np.float64)
					else:
						realsDone = True
						inData = False

def makeROOTfile(rootFileName, data):
	with root_open(rootFileName, "RECREATE") as outFile:
		count = 0
		for wave in waves:
			for L in spins:
				reals, imags = data[count]
				waveName = wave + L
				graph    = ROOT.TGraph(len(reals), reals, imags)
				graph.SetTitle(waveName)
				graph.SetName(waveName)
				graph.Write()
				count   += 1

def getMacroFileData(waveName, changed = False):
	folder   = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/changedValues/"
	fileName = folder + waveName
	if changed:
		fileName += "_changed.C"
	else:
		fileName += "_unchanged.C"
	return parseRootMacro(fileName)

def getAllChangedData(changed = False):
	retVals = []
	for wave in waves:
		for L in spins:
			waveName = wave + L
			retVals.append(getMacroFileData(waveName, changed))
	return retVals

def getFactors(data, corr, changedCorr):
	retVals = []
	count   = 0
	for wave in waves:
		for L in spins:
			realFakks = 1. + (changedCorr[count][0] - corr[count][0])/data[count][0]
			imagFakks = 1. + (changedCorr[count][1] - corr[count][1])/data[count][1]
			realFakks[np.isinf(realFakks)] = 1.
			imagFakks[np.isinf(imagFakks)] = 1.
			retVals.append((realFakks, imagFakks))
			count += 1
	return retVals

def modifyFitResult(inFileName, outFileName, factors):
	ampls = parseResultFile(inFileName)
	countA   = 0
	countW	 = 0
	for wave in waves:
		for L in spins:
			factorsRe = factors[countW][0]
			factorsIm = factors[countW][1]
			for a in range(len(factorsRe)):
				ampls[countA] = ampls[countA].real * factorsRe[a] + 1.j * ampls[countA].imag * factorsIm[a]
				countA += 1
			countW += 1
	countW	 = 0 # Just do it twice to hit the CP model
	for wave in waves:
		for L in spins:
			factorsRe = factors[countW][0]
			factorsIm = factors[countW][1]
			for a in range(len(factorsRe)):
				ampls[countA] = ampls[countA].real * factorsRe[a] + 1.j * ampls[countA].imag * factorsIm[a]
				countA += 1
			countW += 1
	if not countA == len(ampls) - 1:
		raise ValueError("Something's wrong")
	with open(outFileName, 'w') as outFile:
		for a in ampls:
			outFile.write("("+str(a.real)+","+str(a.imag)+")\n")


def main():
	data        = loadPlotPoints(False)
	corr        = loadPlotPoints(True)
	changedCorr = getAllChangedData(True)

	amplFileName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/build/BELLE_fit_results/BELLE_fit_111111111_-16951414.880645_1543656150.CPeff2"
	factors      = getFactors(data, corr, changedCorr)
	outFileName  = "./changedFileName.lol"
	modifyFitResult(amplFileName, outFileName, factors)
	
if __name__ == "__main__":
	main()
