#!/usr/bin/python
# summarizeLikelihood.py
# Created: 2018-11-09 10:09:40.113776
# Author: Fabian Krinner
import os, sys
import ROOT
from analyzeBELLE_version2 import getLL
from getBranchFileEnding import getBranchFileEnding
from getMasterDirectory import getMasterDirectory

def print_stats(ll, name = None):
	dim  = len(ll)
	mean = 0.
	for val in ll:
		mean += val
	mean /= dim
	err = 0.
	for val in ll:
		err += (val - mean)**2
	err **= .5
	err  /= dim - 1

	print "-------------------"
	if name is not None:
		print name
	print "dim",dim
	print "min",min(ll)
	print "mean",mean,"+-",err
	print "max",max(ll)
	print "-------------------"
	return min(ll),mean,err,max(ll)

def main():
	bfe = getBranchFileEnding()

	allLL   = []
	FFTbyLL = {}
	LLbyFFT = {}

	priorVal = 0.

	inFileName = "./build_CP_MC/summary.dat"
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			chunks = line.split()
			fft    = float(chunks[0])
			LL     = float(chunks[1])
			prior  = float(chunks[2])
			if not prior == priorVal:
				continue	
			allLL.append(LL)
			if not LL in FFTbyLL:
				FFTbyLL[LL] = []
			FFTbyLL[LL].append(fft)
			if not fft in LLbyFFT:
				LLbyFFT[fft] = []
			LLbyFFT[fft].append(LL)
			

#	singleFolder = getMasterDirectory() + "/MC_fits"
#	for fn in os.listdir(singleFolder):
#		if not "prior1" in fn:
#			continue
#		if fn.endswith(bfe):
#			chunks = fn.split("_")
#			LL = float(chunks[2])
#			fft =float(chunks[3][3:])
#			if not LL in FFTbyLL:
#				FFTbyLL[LL] = []
#			FFTbyLL[LL].append(fft)
#			if not fft in LLbyFFT:
#				LLbyFFT[fft] = []
#			LLbyFFT[fft].append(LL)
#			allLL.append(LL)

#	keys = []
#	for key in LLbyFFT:
#		keys.append(key)
#	keys.sort()

#	histMin  = ROOT.TH1D("min" ,"min" ,10, 0., 1.)
#	histMean = ROOT.TH1D("mean","mean",10, 0., 1.)
#	histMax  = ROOT.TH1D("max" ,"max" ,10, 0., 1.)

#	for i,coeff in enumerate(keys):
#		mi, me, er, ma = print_stats(LLbyFFT[coeff],str(coeff))
#		histMin.SetBinContent(i+1, mi)
#		histMean.SetBinContent(i+1, me)
#		histMean.SetBinError(i+1, er)
#		histMax.SetBinContent(i+1, ma)

#	histMax.SetMinimum(histMin.GetMinimum())
#	histMax.Draw()
#	histMin.Draw("SAME")
#	histMean.Draw("SAME")
#	raw_input()

#	return


	allLL.sort()
	minLL = allLL[0]

	nBins    = 101
	minLL = min(allLL)
	hist     = ROOT.TH1D("LL-dist","LL-dist", nBins, -1, 100)

	for LL in allLL:
		hist.Fill(LL - minLL)
	c1 = ROOT.TCanvas()
	hist.Draw()
	c1.Print('LL_summary.png')

	allLL.sort()
	for i in xrange(5):
		print i, allLL[i], FFTbyLL[allLL[i]]
	raw_input()

if __name__ == "__main__":
	main()
