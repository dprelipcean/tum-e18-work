#!/usr/bin/python
# utils.py
# Created: 2019-11-22 09:12:34.548459
# Author: Fabian Krinner
import os, sys
import numpy as np

def abs2(x):
	return x.real**2 + x.imag**2

def q2(M2, m12, m22):
	return (M2**2 + m12**2 + m22**2 - 2.*(M2*m12 + M2*m22 + m12 * m22))/4/M2

def isValidDalitzPoint(dalitzPoints,s, m12,m22,m32):
	s23     = s + m12 + m22 + m32 - dalitzPoints[:,0] - dalitzPoints[:,1];

	retVal  = s23 >= (m22**.5+m32**.5)**2
	retVal &= s23 <= (s**.5-m12**.5)**2
	p12 = q2(s, s23, m12)
	p22 = q2(s, dalitzPoints[:,1], m22)
	p32 = q2(s, dalitzPoints[:,0], m32)
	q12 = q2(dalitzPoints[:,0], m12, m22)
	q13 = q2(dalitzPoints[:,1], m12, m32)
	q23 = q2(s23, m22, m32)

	retVal &= p12 >= 0.
	retVal &= p22 >= 0.
	retVal &= p32 >= 0.
	retVal &= q12 >= 0.
	retVal &= q13 >= 0.
	retVal &= q23 >= 0.

	p1 = abs(p12)**.5 # can take the abs() here, since p < 0. is already set to false in retVal
	p2 = abs(p22)**.5
	p3 = abs(p32)**.5

	retVal &= p1 + p2 > p3
	retVal &= p1 + p3 > p2
	retVal &= p2 + p3 > p1

	return retVal.flatten()

def main():
	import ROOT
	mDc = 1.86958
	mPi = 0.13957
	mKc = 0.493677
	nBins = 1000
	hist = ROOT.TH2D("f","f",nBins,0., mDc**2, nBins,0.,mDc**2)
	binCenters = []
	for i in range(nBins):
		for j in range(nBins):
			binCenters.append([hist.GetXaxis().GetBinCenter(i+1), hist.GetYaxis().GetBinCenter(j+1)])
	binCenters = np.array(binCenters)

	valid        = isValidDalitzPoint(binCenters, mDc**2, mKc**2, mPi**2, mPi**2)
	count = 0
	for v in valid:
		if v:
			count += 1
	print 'sdfsfa', count
	centersToSet = binCenters[valid,:]
	print centersToSet.shape
	for cts in centersToSet:
		i = hist.GetXaxis().FindBin(cts[0])
		j = hist.GetYaxis().FindBin(cts[1])
		hist.SetBinContent(i,j,1.)
	hist.Draw("col")
	raw_input()


if __name__ == "__main__":
	main()
