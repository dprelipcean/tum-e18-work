#!/usr/bin/python
# chi2model.py
# Created: 2019-11-22 09:11:21.158962
# Author: Fabian Krinner
import os, sys
import numpy as np
import ROOT

from utils import isValidDalitzPoint, abs2

import bins

class dalitzChi2model:
	def __init__(self, mMother, fsMasses, bins, meshWidth, functions, nTermInBin = 10):
		"""
		@param mMother    mass of the mother particle
		@param fsMasses   list of final state masses: [m1, m2, m3] with m12**2 on the X-axis and m13**2 on the Y-axis
		@param bins       list of bins to use. Format [(m12_low, m12_high, m13_low, m13_high), ... ]
		@param meshWidth  distance of two points in one direction on the integral mesh
		@param functions  list of functions ( = pratialWaves)
		@param nTermInBin (estimated) number of non-zero terms in every bin ( = nonZeroAmpls**2)
		"""

		# Physics
		self.mMother    = mMother
		self.s          = mMother**2
		self.fsMasses   = fsMasses
		self.fsSquare   = [m**2 for m in fsMasses]

		# Binning
		self.nBins      = len(bins)
		self.bins       = bins

		# Functions
		self.functions  = functions
		self.nFunc      = len(functions)

		# Options
		self.meshWidth       = meshWidth
		self.nTermInBin = nTermInBin

		# Kinematic validity of bins
		self.nValid          = None
		self.validBins       = None # List of kinematically allowed bins

		# Data
		self.data            = None
		self.errors          = None

		# Integrals
		self.totalNmesh      = None
		self.nMeshBin        = None
		self.integrals       = None
		self.integralIndices = None
		self.norms           = None

		if not self.checkContiguous():
			print("dalitzChi2model.__init__(...): WARNING: The rewceived bins do not contiguously cover the Dalitz plot.")

		# Plotting
		self.binningX = None
		self.binningY = None
		self.makeMaximalBinning()

	def makeMaximalBinning(self):
		"""
		Creates all bin borders 
		"""
		binningX = []
		binningY = []
		for b in self.bins:
			allBorders = b.getAllBorders()
			for borders in allBorders:
				if not borders[0] in binningX:
					binningX.append(borders[0])
				if not borders[1] in binningX:
					binningX.append(borders[1])
				if not borders[2] in binningY:
					binningY.append(borders[2])
				if not borders[3] in binningY:
					binningY.append(borders[3])
		binningX.sort()
		binningY.sort()
		self.binningX = np.array(binningX, dtype = np.float64)
		self.binningY = np.array(binningY, dtype = np.float64)

	def getOverallFunctionIndex(self,f,g):
		"""
		Allows one-dimensional indexing of two functions
		@param f first function index
		@param g second funtcion index
		"""
		return f * self.nFunc + g

	def getSingleFunctionIndices(self, t):
		"""
		Allows one-dimensional indexing of two functions
		@param t combined function index
		"""
		f = int(t/self.nFunc)
		g = t - self.nFunc * f
		return f,g

	def findBin(self, m2_12, m2_13, breakFound = True):
		"""
		Find the bin index to a given point on the Dalitz Plot
		@param m2_12 invariant mass square of particles 12
		@param m2_13 invariant mass square of particles 13
		@param breekFound flag if break after first found bin (If False, returns the last valid bin found)
		"""
		nFound = 0
		nBin   = None
		for b, binn in enumerate(self.bins):
			if binn.contains(m2_12, m2_13):
				nFound += 1
				nBin    = b
				if breakFound:
					break
		return nBin, nFound

	def checkContiguous(self):
		"""
		Checks, if the received binning is contiguous
		"""
		nGrid   = int(self.s/self.meshWidth)
		linGrid = np.array([i * self.meshWidth for i in range(nGrid+1)])
		grid    = np.array(np.meshgrid(linGrid, linGrid)).T.reshape(-1,2)
		grid    = grid[isValidDalitzPoint(grid, self.s, self.fsSquare[0], self.fsSquare[1], self.fsSquare[2])]
		retVal  = True
		for i in range(grid.shape[0]):
			t,n = self.findBin(grid[i,0], grid[i,1], breakFound = False)
			if n == 0:
				print( "dalitzChi2model.checkContiguous(): WARNING: No valid bin found for valid grid point: ("+str(grid[i,0]) + ", " + str(grid[i,1]) + ").")
				retVal = False
			if n > 1:
                                print( "dalitzChi2model.checkContiguous(): WARNING: " + str(n) + " valid bins found for valid grid point: ("+str(grid[i,0]) + ", " + str(grid[i,1]) + ").")
				retVal = False
		return retVal
		
	def makeGrid(self, nBin):
		"""
		Makes the integral meshgrid for a 1-D bin index
		@param nBin bin index to make the grid for
		"""
		return self.bins[nBin].makeGrid(self.meshWidth)

#		binBorders = self.bins[nBin]

#		iMin = int(binBorders[0]/self.meshWidth) + 1
#		iMax = int(binBorders[1]/self.meshWidth) + 1

#		jMin = int(binBorders[2]/self.meshWidth) + 1
#		jMax = int(binBorders[3]/self.meshWidth) + 1
		
#		pX = np.array([i*self.meshWidth for i in range(iMin, iMax)])
#		pY = np.array([j*self.meshWidth for j in range(jMin, jMax)])

#		grid = np.array(np.meshgrid(pX,pY)).T.reshape(-1,2)
#		return grid

	def makeValidBins(self):
		"""
		Produces the list of valid bins i.e. the list of bins, where at least one integral point is valid
		"""
		self.validBins = np.zeros(self.nBins, dtype = bool)
		for b in range(self.nBins):
			grid  = self.makeGrid(b)
			valid = isValidDalitzPoint(grid, self.s, self.fsSquare[0], self.fsSquare[1], self.fsSquare[2])
			if np.sum(valid) > 0:
				self.validBins[b] = True
		for i,v in enumerate(self.validBins):
			if v:
				print i,
		print
		self.nValid = np.sum(self.validBins)

	def increaseIntegralSize(self):
		"""
		Increases the sorage swize of the intergal values for each bin by one
		"""
		print( "dalitzChi2model.increaseIntegralSize(): WARNING: The estimated number of non-zero ampliudes ("+str(self.nTermInBin)+") is too small. Increase by one.")
		for i in range(len(self.integrals)):
			self.integrals[i] = np.resize(self.integrals[i], self.nTermInBin+1)
			self.integralIndices[i] = np.resize(self.integralIndices[i], self.nTermInBin+1)
		self.nTermInBin += 1

	def makeIntegrals(self):
		"""
		Make the integral for all interference terms
		"""
		self.integrals       = []
		self.integralIndices = []
		self.totalNmesh      = 0
		self.nMeshBin        = np.zeros(self.nBins, dtype = int)
		self.norms = np.zeros(self.nFunc)
		for t in range(self.nBins):
			if not self.validBins[t]:
				continue
			indices, integral, nMesh = self.makeBinIntegrals(t)
			self.totalNmesh += nMesh
			self.nMeshBin[t] = nMesh
			self.integrals.append(integral)
			self.integralIndices.append(indices)
			for i, I in enumerate(indices):
				f,g = self.getSingleFunctionIndices(I)
				if f == g:
					self.norms[f] += integral[i]
		# Finally do the normalization
		for i,I in enumerate(self.integralIndices):
			for j,J in enumerate(I):
				f,g = self.getSingleFunctionIndices(J)
				self.integrals[i][j] /= (self.norms[f] * self.norms[g])**.5

	def makeBinIntegrals(self, nBin):
		"""
		Make the integral for all interference terms for a single bin
		@param nBin bin index to calculate the integral for
		"""
		if not self.validBins[nBin]:
			return None
		grid      = self.makeGrid(nBin)
		validGrid = isValidDalitzPoint(grid, self.s, self.fsSquare[0], self.fsSquare[1], self.fsSquare[2])
		nMesh     = np.sum(validGrid)
		grid      = grid[validGrid,:]
		fEvals    = np.array([f(grid) for f in self.functions])
		integral  = np.dot(np.conjugate(fEvals), fEvals.T)
		fIndices  = np.full(self.nTermInBin, -1, dtype = int)
		lintegral = np.zeros(self.nTermInBin, dtype = complex)
		count     = 0
		for f in range(self.nFunc):
			for g in range(f,self.nFunc):
				if integral[f,g] == 0.+0.j:
					continue
				if count == self.nTermInBin:
					fIndices.resize(self.nTermInBin+1) # Do not set the last value here, since it will be overwritten in any case (won't stay empty)
					lintegral.resize(self.nTermInBin+1)
					self.increaseIntegralSize()

				sfi              = self.getOverallFunctionIndex(f,g)
				fIndices[count]  = sfi
				lintegral[count] = integral[f,g]
				count += 1
		return fIndices, lintegral, nMesh

	def getIntentsity(self, prodAmps):
		"""
		Get the model intensity for all valid bins.
		@param prodAmps production amplotudes as complex numpy array
		"""
		retVal       = np.zeros(self.nValid, dtype = float)
		prodProdAmps = np.zeros(self.nFunc**2)
		for f in range(self.nFunc):
			for g in range(f,self.nFunc):
				sfi = self.getOverallFunctionIndex(f,g)
				prodProdAmps[sfi] = np.conjugate(prodAmps[f])*prodAmps[g]
				if f == g:
					prodProdAmps[sfi] /= 2 # Since below, we sum over 2*Re(...), the diagonal part need a factor of 1/2
		for t in range(self.nValid):
			val = 2*np.dot(prodProdAmps[self.integralIndices[t]], self.integrals[t]).real
			retVal[t] = val
		return retVal

	def eval(self, params):
		"""
		Evaluate the chi2 function
		@param params production amplutdes as real-values numpy array [re_0, im_0, re_1, ...]
		"""
		prodAmps = params[::2] + 1.j * params[1::2]
		model    = self.getIntentsity(prodAmps)
		return np.sum((model - self.data)**2 * self.invErrors)

	def loadData(self, m2_12, m2_13):
		"""
		Load the data points
		@param m2_12 numpy array for kinematic values for m2_12 of the events
		@param m2_13 numpy array for kinematic values for m2_13 of the events
		"""
		nDat = len(m2_12)
		if not len(m2_13) == nDat:
			raise ValueError("Number of data points does not match.")
		data = np.zeros(self.nBins, dtype = float)
		for d in range(nDat):
			t, n = self.findBin(m2_12[d], m2_13[d], True)
			if n == 0:
				print ("dalitzChi2model.loadData(...): WARNING: Could not find a valid bin for the data point: (" + str(m2_12[d]) + ", " + str(m2_13[d]) + ").")
				continue
			data[t] += 1.
		countValid = 0
		self.data  = np.zeros(self.nValid, dtype = float)
		for t in range(self.nBins):
			if not self.validBins[t]:
				if data[t] > 0:
					print("dalitzChi2model.loadData(...): WARNING: " + str(data[t]) + " events laying in invalid bin " + str(t) + ": " + str(self.bins[t]))
				continue
			self.data[countValid] = data[t]
			countValid += 1
		self.invErrors = 1./self.data**.5
		self.invErrors[np.isinf(self.invErrors)] = 0.

	def makeTheoHist(self, params):
		"""
		Produces a histogram for the theory
		@param params production amplitudes to be used (real-valued numpy array: [re_0, im_0, re_1, ...])
		"""
		prodAmps = params[::2] + 1.j* params[1::2]
		return self.makeHist(self.getIntentsity(prodAmps), "_intens_theo")

	def makeDataHist(self):
		"""
		Produces a histogram for the data points
		"""
		return self.makeHist(self.data, "_data")

	def makeHist(self, data, tag = ''):
		"""
		Produces a histogram from given data
		@param data numpy array of the data to be plotted
		@param tag  string to name the hsitogram after
		"""
		hist       = ROOT.TH2D("Dalitz" + tag,"Dalitz" + tag,len(self.binningX)-1, self.binningX,len(self.binningY)-1, self.binningY)
		delta      = self.meshWidth # offset from the actual border
		countValid = 0	
		for b in range(self.nBins):
			if self.validBins[b]:
				if not isinstance(self.bins[b], bins.rectangularBin):
					print "WARNING: The plotting is not set up to handle non-rectangular bins."

				borders = self.bins[b].getBorders()
				iMin    = hist.GetXaxis().FindBin(borders[0] + delta)
				iMax    = hist.GetXaxis().FindBin(borders[1] - delta)
				jMin    = hist.GetYaxis().FindBin(borders[2] + delta)
				jMax    = hist.GetYaxis().FindBin(borders[3] - delta)
				for i in range(iMin, iMax+1):
					for j in range(jMin, jMax+1):
						hist.SetBinContent(i,j,data[countValid]/self.nMeshBin[b])
				countValid += 1
		return hist

def uesen(g):
	return g[:,0]**2

def osh(g):
	return np.zeros(g.shape[0]) + 1.1

def heb(g):
	return np.full(g.shape[0],1.)

def hub(g):
	return g[:,0] + g[:,1]

def hob(g):
	return g[:,0] * g[:,1]

def halfBudalf(g):
	retVal = g[:,0] + g[:,1]
	retVal[g[:,0] < 1.] = 0.
	return retVal

def main():
	mDc = 1.86958
	mPi = 0.13957
	mKc = 0.493677

	fsMasses = [mKc, mPi, mPi]

	nBins = 20

	binningX = np.linspace(0.,mDc**2,nBins)
	binningY = np.linspace(0.,mDc**2,nBins)

	binList = []
	for i in range(len(binningX)-1):
		for j in range(len(binningY)-1):
			binList.append(bins.rectangularBin(binningX[i], binningX[i+1], binningY[j], binningY[j+1]))
#	fcns = [osh,heb,hub,hob,uesen, halfBudalf]
	fcns = [heb ,hub]
	c2model = dalitzChi2model(mDc, fsMasses, binList , 0.001, fcns, 10)
	c2model.makeValidBins()
	c2model.makeIntegrals()

	nData = 10000
	m2s = np.random.uniform((mPi+mKc)**2, (mDc - mPi)**2, 2*nData)
	m2s.shape = (nData, 2)
	valid = isValidDalitzPoint(m2s, mDc**2, fsMasses[0]**2, fsMasses[1]**2, fsMasses[2]**2)
	print "nGOOD =",np.sum(valid)

	m2s = m2s[valid,:]
	weights = 1. + hub(m2s)**2

	mxx = np.max(weights)

	weights -= np.random.random(weights.shape) * mxx

	m2s = m2s[weights > 0.,:]

	print "nDeWeight =",np.sum(weights > 0)

	c2model.loadData(m2s[:,0], m2s[:,1])

#	params = np.random.uniform(-1.,1.,2*len(fcns))
	params = np.zeros(2*len(fcns))		

	from iminuit import Minuit

	m = Minuit.from_array_func(c2model.eval, params, error=0.5)
	m.migrad()

	vals = np.array([m.values['x' + str(i)] for i in range(2*len(fcns))])

	h1 = c2model.makeTheoHist(vals)
	h1.Draw('col')
	raw_input()
	h2 = c2model.makeDataHist()
	h2.Draw('col')
	raw_input()
	print vals

	ntfrr =  ((vals[0] + 1.j*vals[1]) * (vals[2] - 1.j * vals[3]))**2
	print ntfrr/abs(ntfrr), "should be real (phase of pi/2)"

if __name__ == "__main__":
	main()

