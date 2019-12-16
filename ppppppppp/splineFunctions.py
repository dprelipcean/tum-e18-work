#!/usr/bin/python
# functions.py
# Created: 2019-09-17 15:54:55.856414
# Author: Fabian Krinner
import os, sys

class splineDeg:
	def __init__(self, minn,maxx,deg):
		self.min = minn
		self.max = maxx
		self.deg = deg

	def __call__(self, var, nDiff = 0, borderEval = 'all'):
		if nDiff > self.deg:
			return 0.
		x = (var-self.min)/(self.max-self.min)
		dxds = 1./(self.max-self.min)

		if not borderEval in ['all','upper','lower','none']:
			raise ValueError("Invalid option: " + borderEval)

		if x < 0. or x > 1.:
			return 0.
		if x == 1.:
			if borderEval == 'lower' or borderEval == 'none':
				return 0.
		if x == 0.:
			if borderEval == 'upper' or borderEval == 'none':
				return 0.

		if nDiff == 0:
			if self.deg == 0:
				retVal = 1.
			elif self.deg == 1:
				retVal = x -.5
			elif self.deg == 2:
				retVal = x**2 - x + 1./6.
			elif self.deg == 3:
				retVal = x**3 - 3./2.*x**2 + 3./5. *x -1./20
			elif self.deg == 4:
				retVal = x**4 -2. *x**3 + 9./7. *x**2 - 2./7.*x +1./70.
			elif self.deg == 5:
				retVal = x**5 - 5./2. *x**4 + 20./9. *x**3 - 5./6.*x**2 + 5./42.*x - 1./252.
		elif nDiff == 1:
			if self.deg == 1:
				retVal = 1.
			elif self.deg == 2:
				retVal = 2.*x - 1.
			elif self.deg == 3:
				retVal = 3.*x**2 - 3.*x + 3./5.
			elif self.deg == 4:
				retVal = 4.*x**3 - 6. *x**2 + 18./7. *x - 2./7. 
			elif self.deg == 5:
				retVal = 5.*x**4 - 10. *x**3 + 20./3. *x**2 - 5./3.*x + 5./42.
		elif nDiff == 2:
			if self.deg == 2:
				retVal = 2.
			elif self.deg == 3:
				retVal = 6.*x - 3.
			elif self.deg == 4:
				retVal = 12.*x**2 - 12. *x + 18./7.
			elif self.deg == 5:
				retVal = 20.*x**3 - 30. *x**2 + 40./3. *x - 5./3.
		elif nDiff == 3:
			if self.deg == 3:
				retVal = 6.
			elif self.deg == 4:
				retVal = 24.*x - 12. 
			elif self.deg == 5:
				retVal = 60.*x**2 - 60. *x + 40./3.
		elif nDiff == 4:
			if self.deg == 4:
				retVal = 24.
			elif self.deg == 5:
				retVal = 120.*x - 60. 
		elif nDiff == 5:
			if self.deg ==5:
				retVal = 120.
		return retVal * dxds**nDiff

def makeBinnedSplineFunctions(deg, binning):
	retVal = []
	for i in range(len(binning) -1):
		sMin = binning[i ]
		sMax = binning[i+1]
		for d in range(deg+2):
			retVal.append(splineDeg(sMin,sMax,d))
	return retVal

def main():
	binning  = [0.,1.,2.,3.,4.,5.,6.,7.,8.,9.]
	deg      = 5

	spliners = makeBinnedSplineFunctions(deg,binning)

	delta = 1.e-6
	xMin = 1.02
	xMax = 2.34
	xEvl = 1.12
	for deg in range(6):
		fcn = splineDeg(xMin, xMax, deg)
		for i in range(1,6):
			num = (fcn(xEvl + delta, i-1) - fcn(xEvl,i-1))/delta
			ana = fcn(xEvl,i)
			print deg,i,ana,num

if __name__ == "__main__":
	main()
