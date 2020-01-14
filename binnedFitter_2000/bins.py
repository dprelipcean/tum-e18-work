#!/usr/bin/python
# bins.py
# Created: 2020-01-03 09:10:40.360767
# Author: Fabian Krinner
import os, sys

import numpy as np

class rectangularBin:
	def __init__(self, xMin, xMax, yMin, yMax):
		self.borders = [xMin, xMax, yMin, yMax]

	def getAllBorders(self):
		return [self.borders]

	def getBorders(self):
		return self.borders

	def contains(self, x, y):
		return x > self.borders[0] and x <= self.borders[1] and y > self.borders[2] and y <= self.borders[3]

	def makeGrid(self, meshWidth):

		iMin = int(self.borders[0]/meshWidth) + 1
		iMax = int(self.borders[1]/meshWidth) + 1

		jMin = int(self.borders[2]/meshWidth) + 1
		jMax = int(self.borders[3]/meshWidth) + 1
		
		pX = np.array([i*meshWidth for i in range(iMin, iMax)])
		pY = np.array([j*meshWidth for j in range(jMin, jMax)])

		grid = np.array(np.meshgrid(pX,pY)).T.reshape(-1,2)
		return grid

def main():
	pass

if __name__ == "__main__":
	main()
