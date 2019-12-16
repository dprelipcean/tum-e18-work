#!/usr/bin/python
# plotMCdata.py
# Created: 2018-12-13 17:20:16.780501
# Author: Fabian Krinner
import os, sys
import ROOT

def loadFile(inFileName, indices):
	retVal = []
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			chunks = line.split()
			retVal.append([float(chunks[i]) for i in indices])
	return retVal

def main():
	inFileName = "./build2/MC_data"
	data       = loadFile(inFileName, [1,2])
	hist = ROOT.TH2D("dalitz","dalitz", 300, 0.,3.,300,0.,3.)
	for p in data:
		hist.Fill(p[0], p[1])
	hist.Draw("col")
	raw_input()

if __name__ == "__main__":
	main()
