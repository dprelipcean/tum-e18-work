#!/usr/bin/python
# loadConstants.py
# Created: 2019-02-20 15:03:50.928316
# Author: Fabian Krinner
import os, sys

def loadConstants(constantFileName = os.path.dirname(os.path.realpath(__file__)) + os.sep + "constants.h"):
	constMap = {}
	with open(constantFileName, 'r') as inFile:
		for line in inFile.readlines():
			if line.startswith('#'):
				continue
			chunks = line.split()
			name   = chunks[2]
			val    = chunks[4].strip()
			constMap[name] = float(val.replace(";",""))
	return constMap

def main():
	pass

if __name__ == "__main__":
	main()
