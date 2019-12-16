#!/usr/bin/python
# getMasterDirectory.py
# Created: 2019-02-20 14:55:16.496690
# Author: Fabian Krinner
import os, sys

def getMasterDirectory():
	masterDirectFile = os.path.dirname(os.path.realpath(__file__)) + os.sep + "masterDirectory.h"
	with open(masterDirectFile, 'r') as inFile: # Parse the C++ file, so it can be changed globally at a single place.
				for line in inFile.readlines():
					if line.strip().startswith('//'): # ignore commented lines.
						continue
					if "const" in line and "std::string" in line and "masterDirectory" in line:
						chunks = line.split('"')
						return chunks[1]

def main():
	direct = getMasterDirectory()
	print direct

if __name__ == "__main__":
	main()
