#!/usr/bin/python
# cpIntegralFilesToDifferentBranch.py
# Created: 2018-11-30 16:32:55.658491
# Author: Fabian Krinner
import os, sys
from shutil import copyfile

def copyIntegrals(fromBranch, toBranch, folder, hasString = None):
	for fn in os.listdir(folder):
		if fn.endswith(fromBranch):
			if hasString is not None:
				if not hasString in fn:
					continue
			copyfile(folder + os.sep + fn, folder + os.sep + fn.replace(fromBranch, toBranch))

def main():
	folder     = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/build/integralFiles"
	fromBranch = "CPeff"
	toBranch   = "CPeff2"
#	hasString  = "Dalitz_model_efficiencyCP"
	hasString  = None
	copyIntegrals(fromBranch, toBranch, folder, hasString)

if __name__ == "__main__":
	main()
