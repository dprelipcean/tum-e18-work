import sys
sys.path.append("/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/zeroModeFitter_2D")
import numpy as np
import numpy.linalg as la
from integral import toComplex, integral
import pyRootPwa
from rootfabi import root_open
import utils

binning3Pi = [.5]
while binning3Pi[-1] < 2.5:
	binning3Pi.append(binning3Pi[-1] + .04)
binning3Pi = np.asarray(binning3Pi, dtype = np.float64)

def loadBinning(inFileName):
	binning = []
	with open(inFileName, 'r') as inin:
		for line in inin.readlines():
			binning += [float(chunk)**.5 for chunk in line.split()]
	return np.asarray(binning)

def loadRealMatrixFile(fileName):
	matrix = []
	with open(fileName, 'r') as inin:
		for line in inin.readlines():
			vals = [float(chunk) for chunk in line.split()]
			matrix.append(vals)
	return np.asarray(matrix)

def loadAmplitudeFile(fileName, conjugate = False):
	ampls = []
	with open(fileName, 'r') as inin:
		for line in inin.readlines():
			aa = toComplex(line.strip())
			if conjugate:
				ampls.append(aa.conjugate())
			else:
				ampls.append(aa)
	return np.asarray(ampls, dtype = complex)

def loadN():
	ns = []
	with open('./build/largeErrorStudy_nDp.dat') as inin:
		for line in inin.readlines():
			ns += [int(v) for v in line.split()]
	if not len(ns) == 3:
		raise ValueError("Number of values does not match")
	return ns[0], ns[1], ns[2]

def projectOutPhaseDirection(matrix, amplitudes):
	phaseDirection = []
	norm           = 0.
	dim            = len(matrix)
	for a in amplitudes:
		phaseDirection.append(-a.imag)
		phaseDirection.append( a.real)
		norm += abs(a)**2
	norm **= .5
	for i in range(dim):
		phaseDirection[i] /= norm
	projector = np.zeros((dim,dim))
	for i in range(dim):
		projector[i,i] += 1.
		for j in range(dim):
			projector[i,j] -= phaseDirection[i] * phaseDirection[j]
	return np.dot(projector, np.dot(matrix, projector))			

def removeIndex(n, ampls, hessian):
	dim       = len(ampls)
	amplsNew  = []
	matrixNew = []
	for i in range(dim):
		if i == n:
			continue
		amplsNew.append(ampls[i])
		line1 = []
		line2 = []
		for j in range(dim):
			if j == n:
				continue
			line1.append(hessian[2*i  ,2*j  ])
			line1.append(hessian[2*i  ,2*j+1])
			line2.append(hessian[2*i+1,2*j  ])
			line2.append(hessian[2*i+1,2*j+1])
		matrixNew.append(line1)
		matrixNew.append(line2)
	return np.asarray(amplsNew, dtype = complex), np.asarray(matrixNew)

def main():
	nF0, nRho, nF2 = loadN()

	bin            = 34
#				Seed    # neg log like
#	seedAppendix   =  "_1498549364"
#
#	seedAppendix   =  "_1498570556" # -1.31183e+07
#	seedAppendix   =  "_1498570525" # -1.31178e+07
#	seedAppendix   =  "_1498570546" # -1.31174e+07
#	seedAppendix   =  "_1498570532" # -1.31184e+07
#	seedAppendix   =  "_1498570552" # -1.31156e+07

#	seedAppendix   =  "_1498652006" # -1.34472e+07
#	seedAppendix   =  "_1498652015" # -1.34515e+07
#	seedAppendix   =  "_1498652011" # -1.34499e+07
#	seedAppendix   =  "_1498652008" # -1.34471e+07
#	seedAppendix   =  "_1498652038" # -1.34464e+07

#	seedAppendix   =  "_1498815114" # -1.34505e+07
#	seedAppendix   =  "_1498815131" # -1.34466e+07
#	seedAppendix   =  "_1498815125" # -1.34459e+07
#	seedAppendix   =  "_1498815128" # -1.34462e+07
#	seedAppendix   =  "_1498815121" # -1.34475e+07

#	seedAppendix = "_1512582427"
	seedAppendix = "_1512637090"

	amplFileName   = "./build/largeErrorStudy_amplitudes"+seedAppendix+".dat"
	hessFileName   = "./build/largeErrorStudy_hessian"+seedAppendix+".dat"
	inteFileName   = "./build/largeErrorStudy_integral.dat"
	binFfileName   = "./build/largeErrorStudy_binningF0.dat"
	binRfileName   = "./build/largeErrorStudy_binningRho.dat"
	binF2fileName  = "./build/binningF2.dat"

	conjugate      = True

	binningF0      = loadBinning(binFfileName)
	binningRho     = loadBinning(binRfileName)
	binningF2      = loadBinning(binF2fileName)

	nBins          = 0
	eigenString    = ""
	sectors        = {}
	if nF2 > 1:
		dnb = len(binningF2)-1
		sectors["Dp[pi,pi]2++PiD"] = (nBins, nBins+dnb)
		nBins += dnb
		eigenString += "Dp[pi,pi]2++PiD"
	if nF0 > 1:
		dnb = len(binningF0)-1
		sectors["Dp[pi,pi]0++PiS"] = (nBins, nBins+dnb)
		nBins += dnb
		if len(eigenString) > 0:
			eigenString += "<|>"
		eigenString += "Dp[pi,pi]0++PiS"
	if nRho > 1:
		dnb = len(binningRho)-1
		sectors["Dp[pi,pi]1--PiP"] = (nBins, nBins+dnb)
		nBins += dnb		
		if len(eigenString) > 0:
			eigenString += "<|>"
		eigenString += "Dp[pi,pi]1--PiP"

	print nF0,nRho,nF2,nBins

	ampl    = loadAmplitudeFile(amplFileName, conjugate)
	hessian = loadRealMatrixFile(hessFileName)
	inte    = integral(inteFileName)

	if nRho == 1: # Remove fixed waves from the integral matrix
		inte.removeIndices(nF0 + nF2)
		ampl, hessian = removeIndex(nF0 + nF2, ampl, hessian)
	if nF0 == 1:
		inte.removeIndices(nF2)
		ampl, hessian = removeIndex(nF2, ampl, hessian)
	if nF2 == 1:
		inte.removeIndices(0)
		ampl, hessian = removeIndex(0, ampl, hessian)

	inte.norm()
	inte.eigen()
	norm = inte.norms

	svals, sev = inte.getSmallVectors(0.001)

	hasZeroMode = False
	if len(sev) > 0:
		hasZeroMode = True
	else:
		print "No zero-mode found"

	comaName   = "COMA_0_"+str(bin)
	histReName = "INTEGRAL_r_0_"+str(bin)
	histImName = "INTEGRAL_i_0_"+str(bin)
#	COMA     = la.inv(hessian)
	COMA     = utils.pinv(hessian, 1.e-4)
	if conjugate:
		for i in range(2*nBins):	
			for j in range(2*nBins):
				iNdex = i
				jNdex = j
				sign  = (-1)**(iNdex + jNdex)
				COMA[iNdex, jNdex] *= sign

#	COMA     = projectOutPhaseDirection(COMA, ampl)

	with root_open("DpPiPiPi_largeErrorStudy.root", "RECREATE") as outFile:
		COMAhist = pyRootPwa.ROOT.TH2D(comaName, comaName, 2*nBins, 0.,1., 2*nBins, 0.,1.)
		intReHist = pyRootPwa.ROOT.TH2D(histReName, histReName, nBins, 0.,1.,nBins,0.,1.)
		intImHist = pyRootPwa.ROOT.TH2D(histImName, histImName, nBins, 0.,1.,nBins,0.,1.)
		for i in range(nBins):
			for j in range(nBins):
					intReHist.SetBinContent(i+1,j+1,inte[i,j].real)
					intImHist.SetBinContent(i+1,j+1,inte[i,j].imag)
		intReHist.Write()
		intImHist.Write()

	
		for i in range(2*nBins):
#			COMAhist.SetBinContent(i+1, i+1, 1.) # first two are fixed isobar
			for j in range(2*nBins):
				iNdex = i
				jNdex = j
				COMAhist.SetBinContent(i+1, j+1, COMA[iNdex, jNdex]) # first two are fixed isobar
				pass
		COMAhist.Write()

		for i, zeroMode in enumerate(sev):
			eigenHist = pyRootPwa.ROOT.TH1D("eigen"+str(i)+"_0", eigenString, 50, .5 , 2.5)
			eigenHist.SetBinContent(bin+1, svals[i])
			eigenHist.Write()
	
			zeroHist = pyRootPwa.ROOT.TH2D("zero"+str(i)+"_0", eigenString, 50, .5 , 2.5, nBins, 0.,1.)
			for i in range(nBins):
				zeroHist.SetBinContent(bin+1, i+1, zeroMode[i].real)
			zeroHist.Write()

		for sector in sectors:
			if "1--" in sector:
				isobarBinning = binningRho
			elif "0++" in sector:
				isobarBinning = binningF0
			elif "2++" in sector:
				isobarBinning = binningF2
			else:
				"Sector not valid: '"+sector+"'"
			histIntens = pyRootPwa.ROOT.TH2D(sector+"_0_intens", sector+"_0_intens", len(binning3Pi)-1, binning3Pi, len(isobarBinning)-1, isobarBinning)
			histReal   = pyRootPwa.ROOT.TH2D(sector+"_0_real"  , sector+"_0_real"  , len(binning3Pi)-1, binning3Pi, len(isobarBinning)-1, isobarBinning)
			histImag   = pyRootPwa.ROOT.TH2D(sector+"_0_imag"  , sector+"_0_imag"  , len(binning3Pi)-1, binning3Pi, len(isobarBinning)-1, isobarBinning)
			histNorm   = pyRootPwa.ROOT.TH2D(sector+"_0_norm"  , sector+"_0_norm"  , len(binning3Pi)-1, binning3Pi, len(isobarBinning)-1, isobarBinning)
			histIndex  = pyRootPwa.ROOT.TH2D(sector+"_0_index" , sector+"_0_index" , len(binning3Pi)-1, binning3Pi, len(isobarBinning)-1, isobarBinning)
			for b,B in enumerate(range(sectors[sector][0],sectors[sector][1])):
				histIntens.SetBinContent(bin+1, b+1, abs(ampl[B])**2)
				histReal.SetBinContent(  bin+1, b+1, ampl[B].real)
				histImag.SetBinContent(  bin+1, b+1, ampl[B].imag)
				histNorm.SetBinContent(  bin+1, b+1, norm[B])
				histIndex.SetBinContent( bin+1, b+1, B)
			histIntens.Write()
			histReal.Write()
			histImag.Write()
			histNorm.Write()
			histIndex.Write()


if __name__ == "__main__":
	main()
