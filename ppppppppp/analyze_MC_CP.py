#!/usr/bin/python
# plot_r sult_differences.py
# Created: 2019-03-04 09:44:27.460364
# Author: Fabian Krinner
import os, sys

import ROOT
import numpy as np
import numpy.linalg as la

from utils import to_complex
from regularize_integral_matrix import parseMatrixFile, isHermitian, regulatrizeMatrix

from loadConstants import loadConstants

from analyzeBELLE_version2 import getSingleWaveHists, makeAllWaveHists, makeZeroModeHists, getNbin

from rootfabi import root_open

def uglyPlot(pts, fileName):
	c = ROOT.TCanvas()
	h = ROOT.TH1D('w','w',len(pts), 0.,1.)
	for i,p in enumerate(pts):
		h.SetBinContent(i+1, p)
	h.Draw()
	c.Print(fileName)

def sp(v1,v2):
	if not len(v1) == len(v2):
		raise ValueError("Size mismatch.")
	val = 0.
	for i in range(len(v1)):
		val += v1[i] * v2[i]
	return val

def normalize(v):
	N = sp(v,v)
	N**=.5
	for i in range(len(v)):
		v[i] /= N

def printSPmatrix(vecs1, vecs2):
	spm = np.zeros((len(vecs1), len(vecs2)))
	for i,v1 in enumerate(vecs1):
		for j,v2 in enumerate(vecs2):
			spm[i,j] = sp(v1,v2)
	print la.det(spm)

def getUnitComaHist(dim, mass, conj = False):
	"""
	Creates a histogram for the covariance matrix
	"""
	nBin            = getNbin(mass)
#	print "Bin number",nBin
#	print "No hessian file -> using unit matrix"
	dim *= 2 # re and im
	coma = np.zeros(dim)
	COMA = [[0.]*dim for _ in xrange(dim)]
	for i in xrange(dim):
		COMA[i][i] = 1.

	hist = ROOT.TH2D('COMA_0_'+str(nBin),'COMA_0_'+str(nBin),len(COMA), 0.,1., len(COMA), 0.,1.)
	for i in xrange(len(coma)):
		for j in xrange(len(coma)):
			fakk = 1.
			if conj:
				fakk = (-1)**(i+j)
			hist.SetBinContent(i+1, j+1, fakk* COMA[i][j])
	return hist

def hessToComa(hess):
	val,vec = la.eig(hess)
	for v in xrange(len(val)):
		if val[v] < -0.01:
			raise ValueError("Too small eigenvalue in HESS " + str(val[v]))
		elif val[v] < 0.:
			val[v] *= -1
		val[v] = 1./val[v]
	return np.dot(vec, np.dot(np.diag(val), np.transpose(vec)))

def getCOMAhist(inFileName, mass, numLim = 1.e-7):
        nBin         = getNbin(mass)
	COMAfileName = inFileName + ".COMA"
	HESS         = False
	HESSfileName = inFileName + ".HESSIAN"
	if os.path.isfile(HESSfileName):
		print "Using HESSIAN file '" + HESSfileName + "'."
		COMAfileName = HESSfileName
		HESS = True
	if not(os.path.isfile(COMAfileName)):
		raise IOError("COMA file name '" + COMAfileName + " does not exist")
	with open(COMAfileName, 'r') as inFile:
		COMA = []
		for line in inFile:
			vals = [float(v) for v in line.split()]
			if len(COMA) == 0:
				dim = len(vals)
			else:
				if not len(vals) == dim:
					raise ValueError("Number of values in the lines do not match")
			COMA.append(vals)
	if HESS:
		COMA = hessToComa(COMA)
	if not len(COMA) == dim:
		raise ValueError("COMA is not quadratic")
	hist = ROOT.TH2D('COMA_0_'+str(nBin),'COMA_0_'+str(nBin),dim, 0.,1., dim, 0.,1.)
	for i in xrange(dim):
		if COMA[i][i] < 0.:
			raise ValueError("Negative entry in the COMA: " + str(COMA[i][i]) + " at  index " + str(i) + ".")
		hist.SetBinContent(i+1, i+1,COMA[i][i])
		for j in xrange(i):
			if not abs(COMA[i][j]-COMA[j][i]) < numLim:
				raise ValueError("COMA not symmetric: " + str(COMA[i][j]) + " vs. " + str(COMA[j][i]) + " at indices " + str((i,j)) + ".")
			hist.SetBinContent(i+1, j+1, COMA[i][j])
			hist.SetBinContent(j+1, i+1, COMA[j][i])
	print("COMA hist sucessfully created from '" + COMAfileName + ".")
	return hist

class integralClass:
	def __init__(self, inFileName):
		self.fileName = inFileName
		self.matrix   = parseMatrixFile(inFileName)

		if not isHermitian(self.matrix):
			raise ValueError("matrix is not hermitian")

		self.dim      = len(self.matrix)
		self.fixMap   = None

		norm          = []
		for i in xrange(self.dim):
			nn = self.matrix[i,i]
			if nn == 0.:
				nn = 0.
			else:
				nn = 1./nn.real**.5
			norm.append(nn)
		self.norm   = np.asarray(norm)
		self.normed = False
		self.eigen  = False
		self.cut    = False
		self.numLim = 1.e-14
		self.zms    = None

	def indexToUnity(self, index):
		"""
		Sets the columns and lines to unity: matrix[i,j] = matrix[j,i] = delta[i,j] with i = @param index.
		This is to remove these indices from the zero-mode calculation.
		"""
		if index >= self.dim or index < 0:
			raise ValueError("Invalid index = " + str(index) + " received (max is " + str(self.dim) + ").")
		for i in xrange(self.dim):
			self.matrix[i,index] = 0.
			self.matrix[index,i] = 0.
		self.matrix[index, index]    = 1.
		self.norm[index] = 1.

	def makeZMs(self, maxZMev):
		val, vec = la.eig(self.matrix)
		self.allEVs = np.real(val)
		zms = []
		evs = []
		allEvs = []
		for i,v in enumerate(val):
			allEvs.append(v)
			if abs(v.imag) > self.numLim:
				raise ValueError("Non vanishing imaginary part of eigenvalue: " + str(v))
			if v.real < 0.:
				raise ValueError("Negative eigenvalue encountered: "+ str(v))
			if v.real < maxZMev:
				zms.append(np.real(vec[:,i]))
				evs.append(v.real)
#				print v.real,"as a zero-mode EV"
		allEvs.sort()

		self.zms       = np.asarray(zms)
		self.evs       = np.asarray(evs)
		self.zmsNormed = self.normed
		self.zmsCut    = self.cut
		self.eigen     = True

	def normalize(self):
		if self.normed:
			print "integralClass.normalize(): WARNING: Already normed... no nothing"
			return
		normmatrix = np.diag(self.norm)
		self.matrix = regulatrizeMatrix(np.dot(normmatrix, np.dot(self.matrix, normmatrix)), silent = True)
		if self.zms is not None:
			self.zms = np.dot(self.zms, normmatrix)
		self.normed = True

	def unnormalize(self):
		if not self.normed:
			print "integralClass.unnormalize(): WARNING: Not normed... no nothing"
			return
		normmatrix = np.diag(1./self.norm)
		self.matrix = regulatrizeMatrix(np.dot(normmatrix, np.dot(self.matrix, normmatrix)), silent = True)
#		if self.zms is not None:
#			self.zms = np.dot(self.zms, normmatrix)
		self.normed = False

def makePlots(vals, binnings):
	count  = 0
	retVal = []
	for b in binnings:
		hist = ROOT.TH1D("h_" + str(count), "h_" + str(count), len(b)-1, b)
		for i in xrange(len(b)-1):
			hist.SetBinContent(i+1, abs(vals[count + i])**2/(b[i+1]-b[i]))
		retVal.append(hist)
	return retVal

def loadResult(inFileName, conjugate = False):
	retVal = []
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			chunks = line.split()
			retVal.append(float(chunks[2]) + 1.j * float(chunks[3]))
	retArr = np.array(retVal, dtype = complex)
	if conjugate:
		return np.conjugate(retArr)
	return retArr

def loadCppResults(inFileName, conjugate = False):
	retVal = []
	with open(inFileName, 'r') as inFile:
		for line in inFile:
			retVal.append(to_complex(line))
	retArr = np.array(retVal, dtype = complex)
	if conjugate:
		return np.conjugate(retArr)
	return retArr

def getBinnings(inFileName):
	binnings = []
	binning  = []
	oldMmax  = None
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			chunks = line.split("|")
			mMin = float(chunks[1])**.5
			mMax = float(chunks[2])**.5
			if oldMmax is None:
				binning.append(mMin)
				binning.append(mMax)
				oldMmax = mMax
			elif mMin == oldMmax:
				binning.append(mMax)
				oldMmax = mMax
			elif mMin < oldMmax:
				binnings.append(np.array(binning))
				binning = [mMin, mMax]
				oldMmax = mMax
	binnings.append(np.array(binning))
	return binnings

def removeMode(val, mode_to_remove):
	norm = (np.sum(np.conjugate(mode_to_remove) * mode_to_remove))**.5
	sp   =  np.dot(np.conjugate(mode_to_remove), val)/norm
	val -= sp * mode_to_remove/norm

def getPriorBestFile(prior, dirr):
	for fnn in os.listdir(dirr):
		prior_fn = float(fnn.split('_')[4][5:])
		if prior == prior_fn:
			return dirr + os.sep + fnn
	raise IOError("No best file found for prior '" + str(prior) + "'")

def getChebyFileName(dirr):
	retVal = None
	for fn in os.listdir(dirr):
		if "dataLstring" in fn:
			continue
		if 'cheby' in fn:
			if retVal is not None:
				raise IOError("Found two 'cheby' file names '" + fn + "' and '" + retVal + "'.")
			retVal = fn
	if retVal is None:
		raise IOError("Could not find a cheby file in '" + dirr + "'")
	return dirr + os.sep + retVal

def getLletter(L):
	return "SPDFGH"[L]

def getLinLoutFileName(Ls_in, Ls_out, filters = []):
	folder = "./build/best_MC_fit/"
	prefix = ''
	for L in Ls_out:
		prefix += str(L)
	prefix += "_"
	suffix  = "_dataLstring"
	for L in Ls_in:
		suffix += getLletter(L)
	suffix += "_0_1_2.CP_MC"
	retVal = []
	for fn in os.listdir(folder):
		if fn.endswith(".COMA"):
			continue
		cntn = False
		for filt in filters:
			if not filt in fn:
				cntn = True
		if cntn:
			continue
		if fn.startswith(prefix) and fn.endswith(suffix):
			retVal.append(fn)
	if not len(retVal) == 1:
		print retVal
		minLL  = float('inf')
		bestFN = ''
		for fn in retVal:
			LL = float(fn.split("_")[2])
			if LL < minLL:
				minLL  = LL
				bestFN = fn
		print "WARNING: More than one file found!"
		return folder + bestFN

	return folder + retVal[0]


def main():
	conj  = True
#	prior = 10000.
	prior = .0
	Ls_in = []
	for i in sys.argv[1]:
		Ls_in.append(int(i))
	Ls_out = []
	for i in sys.argv[2]:
		Ls_out.append(int(i))


	constants        = loadConstants()
	dirr             = "./build/best_MC_fit/"
	inFileName       = getLinLoutFileName(Ls_in, Ls_out, ["jumpfft"])

#	inFileName       = "/home/iwsatlas1/fkrinner/Own_Software/ppppppppp/build/best_MC_fit/012_1564408430_-13098923.817761_jumpfft1.000000_dataLstringSPD_0_1_2.CP_MC"
#	inFileName       = "/home/iwsatlas1/fkrinner/Own_Software/ppppppppp/build/best_MC_fit/012_1558507036_-13098661.674814_jumpfft0.700000_dataLstringSPD_0_1_2.CP_MC"
	inFileName       = "/home/iwsatlas1/fkrinner/Own_Software/ppppppppp/build/best_MC_fit/012_1567075454_-13097571.295910_jumpfft1.000000_dataLstringSPD_0_1_2.CP_MCSTEP"

	print "Best file name '" + inFileName + "'"

	inFileName_bin   = inFileName
	integralFile     = "./build/integralFiles/ps_integral_CP_MCK-pi+pi+_L-012_0_1_2.CP_MC"
	integral         = integralClass(integralFile)

	binBorders       = [0,28,72,104]
#	for i in xrange(binBorders[2], binBorders[3]):
#		integral.indexToUnity(i)

	norms_all        = [integral.matrix[i,i].real for i in xrange(integral.dim)]
	norms            = []
	for L in Ls_out:
		for b in xrange(binBorders[L], binBorders[L+1]):
			norms.append(norms_all[b])

	for L in xrange(3):
		if L in Ls_out:
			continue
		for b in xrange(binBorders[L], binBorders[L+1]):
			integral.indexToUnity(b) # To remove them from the zero-modes.
			pass

	vals             = loadResult(    inFileName, conjugate = conj)
#	with open("012_amps", 'w') as outFile:
#		for v in vals:
#			outFile.write("("+str(v.real)+","+str(v.imag)+")\n")
#			print v
#	return
#	vals             = loadCppResults(inFileName, conjugate = conj)
	all_sector_names = ["[Kpi]S","[Kpi]P","[Kpi]D"]
	sector_names     = [all_sector_names[L] for L in Ls_out]
	all_binnings     = getBinnings(inFileName_bin)

	binnings         = [all_binnings[L] for L in Ls_out]
	import fine_binning as fine

	count            = 0
	nBin             = 34
	hists            = makeAllWaveHists(sector_names, vals, binnings, norms, constants["mDc"])

	integral.normalize()
	integral.makeZMs(0.02)
	integral.unnormalize()


	for z, zm in enumerate(integral.zms):
		ev = integral.evs[z]
#			for val in zm:
#				outFile.write(str(val) + "\n")
		h_zm, h_ev = makeZeroModeHists(sector_names, zm, ev, binnings, z, constants["mDc"], minAmount = 0.01)
		hists.append(h_zm)
		hists.append(h_ev)

	if True:
		coma_hist = getCOMAhist(inFileName, constants["mDc"])
#	except Exception as e:
#		print ("No valid covariance matrix found, due to:")
#		print (e)
#		print ("Using a unit matrix instead")
#		coma_hist = getUnitComaHist(len(vals), constants["mDc"])
	hists.append(coma_hist)

	with root_open("outFile.root", "RECREATE") as outRoot:
		for h in hists:
			h.Write()

if __name__ == "__main__":
	main()

