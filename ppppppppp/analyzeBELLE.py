#!/usr/bin/python
# analyzeBELLE.py
# Created: 2018-06-20 09:11:19.799418
# Author: Fabian Krinner
import os, sys
import numpy as np
import numpy.linalg as la
import freeDdecayBinnings as fDdb
import scipy.optimize
import ROOT

from regularize_integral_matrix import parseMatrixFile

def loadConstants(constantFileName = "./constants.h"):
	with open(constantFileName, 'r') as inFile:
		for line in inFile.readlines():
			if line.startswith('#'):
				continue
			chunks = line.split()
			name   = chunks[2]
			val    = chunks[4].strip()
			exec "global "+name+";"+name +'='+val

def getBestSeed(folder, piSign = "_pi0_", prefix = "BELLEfit_"):
	bestLL = float('inf')
	lls = []
	for fn in os.listdir(folder):
#	#	if not "fabili" in fn:
#	#		continue
		if not fn.startswith(prefix):
			continue
		if not piSign in fn:
			continue
		llString = fn.split("ll")[-1][:-4].split('_')[0]
		ll = float(llString)
		lls.append(ll)
		if ll < bestLL:
			bestLLstring = llString
			bestLL = ll
			bestFn = fn
	print  bestFn,"is the best fit"
	return bestFn.split("seed")[1].split("_")[0],bestFn, lls

def readResultFile(inFileName, cutLast = False):
	results = []
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			s = line.strip().replace(",","+").replace(")","j)").replace("+-","-")
			results.append(complex(s))
	if cutLast:
		results = results[:~0]
	return np.asarray(results, dtype = complex)

def normalize(matrix):
	norm = np.zeros((len(matrix)), dtype = matrix.dtype)
	for i in range(len(matrix)):
		norm[i] = matrix[i,i]**.5
	for i in range(len(matrix)):
		for j in range(len(matrix)):
			matrix[i,j] /= norm[i]*norm[j]
	return norm

def applyNorm(matrix, norm):
	for i, iN in enumerate(norm):
		for j, jN in enumerate(norm):
			matrix[i,j]/= iN*jN

def getZeroModes(matrix, maxVal):
	val, vec = la.eig(matrix)
	iVecs    = []
	valPy = []
	for i,v in enumerate(val):
		if not v.imag == 0.:
			raise ValueError("Complex eigenvalue: "+str(v))
		if v.real < 0.:
			raise ValueError("Negative eigenvalue: "+str(v))
		if v.real < maxVal:
			iVecs.append(i)
		valPy.append(v.real)
	valPy.sort()
	with open("vals.dat",'w') as outFile:
		for i,v in enumerate(valPy):
			outFile.write(str(i) + ' ' + str(v) + '\n')
	vecs = []
	for i in iVecs:
		vector = []
		for j in range(len(matrix)):
			vector.append(vec[j,i].real)
		vecs.append(np.asarray(vector))
	return vecs

def simplestZeroModeFixing(ampls, zms, comaInv, binRange = None):
	if binRange is None:
		binRange = range(len(ampls))
	def chi2(params):
		cpls = [params[2*z] + 1.j*params[2*z+1] for z in range(len(zms))]
		amplCorr = np.zeros(2*len(ampls))
		for i in binRange:
			a = ampls[i]
			for z, zm in enumerate(zms):
				a += cpls[z] * zm[i]
			amplCorr[2*i  ] = a.real
			amplCorr[2*i+1] = a.imag
		c2 = np.dot(amplCorr, np.dot(comaInv, amplCorr))
#		c2 = np.dot(amplCorr, amplCorr)
		return c2
	res = scipy.optimize.minimize(chi2, [0.]*2*len(zms))
	return res.x

def getABC(ampls, zms, comaInv, binRange = None):
	if binRange is None:
		binRange = range(len(ampls))
	dim      = 2*len(ampls)
	ncf      = 2*len(zms)
	linAmpls = np.zeros((dim))
	linZms   = np.zeros((dim, ncf))
	for i in range(dim/2):
		if not i in binRange:
			continue
		linAmpls[2*i  ] = ampls[i].real
		linAmpls[2*i+1] = ampls[i].imag
		for j in range(ncf/2):
			linZms[2*i  ,2*j  ] = zms[j][i].real
			linZms[2*i  ,2*j+1] =-zms[j][i].imag
			linZms[2*i+1,2*j  ] = zms[j][i].imag
			linZms[2*i+1,2*j+1] = zms[j][i].real
	C = np.dot(linAmpls, np.dot(comaInv, linAmpls))
	B = 2*np.dot(np.dot(linAmpls, comaInv), linZms)
	A = np.dot(linZms.T,np.dot(comaInv, linZms))

	return A,B,C

def getModel(bngs, borders, conjugateMap = [False, False, False, False, False, False, False, False, False], useMap = [True, True, True, True, True, True, True, True, True], norms = None, masses = None, widths = None):
	functions = []
#	masses    = [0.824, 0.895, 1.430, 0.824, 0.895,  1.430, 0.500, 0.77526, 1.2755]
#	widths    = [0.478, 0.047, 0.109, 0.478, 0.047 , 0.109, 0.500, 0.1491 , 0.1867]

	if masses is None:
		masses    = [0.824, 0.895, 1.430, 0.824, 0.895,  1.430, 0.980, 0.77526, 1.2755]
	if widths is None:
		widths    = [0.478, 0.047, 0.109, 0.478, 0.047 , 0.109, 0.110, 0.1491 , 0.1867]

	for i in range(9):
		if not useMap[i]:
			continue
		m0 = masses[i]
		G0 = widths[i]

		num = m0*G0
		den = m0**2 - 1.j*m0*G0

		vals = np.zeros((borders[-1]), dtype = complex)

		for j in range(len(bngs[i])-1): # sixth isobar is pipi f0
			m  = (bngs[i][j] + bngs[i][j+1])/2
			vals[borders[i]+j] = num/(den - m**2)
			if conjugateMap[i]:
				vals[borders[i]+j] = vals[borders[i]+j].conjugate()
			if norms[i] is not None:
				vals[borders[i]+j] *= norms[borders[i]+j]
		functions.append(vals)
	return functions

def getEquidistandBinning(mMother, mDaughter1, mDaughter2, mBatchelor, binWidth):
	m = mDaughter1 + mDaughter2
	binning = [m]
	while m < mMother - mBatchelor:
		m += binWidth
		binning.append(m)
	return np.asarray(binning, dtype = np.float64)

def readHessianFile(inFileName, cutLast = False):
	hessian = []
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			vals = [float(v) for v in line.split()]
			for v in vals:
				if np.isnan(v):
					raise ValueError("NaN encountered in parsing the hessian")
			if cutLast:
				vals = vals[:~1] # cut 2 because of re and im
			hessian.append(vals)
	if cutLast:
		hessian = hessian[:~1] # cut 2 because of re and im
	return np.asarray(hessian)

def loadComaInv(inFileName, cutLast = False):
	hessian = readHessianFile(inFileName, cutLast = cutLast)
	return hessian

def removeZeroModeFromComa(coma, zeroModes, acv = None): # acv = Artificial coma value
	"""
	Projects the coma such, that all directions of zero modes are removed. 
	(i.e. no uncertainties remain in these directions if an acv (artificial 
	coma value) is given. The zero mode directions are then added again with
	this eigenvalue. This forces the coma to have eigenvectors in ZM 
	direction, which is a very useful property when fixing the corresponding
	ambiguities
	"""
	dim = len(coma)
	transformatrionMatrix = np.identity(dim)
	print len(zeroModes[0]),"'''''''''''''''''''''''''''''''''''''''''", dim
	for zm in zeroModes:
		for i in range(dim/2):
			for j in range(dim/2):
				transformatrionMatrix[2*i  ,2*j  ] -= zm[i]*zm[j]
				transformatrionMatrix[2*i+1,2*j+1] -= zm[i]*zm[j]
	transformed = np.dot(transformatrionMatrix, np.dot(coma,transformatrionMatrix))
	if acv is not None:
		for zm in zeroModes:
			for i in range(dim/2):
				for j in range(dim/2):
					transformed[2*i  ,2*j  ] += acv * zm[i] * zm[j]
					transformed[2*i+1,2*j+1] += acv * zm[i] * zm[j]
	return transformed

class chi2:
	def __init__(self, ampls, zms, params, bngs, norms, comaInv):
		self.ampls        = ampls
		self.zms          = zms
		self.params       = params
		self.binRange     = None
		self.bngs         = bngs
		self.norms        = norms
		self.useMap       = [True]*9
		self.conjugateMap = [False]*9
		self.comaInv      = comaInv
		self.initBorders()

	def initBorders(self):
		self.borders = [0]
		for b in self.bngs:
			self.borders.append(self.borders[~0] + len(b) -1)

	def __call__(self, params):
		funcs = getModel(self.bngs, self.borders, self.conjugateMap, self.useMap, self.norms, params[:9], params[9:])
		A,B,C = getABC(self.ampls, self.zms + funcs, self.comaInv, self.binRange)
		params    = -np.dot(la.inv(A+A.T), B)
		c2 = np.dot(params,np.dot(A, params)) + np.dot(B,params) + C
		return c2

	def getCpls(self):
		funcs = getModel(self.bngs, self.borders, self.conjugateMap, self.useMap, self.norms, self.params[:9], self.params[9:])
		A,B,C = getABC(self.ampls, self.zms + funcs, self.comaInv, self.binRange)
		params    = -np.dot(la.inv(A+A.T), B)
		return params

def main():
	loadConstants()
	print mD0, mKs, mPi

	freeString = '111111111'

	resultFolder   = "./build/BELLE_fit_results/"


	bestFileName             = resultFolder + "BELLE_fit_111111111_-23163380.140807_1538140262.dat"

	intFolder = "./build/integralFiles/"

	ps_name  = intFolder + "ac_integral_model_"+freeString+"_regular.dat"
	ac_name  = intFolder + "ps_integral_model_"+freeString+"_regular.dat"


	result   = readResultFile(bestFileName, cutLast = True)

	nTot = 2*(len(fDdb.binningKpiS) + len(fDdb.binningKpiP) + len(fDdb.binningKpiD) - 3) + (len(fDdb.binningPiPiS) + len(fDdb.binningPiPiP) + len(fDdb.binningPiPiD) -3)
	result = result[:nTot]

	comaInv = np.identity(2*nTot)


	acMatrix = parseMatrixFile(ac_name)
	psMatrix = parseMatrixFile(ps_name)

	norm = normalize(psMatrix)
	applyNorm(acMatrix, norm)

	binCount  = 0
	bngs      = [fDdb.binningKpiS,fDdb.binningKpiP,fDdb.binningKpiD,fDdb.binningKpiS,fDdb.binningKpiP,fDdb.binningKpiD,fDdb.binningPiPiS,fDdb.binningPiPiP,fDdb.binningPiPiD]
	for i in range(len(bngs)):
		bngs[i] = np.asarray(bngs[i], dtype = np.float64)
	

	borders   = [0]
	for bng in bngs:
		borders.append(borders[-1] + len(bng) - 1)

	zms       = getZeroModes(psMatrix, 0.01)
	writZeroModes = False
	if writZeroModes:
		for z,zm in enumerate(zms):
			with open("./build/zm_"+str(z)+".dat", 'w') as outFile:
				for v in zm:
					outFile.write("("+str(v)+",.0)\n")

	print len(zms), "zero modes found"

	coma      = removeZeroModeFromComa(la.pinv(comaInv), zms, 1.)
	comaInv   = la.pinv(coma)

	conjugateMap = [False]*9

	resonanceParams = [1.4405,.8937,1.4256,1.4405,.8937,1.4256,0.98, .77526, 1.2751,0.1926,.0472,.0985,0.1926,.0472,.0985,.1,.1478,.1851]

	functions = getModel(bngs, borders, norms = norm, conjugateMap = conjugateMap, masses = resonanceParams[:9], widths = resonanceParams[9:])
#	print functions	
	useMap = [True]*9
	newFunctions = []
	for i,f in enumerate(functions):
		if useMap[i]:
			newFunctions.append(f)
	functions = newFunctions

	binRange  = range(borders[~0])

	laClasse  = chi2(result, zms, resonanceParams, bngs, norm, comaInv)
	laClasse.conjugateMap = conjugateMap
#	laClasse.useMap = useMap
	params = laClasse.getCpls()

	params = [0.]*len(params)

	cpls      = [params[2*z] + 1.j*params[2*z+1] for z in range(len(zms))]
	fCpls     = [params[2*(len(zms)+f)] + 1.j*params[2*(len(zms)+f)+1] for f in range(len(functions))]
	amplCorr  = np.zeros((len(result)), dtype = result.dtype)
	for i in range(len(result)):
		a = result[i]
		for z,zm in enumerate(zms):
			a  += cpls[z]*zm[i]
		amplCorr[i] = a

	c1        = ROOT.TCanvas()

	histNames = iter(["KpiS","KpiP","KpiD","KpiS","KpiP","KpiD","piPiS","piPiP","piPiD"])

	for w,bng in enumerate(bngs):
		n     = len(bng) - 1
		re    = np.zeros((n), dtype = np.float64)		
		im    = np.zeros((n), dtype = np.float64)

		reerr = np.zeros((n), dtype = np.float64)		
		imerr = np.zeros((n), dtype = np.float64)

		hist  = ROOT.TH1D('h','h',n,bng)
		hist2 = ROOT.TH1D('h','h',n,bng)

		hist_re  = ROOT.TH1D('h','h',n,bng)
		hist2_re = ROOT.TH1D('h','h',n,bng)

		hist_im  = ROOT.TH1D('h','h',n,bng)
		hist2_im = ROOT.TH1D('h','h',n,bng)

		hist.SetTitle(histNames.next())
		for i in range(n):
			binWidth = bng[i+1] - bng[i]
			hist.SetBinContent( i+1, abs(amplCorr[binCount])**2/binWidth)
			hist2.SetBinContent(i+1, abs(fCpls[w]*functions[w][binCount])**2/binWidth)
			re[i]    = amplCorr[binCount].real/binWidth**.5
			im[i]    = amplCorr[binCount].imag/binWidth**.5
			reerr[i] = coma[2*binCount  , 2*binCount  ]**.5/binWidth**.5
			imerr[i] = coma[2*binCount+1, 2*binCount+1]**.5/binWidth**.5
			jac      = [2*re[i], 2*im[i]]
			var      = 0.
			for k,K in enumerate([2*binCount, 2*binCount+1]):
				for l,L in enumerate([2*binCount, 2*binCount+1]):
					var += jac[k]*jac[l]*coma[K,L]
			hist.SetBinError(i+1, var**.5)
#			hist.SetBinError(i+1,  coma[2*binCount  ,2*binCount  ]**.5 )
#			hist2.SetBinError(i+1, coma[2*binCount+1,2*binCount+1]**.5)

			hist_re.SetBinContent(i+1, amplCorr[binCount].real/binWidth**.5)
			hist_im.SetBinContent(i+1, amplCorr[binCount].imag/binWidth**.5)

			hist_re.SetBinError(i+1, coma[2*binCount  ,2*binCount  ]**.5/binWidth**.5 )
			hist_im.SetBinError(i+1, coma[2*binCount+1,2*binCount+1]**.5/binWidth**.5 )

			hist2_re.SetBinContent(i+1, (fCpls[w]*functions[w][binCount]).real/binWidth**.5)
			hist2_im.SetBinContent(i+1, (fCpls[w]*functions[w][binCount]).imag/binWidth**.5)


			binCount += 1
		gr = ROOT.TGraphErrors(n,re,im, reerr, imerr)
		hist2.SetLineColor(2)
		hist2_re.SetLineColor(2)
		hist2_im.SetLineColor(2)

#		hist_re.Draw()
#		hist_im.Draw("SAME")
#		hist2_re.Draw("SAME")
#		hist2_im.Draw("SAME")

		gr.Draw()
#		hist.Draw()
#		hist2.Draw("SAME")
		c1.Update()
		raw_input()
		c1.Clear()

	print binCount, len(amplCorr)

if __name__ == "__main__":
	main()
