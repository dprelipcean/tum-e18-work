#!/usr/bin/python
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import modernplotting.mpplot       as mpmp
import modernplotting.colors       as mpco
import modernplotting.toolkit      as mptk
import modernplotting.specialPlots as mpsp
import modernplotting.root         as mpro
import ROOT
import os,sys
def makeHist(inFileName, m2Min = 0.3, m2Max = 3.2, nBin = 250):
	hist = ROOT.TH2D("hist", "hist", nBin, m2Min, m2Max, nBin, m2Min, m2Max)
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			vals = [float(c) for c in line.split()]
			hist.Fill(vals[1], vals[2])
	return hist

def main():
#	inFileName = "./build/CP_MC_data.CP_MC"
	inFileName = sys.argv[1]

#	inFileName2 = "./build/CP_MC_data_SPD.CP_MC"

	outFileName = inFileName + ".pdf"

	print "Using:",inFileName
	hist             = makeHist(inFileName)
#	hist2 = makeHist(inFileName2)

	style            = mpmp.PlotterStyle()
	style.titleLeft  = r"D$^+\to$K$^-\pi^+\pi^+$"
	style.titleRight = r"$10^6$ simulated events"
	style.setFixed2DAxisPos(0.19*4./4., 0.16, 0.655, 0.77)
	with mptk.PdfWriter(outFileName) as pdfOutput:
		plot = style.getPlot2D()
		plot.setXlabel(r"$m^2_{\text{K}^-\pi^+_1}$\,$[\text{GeV}^2]$")
		plot.setYlabel(r"$m^2_{\text{K}^-\pi^+_2}$\,$[\text{GeV}^2]$")
		mpro.plotTH2D(hist, plot, maskValue = 0.)
		pdfOutput.savefigAndClose()
#	hist.Draw("COLZ")
#	raw_input()

#	hist.Add(hist2,-1)
#	hist.Divide(hist2)
#	with mptk.PdfWriter(inFileName + "_rat.pdf") as pdfOutput:
#		plot = style.getPlot2D()
#		plot.setXlabel(r"$m^2_{\text{K}^-\pi^+_1}$\,$[\text{GeV}^2]$")
#		plot.setYlabel(r"$m^2_{\text{K}^-\pi^+_2}$\,$[\text{GeV}^2]$")
#		mpro.plotTH2D(hist, plot, maskValue = 0.)
#		pdfOutput.savefigAndClose()	

	return 0

if __name__ == "__main__":
	sys.exit(main())
