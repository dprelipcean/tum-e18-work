#ifndef DNINE__
#define DNINE__
#include<iostream>
#include<memory>
#include<vector>
#include<chrono>
#include<ctime>
#include<fstream>
#include<iomanip>
#include<limits>

#include"massShape.h"
#include"angularDependence.h"
#include"amplitude.h"
#include"integrator.h"
#include"generator.h"
#include"efficiencyFunction.h"
#include"modelAmplitude.h"
#include"modelGenerator.h"
#include"utils.h"
#include"logLikelihood.h"
#include"constants.h"
#include"getBELLEdata.h"
#include"fabiliLL.h"
#include"fabiliLL_openMP.h"
#include"fabili.h"
#include"ROOT_fit_wrapper.h"
#include"branchFileEnding.h"
#include"masterDirectory.h"

#include"TH1D.h"
#include"TH2D.h"
#include"TFile.h"

TH2D prepareHist(const std::string& name, size_t isobX = 12, size_t isobY = 13, int nBin = 300, double sMin = 0., double sMax = 3.) {
	TH2D hist = TH2D(name.c_str(), name.c_str(), nBin, sMin, sMax, nBin, sMin, sMax);
	if (isobX == isobY) {
		std::cout << "Dnine.h::prepareHist(...): ERROR: Two times same axis encountered" << std::endl;
		throw;
	}
	if (isobX == 12) {
		hist.GetXaxis()->SetTitle("m^{2}_{K#pi_{right}}");
	} else if (isobX == 13) {
		hist.GetXaxis()->SetTitle("m^{2}_{#pi^{#plus}#pi^{#minus}}");
	} else if (isobX == 23) {
		hist.GetXaxis()->SetTitle("m^{2}_{K#pi_{wrong}}");
	} else {
		std::cout << "Dnine.h::prepareHist(...): ERROR: Invalidd isobar combination for X axis: " << isobX << std::endl;
		throw;
	}
	if (isobY == 12) {
		hist.GetYaxis()->SetTitle("m^{2}_{K#pi_{right}}");
	} else if (isobY == 13) {
		hist.GetYaxis()->SetTitle("m^{2}_{#pi^{#plus}#pi^{#minus}}");
	} else if (isobY == 23) {
		hist.GetYaxis()->SetTitle("m^{2}_{K#pi_{wrong}}");
	} else {
		std::cout << "Dnine.h::prepareHist(...): ERROR: Invalid isobar combination for Y axis: " << isobY << std::endl;
		throw;
	}
	return hist;
}

TH2D makeDataDalitz(const realMatrix& data, const realVector& fsMasses,
	size_t isobX = 12, size_t isobY = 13, int nBin = 300, double sMin = 0., double sMax = 3.) {
	if (fsMasses.size() != 3) {
		std::cout << "Dnine.h::makeDataDalitz(...): ERROR: Not exectrly three final state masses encountered" << std::endl;
		throw;
	}
	std::string name = "data";
	if (isobX == 12) {
		name = name + "_KpiRight";
	} else if (isobX == 13) {
		name = name + "_piPi";
	} else if (isobX == 23) {
		name = name + "_KpiWrong";
	}
	if (isobY == 12) {
		name = name + "_KpiRight";
	} else if (isobY == 13) {
		name = name + "_piPi";
	} else if (isobY == 23) {
		name = name + "_KpiWrong";
	}
	TH2D hist = prepareHist(name, isobX, isobY, nBin, sMin, sMax);
	double sMasses = 0.;
	for (double m: fsMasses) {
		sMasses += m*m;
	}
	for (const realVector& kin : data) {
		const double s3 = kin[0] + sMasses - kin[1] - kin[2];
		if (isobX == 12 && isobY == 13) {
			hist.Fill(kin[1], kin[2]);
		} else if (isobX == 13 && isobY == 12) {
			hist.Fill(kin[2], kin[1]);
		} else if (isobX == 12 && isobY == 23) {
			hist.Fill(kin[1], s3);
		} else if (isobX == 23 && isobY == 12) {
			hist.Fill(s3, kin[1]);
		} else if (isobX == 13 && isobY == 23) {
			hist.Fill(kin[2], s3);
		} else if (isobX == 23 && isobY == 13) {
			hist.Fill(s3, kin[2]);
		} else {
			std::cout << "makeDataDalitz(...): ERROR: Invalid combination of axes: " << isobX << " & " << isobY << std::endl;
			throw;
		}
	} 
	return hist;
}

TH2D makeDalizFromModel(const complexVector& prodAmps, 
                        const std::vector<std::shared_ptr<amplitude> >& ampls, 
                        const realVector& norms,
                        const realVector& fsMasses, 
                        std::shared_ptr<efficiencyFunction> efficiency,
                        size_t isobX = 12, size_t isobY = 13, const std::string& tag = "", int nBin = 300, double sMin = 0., double sMax = 3.) {
	if (fsMasses.size() != 3) {
		std::cout << "Dnine.h::makeDalizFromModel(...): ERROR: Not exectrly three final state masses encountered" << std::endl;
		throw;
	}
	const double binWidth = (sMax-sMin)/nBin;
	const size_t nAmpl = prodAmps.size();
	if (ampls.size() != nAmpl) {
		std::cout << "Dnine.h::makeDalizFromModel(...): ERROR: Mismatch in number of amplitudes" << std::endl;
		throw;
	}
//	if (norms.size() != nAmpl) {
//		std::cout << "Dnine.h::makeDalizFromModel(...): ERROR: Mismatch in number of normalizations" << std::endl;
//		throw;
//	}
	std::string name = "dalitz_for_"+std::to_string(nAmpl) + "_waves";
	if (nAmpl == 1) {
		name = ampls[0]->name();
	}
	name = name + tag;
	if (isobX == 12) {
		name = name + "_KpiRight";
	} else if (isobX == 13) {
		name = name + "_piPi";
	} else if (isobX == 23) {
		name = name + "_KpiWrong";
	}
	if (isobY == 12) {
		name = name + "_KpiRight";
	} else if (isobY == 13) {
		name = name + "_piPi";
	} else if (isobY == 23) {
		name = name + "_KpiWrong";
	}
	TH2D hist = prepareHist(name, isobX, isobY, nBin, sMin, sMax);
	double sMasses = mD0*mD0;
	for (double m: fsMasses) {
		sMasses += m*m;
	}
	for (int i = 0; i < nBin; ++i) {
		double sX = hist.GetXaxis()->GetBinCenter(i+1);
		for (int j = 0; j < nBin; ++j) {
			double sY = hist.GetYaxis()->GetBinCenter(j+1);
			double s3 = sMasses - sX - sY;
			realVector kin;
			if (isobX == 12 && isobY == 13) {
				kin = {mD0*mD0, sX, sY};
			} else if (isobX == 13 && isobY == 12) {
				kin = {mD0*mD0, sY, sX};
			} else if (isobX == 12 && isobY == 23) {
				kin = {mD0*mD0, sX, s3};
			} else if (isobX == 23 && isobY == 12) {
				kin = {mD0*mD0, sY, s3};
			} else if (isobX == 13 && isobY == 23) {
				kin = {mD0*mD0, s3, sX};
			} else if (isobX == 23 && isobY == 13) {
				kin = {mD0*mD0, s3, sY};
			} else {
				std::cout << "makeDalizFromModel(...): ERROR: Invalid combination of axes: " << isobX << " & " << isobY << std::endl;
			}
			if (utils::isSanePoint(kin, fsMasses)) {
				complex ampl(0.,0.);
				for (size_t a = 0; a < nAmpl; ++a) {
					ampl += prodAmps[a] * ampls[a]->eval(kin) * norms[a];
				}
				hist.SetBinContent(i+1, j+1, std::norm(ampl)*binWidth*binWidth*efficiency->eval(kin));
			}
		}
	}
	return hist;
}

std::shared_ptr<lookupAmplitudeIntens> getLookupAmplitudeFromRootFile(const std::string& fileName, const std::string& histName, const realVector& fs_masses) {
	utils::makeInfo("Dnine.h::getLookupAmplitudeFromRootFile(...)","Loading BG from '" + fileName + "'.");
	TFile* inFile = new TFile(fileName.c_str(), "READ");
	if (!inFile) {
		std::cout << "Dnineh.h::getLookupAmplitudeFromRootFile(...): ERROR: Could not open '" << fileName << "'" << std::endl;
		throw;
	}
	TH2D*    hist = (TH2D*)inFile->Get(histName.c_str());
	if (!hist) {
		std::cout << "Dnine.h::getLookupAmplitudeFromRootFile(...): ERROR: Could not (TH2D*)Get('" << histName << "' from '" << fileName << "'" << std::endl;
		throw;
	}
	double sMinX = hist->GetXaxis()->GetBinLowEdge(1);
	double widthX = hist->GetXaxis()->GetBinWidth(1);
	
	double sMinY = hist->GetYaxis()->GetBinLowEdge(1);
	double widthY = hist->GetYaxis()->GetBinWidth(1);

	std::shared_ptr<efficiencyFunction> eff = std::make_shared<BELLE_DtoKpipi_efficiency_CP>(fs_masses);

	std::cout << "Dnineh.h::getLookupAmplitudeFromRootFile(...): INFO: Found lookup table of " << hist->GetNbinsX() << " times " << hist->GetNbinsY() << " entries" << std::endl;
	
	realMatrix intensities(hist->GetNbinsX(), realVector(hist->GetNbinsY(), 0.));
	for (int i = 0; i < hist->GetNbinsX(); ++i) {
		for (int j = 0; j < hist->GetNbinsY(); ++j) {
			intensities[i][j] = hist->GetBinContent(i+1, j+1);
		}
	}
	delete hist;
	inFile->Close();
	delete inFile;
	return std::make_shared<lookupAmplitudeIntens_efficiency>(eff, std::make_shared<kinematicSignature>(2), histName + "_efficiencyCP", sMinX, widthX, sMinY, widthY, intensities);
}

template<typename T>
bool fillHistogramFromEvallableType(TH2D* histReal, TH2D* histImag, std::shared_ptr<T> ampl, double mD, size_t samplePoints = 10, bool skipNaN = false) {
	const size_t nSample= samplePoints * samplePoints;

	realVector kin = {mD*mD, 0.,0.};

	for (int i_x = 0; i_x < histReal->GetNbinsX(); ++i_x) {
		const double lowEdgeX  = histReal->GetXaxis()->GetBinLowEdge(i_x+1);
		const double highEdgeX = histReal->GetXaxis()->GetBinUpEdge(i_x+1);
		const double stepX     = (highEdgeX-lowEdgeX)/samplePoints;
		for (int j_y = 0; j_y < histReal->GetNbinsY(); ++j_y) {
			const double lowEdgeY  = histReal->GetYaxis()->GetBinLowEdge(j_y+1);
			const double highEdgeY = histReal->GetYaxis()->GetBinUpEdge(j_y+1);
			const double stepY     = (highEdgeY - lowEdgeY)/samplePoints;
			
			complex val(0.,0.);
			for (size_t i = 0; i < samplePoints; ++i) {
				kin[1] = lowEdgeX + (i+.5) * stepX;
				for (size_t j = 0; j < samplePoints; ++j) {
					kin[2] = lowEdgeY + (j+.5) * stepY;
					complex evl = ampl->eval(kin);
					if (skipNaN && (std::isnan(evl.real()) || std::isnan(evl.imag()))) {
						evl = complex(0.,0.);
					}
					val += evl;
				}
			}

			val /= nSample;
			histReal->SetBinContent(i_x+1, j_y+1, val.real());
			histImag->SetBinContent(i_x+1, j_y+1, val.imag());
		}
	}
	return true;
}

realVector get_binning(const double m_min, const double m_max, const double bin_width, const double finer_min, const double finer_max, const int finer_factor) {
	double m = m_min;
	realVector binning = {m*m};
	while (m < m_max) {
		if (m > finer_min && m < finer_max) {
			m += bin_width/finer_factor;
		} else {
			m += bin_width;
		}
		binning.push_back(m*m);
	}
	utils::makeError("Dnine.h::getBinning(...)","Something messed up during rebase.");
	throw;
	realVector mergeList;
	if (mergeList.size() > 0) {
		std::vector<double> new_binning(1, binning[0]);
		size_t count = 0;
		for (size_t chunk : mergeList) {
			count += chunk;
			if (count >= binning.size()) {
				std::cout << "Dnine.h::get_binning(...): ERROR: Chunks exceed number of bins." << std::endl;
				throw;
			}
			new_binning.push_back(binning[count]);
		}
		if (new_binning[new_binning.size()-1] != binning[binning.size()-1]) {
			std::cout << "Dnine.h::get_binning(...): ERROR: Combined chunks do not mach number of original bins." << std::endl;
			throw;			
		}
		binning = new_binning;
	}
	return binning;
}

const std::vector<std::vector<size_t> > emptyMergeMap = std::vector<std::vector<size_t> >(9, std::vector<size_t>());

realMatrix get_all_binnings(const double mD0, const double mPi, const double mKs, const std::vector<sizeVector> mergeMap = emptyMergeMap) {
	utils::makeError("Dnine.h::get_all_binnings","Something messes up during rebase.");
	throw;
	std::cout << mergeMap.size() << std::endl;// To avoid unused parameter error.	

	const double bin_width = 0.04;
	realMatrix retVal(9);
	retVal[0] = get_binning(  mPi + mKs, mD0 - mPi, bin_width, 20.   , 20.,    1); // No finer binning here (K pi S)
	retVal[1] = get_binning(  mPi + mKs, mD0 - mPi, bin_width,   .789,   .999, 4); // K pi P
	retVal[2] = get_binning(  mPi + mKs, mD0 - mPi, bin_width,  1.329,  1.529, 2); // K pi D
	retVal[3] = get_binning(  mPi + mKs, mD0 - mPi, bin_width, 20.   , 20.,    1); // No finer binning here (K pi S)
	retVal[4] = get_binning(  mPi + mKs, mD0 - mPi, bin_width,   .789,   .999, 4); // K pi P
	retVal[5] = get_binning(  mPi + mKs, mD0 - mPi, bin_width,  1.329,  1.529, 2); // K pi D
	retVal[6] = get_binning(2*mPi      , mD0 - mKs, bin_width,   .919,  1.079, 4); // pi pi S
	retVal[7] = get_binning(2*mPi      , mD0 - mKs, bin_width,   .679,   .919, 2); // pi pi P
	retVal[8] = get_binning(2*mPi      , mD0 - mKs, bin_width,  1.159,  1.39,  2); // pi pi D
	return retVal;
}

std::vector<std::shared_ptr<massShape> > get_isobars(const size_t waveIndex, const double mD0, const double mPi, const double mKs, bool writeHists = false) {
	std::vector<std::shared_ptr<massShape> > isobars;
	if (waveIndex > 8) {
		std::cout << "get_isobars(...): ERROR: Invalid wave index (>8) " << waveIndex << std::endl;
		return isobars;
	}
	double mRho      =  .77526; double Grho      = .1478;
	double mRhoPrime = 1.465;   double GrhoPrime = .4; // This is never on-shell?? How can we treat this?
	double mf2       = 1.2751;  double Gf2       = .1851;
	double mOmega    =  .78265; double Gomega    = .00849;
//	double mK892     =  .89166; double GK892     = .0508;
	double mK892     =  .8937;  double GK892     = .0472;
//	double mK01430   = 1.425;   double GK01430   = .270;
	double mK21430   = 1.4256;  double GK21430   = .0985;
	double mK1410    = 1.414;   double GK1410    = .232;
	double mK1680    = 1.718;   double GK1680    = .322;
//	Ordering: "a", "r", "M0", "G0", "phiF", "phiR", "phiRsin", "F", "R", "MMax" // phiRsin = phiR ("Typo in some paper?? " - S. Wallner); chose MMax
//	realVector LASS_parameters    = {0.113, -33.8, 1.441, 0.193, utils::degToRad(0.1), utils::degToRad(-109.7), 0.96, 1.};
	realVector LASS_parameters    = {0.113, -33.797, 1.4405, 0.1926, utils::degToRad(0.099), utils::degToRad(-109.695), 0.955, 1.};
	realVector KMatrix_parameters = {  8.5212, utils::degToRad(  68.505), // As in the BELLE note
	                                           12.1895, utils::degToRad(  23.950), 
	                                           29.1463, utils::degToRad(  -0.106), 
	                                           10.7457, utils::degToRad( -51.893), 
	                                            0.    , utils::degToRad(   0.   ), 
	                                            8.0441, utils::degToRad(-125.963), 
	                                           26.2981, utils::degToRad(-152.324), 
	                                           33.0346, utils::degToRad( -93.228), 
	                                           26.1735, utils::degToRad(-121.405), 
	                                            0.    , utils::degToRad(   0.   ), 
	                                          -0.07};
//	realVector KMatrix_parameters = {  9.3 , utils::degToRad( -78.7), // As in the BaBar paper [https://arxiv.org/pdf/0804.2089.pdf]
//	                                           10.89, utils::degToRad(-159.1), 
//	                                           24.2 , utils::degToRad( 168.0), 
//	                                            9.16, utils::degToRad(  90.5), 
//	                                            0.  , utils::degToRad(   0. ), 
//	                                            7.94, utils::degToRad(  73.9), 
//	                                            2.0 , utils::degToRad( -18. ), 
//	                                            5.1 , utils::degToRad(  33. ), 
//	                                            3.23, utils::degToRad(   4.8), 
//	                                            0.  , utils::degToRad(   0. ), 
//	                                          -0.07};

	if (waveIndex == 0 || waveIndex == 3) { // Kpi S
		isobars.push_back(std::make_shared<BELLE_LASS>(LASS_parameters, mPi, mKs));
	} else if (waveIndex == 1 || waveIndex == 4) { // K pi P
		isobars.push_back(std::make_shared<BELLEbreitWigner>("K*(892)", mK892, GK892, 1, mD0, mPi, mKs, mPi));
		isobars.push_back(std::make_shared<BELLEbreitWigner>("K*(1410)", mK1410, GK1410, 1, mD0, mPi, mKs, mPi));
		if (waveIndex == 4) { // This is only used in one combination
//		if (waveIndex == 1) {
			isobars.push_back(std::make_shared<BELLEbreitWigner>("K*(1680)", mK1680, GK1680, 1, mD0, mPi, mKs, mPi));
		}
	} else if (waveIndex == 2 || waveIndex == 5) { // K pi D
		isobars.push_back(std::make_shared<BELLEbreitWigner>("K2*(1430)", mK21430, GK21430, 2, mD0, mPi, mKs, mPi));
	} else if (waveIndex == 6) { // pi pi S
		isobars.push_back(std::make_shared<BELLE_KMatrix>(KMatrix_parameters)); 
	} else if (waveIndex == 7) { // pi pi P
		isobars.push_back(std::make_shared<BELLEbreitWigner>("omega", mOmega, Gomega, 1, mD0, mKs, mPi, mPi));
		isobars.push_back(std::make_shared<BELLEbreitWigner>("rho(770)", mRho, Grho, 1, mD0, mKs, mPi, mPi));
		isobars.push_back(std::make_shared<BELLEbreitWigner>("rho'(1450)", mRhoPrime, GrhoPrime, 1, mD0, mKs, mPi, mPi));
	} else if (waveIndex == 8) { // pi pi D
		isobars.push_back(std::make_shared<BELLEbreitWigner>("f2(1270)", mf2, Gf2, 2, mD0, mKs, mPi, mPi));
	} else {
		std::cout << "Dnine.h::get_isobars(...): ERROR: Invalid wave index " << waveIndex << std::endl;
		return isobars;
	}
	if (writeHists) {
		TFile* outFile = new TFile("isobars.root", "UPDATE");
		double mThresh = mPi + mKs;
		if (waveIndex > 5) {
			mThresh = 2*mPi;
		}
		for (std::shared_ptr<massShape> ms : isobars) {
			const size_t nBins = 1000;
			TH1D hist(ms->name().c_str(), ms->name().c_str(), nBins, 0., mD0);
			for (size_t b = 0; b < nBins; ++b) {
				double m = hist.GetXaxis()->GetBinCenter(b+1);
				if (m < mThresh) {
					continue;
				}
				hist.SetBinContent(b+1, std::norm(ms->eval(m*m)));
			}
			hist.Write();
		}
		outFile->Close();
	}
	return isobars;
}

sizeVector getNwaves(const std::vector<bool>& freeMap) {
	if (freeMap.size() != 9) {
		std::cout << "Dnine.h::getNwaves(): ERROR: freeMap has to have .size == 9" << std::endl;
		throw;
	}
	realMatrix binnings = get_all_binnings(mD0, mPi, mKs);
	sizeVector retVal(9,0);
	for (size_t w = 0; w < 9; ++w) {
		if (freeMap[w]) {
			retVal[w] = binnings[w].size() - 1;
		} else {
			retVal[w] = get_isobars(w, mD0, mPi, mKs, false).size();
		}
	}
	return retVal;
}

std::vector<std::shared_ptr<amplitude> > get_model(std::vector<bool> free_map, const double mD0, const double mPi, const double mKs, bool cp_conj = false) {

	const realVector fs_masses = {mPi, mKs, mPi};

	std::vector<std::shared_ptr<amplitude> > model;

	const size_t piKrightIndex = cp_conj ? 23 : 12;
	const size_t piKwrongIndex = cp_conj ? 12 : 23;
	const size_t piPiIndex     = 13;

	std::shared_ptr<BELLE_S> S_right = std::make_shared<BELLE_S>(piKrightIndex, fs_masses);
	std::shared_ptr<BELLE_S> S_pipi  = std::make_shared<BELLE_S>(piPiIndex,     fs_masses);
	std::shared_ptr<BELLE_S> S_wrong = std::make_shared<BELLE_S>(piKwrongIndex, fs_masses);

	std::shared_ptr<BELLE_P> P_right = std::make_shared<BELLE_P>(piKrightIndex, fs_masses);
	std::shared_ptr<BELLE_P> P_pipi  = std::make_shared<BELLE_P>(piPiIndex,     fs_masses);
	std::shared_ptr<BELLE_P> P_wrong = std::make_shared<BELLE_P>(piKwrongIndex, fs_masses);

	std::shared_ptr<BELLE_D> D_right = std::make_shared<BELLE_D>(piKrightIndex, fs_masses);
	std::shared_ptr<BELLE_D> D_pipi  = std::make_shared<BELLE_D>(piPiIndex,     fs_masses);
	std::shared_ptr<BELLE_D> D_wrong = std::make_shared<BELLE_D>(piKwrongIndex, fs_masses);

	realMatrix binnings = get_all_binnings(mD0, mPi, mKs);
	realVector binning;

	sizeVector isob_combinations = {piKrightIndex, piKrightIndex, piKrightIndex, piKwrongIndex, piKwrongIndex, piKwrongIndex, piPiIndex, piPiIndex, piPiIndex};
	std::vector<std::string> prefixes = {"KpiRight[","KpiRight[","KpiRight[","KpiWrong[","KpiWrong[","KpiWrong[", "piPi[", "piPi[", "piPi["};
	std::vector<std::string> suffixes = {"]PiS", "]PiP", "]PiD","]PiS", "]PiP", "]PiD", "]KS", "]KP", "]KD"};
	std::vector<std::shared_ptr<angularDependence> > angular_dependences = {S_right, P_right, D_right, S_wrong, P_wrong, D_wrong, S_pipi, P_pipi, D_pipi};
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	for (size_t d = 0; d < 9; ++d) {
		if (free_map[d]) {
			binning = binnings[d];
			for (size_t b = 0; b < binning.size() - 1; ++b) {
				std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binning[b], binning[b+1]);
				std::shared_ptr<threeParticleIsobaricAmplitudeNoBose> stepWave = std::make_shared<threeParticleIsobaricAmplitudeNoBose>(isob_combinations[d], prefixes[d] + std::to_string(b) + suffixes[d], step, angular_dependences[d], fs_masses);
				model.push_back(stepWave);
			}
		} else {
			std::vector<std::shared_ptr<massShape> > isobars = get_isobars(d, mD0, mPi, mKs, false);
			for (std::shared_ptr<massShape>& isobar : isobars) {
				std::shared_ptr<threeParticleIsobaricAmplitudeNoBose> wave = std::make_shared<threeParticleIsobaricAmplitudeNoBose>(isob_combinations[d], prefixes[d] + isobar->name() + suffixes[d], isobar, angular_dependences[d], fs_masses);
				model.push_back(wave);
			}
		}
	}
	return model;
}

std::shared_ptr<amplitude> get_bg_amplitude(const realVector& fs_masses) {
//	std::shared_ptr<dalitzPolynomialAmplitude> bg_amplitude = std::make_shared<dalitzPolynomialAmplitude>(std::make_shared<kinematicSignature>(2), "dalitzBackground_version2.pol", .405, 1.845, .407, 2.057);
//	return bg_amplitude;
	return getLookupAmplitudeFromRootFile(masterDirectory + "/../polotossos.root","Dalitz_model", fs_masses);
}

complexVector getBelleProdAmps(
	realVector fitFractions = {0.0001 , 0.0055 , 0.0004 , 0.0001 , 0.0699, 0.5987 , 0.0007, 0.0054  , 0.013  , .0, 0.0054 , 0.2042, 0.0059 , 0.0076}, 
	realVector phasesInDeg  = {162.313, -42.164, 150.166, -89.608, 99.417, 136.843, 99.386, -118.157, -44.066, 0., 120.743, 0.0   , 102.107, -36.285}) { // Assume pipiS phase to be zero
	size_t nAmpl = fitFractions.size();
	if (nAmpl != phasesInDeg.size()) {
		std::cout << "Dnine.h::getBelleProdAmps: ERROR: Size mismatch (" << nAmpl << " vs. " << phasesInDeg.size() << ")" << std::endl;
		throw;
	}
	complexVector retVal(nAmpl);
	for (size_t a = 0; a < nAmpl; ++a) {
		retVal[a] = std::polar(pow(fitFractions[a],.5), utils::degToRad(phasesInDeg[a]));
	}
	return retVal;
}

void randomizeStartValues(complexVector ampls, double scaleMin, double scaleMax) {
	double scaleRange = scaleMax - scaleMin;
	for (size_t a = 0; a < ampls.size(); ++a) {
		double scale = scaleMin + utils::random() * scaleRange;
		ampls[a] *= scale;
	}
}

std::vector<std::string> wave_names = {"KpiRightS", "KpiRightP", "KpiRightD", "KpiWrongS", "KpiWrongP", "KpiWrongD", "piPiS", "piPiP", "piPiD"};

complexVector getRandomizedStartValues(const sizeVector& n_waves, 
                                                            double                     CP_scale_factor, 
                                                            bool                       anchorFirst, 
                                                            bool                       copyRightToWrong, 
                                                            complex       bg_amplitude, 
                                                            bool                       fixBG = false) {
	complexVector retVal;
	size_t nKright = n_waves[0] + n_waves[1] + n_waves[2];
	size_t n_model = 0;
	const double expectedDifference = 10.; // wrong waves is expected to be 100 + 10*10 times higher than right wave
	complexVector kScales = {expectedDifference * complex(utils::random2(), utils::random2()),
	                                              expectedDifference * complex(utils::random2(), utils::random2()),
	                                              expectedDifference * complex(utils::random2(), utils::random2())};
	for (size_t w = 0; w < 9; ++w) {
		for (size_t n = 0; n < n_waves[w]; ++n) {
			if ((w == 3 || w == 4 || w == 5) and copyRightToWrong) {
				retVal.push_back(kScales[w-3] * retVal[n_model + n - nKright]);
			} else {
				retVal.push_back(complex(utils::random2(), utils::random2()));	
			}
		}
		n_model += n_waves[w];
	}
	for (size_t n = 0; n < n_model; ++n) {
		retVal.push_back(retVal[n] * CP_scale_factor);
	}
	if (anchorFirst) {
		complex phase = retVal[0] / pow(std::norm(retVal[0]), .5);
		for (size_t n = 0; n < 2*n_model; ++n) {
			retVal[n] /= phase;
		}
	}
	if (!fixBG) {
		retVal.push_back(bg_amplitude);
	}
	return retVal;
}

realVector getParamsFromProdAmps(const complexVector& prodAmps, 
                                          const sizeVector&                n_waves, 
                                          bool                                      anchorFirst, 
                                          bool                                      copyRightToWrong, 
                                          const                                     std::vector<bool>& fixToZero, 
                                          bool                                      fixBG = false) {
	realVector retVal;
	const size_t nKright = n_waves[0] + n_waves[1] + n_waves[2];
	const size_t nKwrong = n_waves[3] + n_waves[4] + n_waves[5];
	const size_t nPiPi   = n_waves[6] + n_waves[7] + n_waves[8];
	const size_t n_model = nKright + nKwrong + nPiPi;
	if (!fixToZero[0]) {
		retVal.push_back(prodAmps[0].real());
		if (!anchorFirst) { // Kright
			retVal.push_back(prodAmps[0].imag());
		}
	}
	size_t count = 0;
	for (size_t w = 0; w < 3; ++w) {
		for (size_t n = 0; n < n_waves[w]; ++n) {
			if (!fixToZero[w]) {
				if (w == 0 && n == 0) { // Handeled above because of anchor
					continue;
				}
				retVal.push_back(prodAmps[count].real());
				retVal.push_back(prodAmps[count].imag());
			}
			++count;
		}
	}
	complexVector scaleFactors;
	if (copyRightToWrong) {
		size_t countRight = 0;
		size_t countWrong = nKright;
		for (size_t w = 0; w < 3; ++w) {
			if (fixToZero[w+3]) {
				std::cout << "Dnine.h::getParamsFromProdAmps(...): ERROR: Cannot copy and fix wave " << w+3 << " to zero" << std::endl;
				throw;
			}
			scaleFactors.push_back(prodAmps[countWrong]/prodAmps[countRight]);
			if (n_waves[w+3] > n_waves[w]) { // Fixed KpiPwrong wave has on
				for (size_t i = 0; i < n_waves[w+3] - n_waves[w]; ++i) {
					complex pa = prodAmps[nKright + countRight + i];
					retVal.push_back(pa.real());
					retVal.push_back(pa.imag());
				}
			}
			countRight += n_waves[w];
			countWrong += n_waves[w+3];
		}
	} else {
		count = nKright;
		for (size_t w = 3; w < 6; ++w) {		
			for (size_t n = 0; n < n_waves[w]; ++n) {
				if (!fixToZero[w]) {
					retVal.push_back(prodAmps[count].real());
					retVal.push_back(prodAmps[count].imag());
				}
				++count;
			}
		}
	}

	count = nKright + nKwrong;
	for (size_t w = 6; w < 9; ++w) { //PiPi
		for (size_t n = 0; n < n_waves[w]; ++n) {
			if (!fixToZero[w]) {
				retVal.push_back(prodAmps[count].real());
				retVal.push_back(prodAmps[count].imag());
			}
			++count;
		}
	}
	if (!fixBG) {
		retVal.push_back(prodAmps[2*n_model].real()); // The background
	}
	for (complex& scale : scaleFactors) { // The scale factors
		retVal.push_back(scale.real());
		retVal.push_back(scale.imag());
	}
	return retVal;
}

complexVector getStartValuesForBelleParameters(const std::vector<bool>& freeMap, const realVector& normsFree, const realVector& normsFixed) {
	if (freeMap.size() != 9) {
		std::cout << "Dnine.h::getStartValuesForBelleParameters(...): ERROR: Invalid size of freeMap" << std::endl;
		throw;
	}
	complexVector retVal;
	complexVector belleValues = getBelleProdAmps();
	realMatrix  binnings    = get_all_binnings(mD0, mPi, mKs);
	size_t isobarCount = 0;
	size_t freeCount   = 0;
	for (size_t wave = 0; wave < 9; ++wave) {
		std::vector<std::shared_ptr<massShape> > isobars = get_isobars(wave, mD0, mPi, mKs, false);
		if (freeMap[wave]) {
			complexVector coefficients(isobars.size());
			for (size_t i = 0; i < isobars.size(); ++i) {
				coefficients[i] = belleValues[isobarCount] * normsFixed[isobarCount];
				++isobarCount;
			}
			for (size_t b = 0; b < binnings[wave].size()-1; ++b) {
				double m2 = (pow(binnings[wave][b],.5) + pow(binnings[wave][b+1],.5))/2.; // bin center 
				m2 *= m2; // Get m^2 = s
				complex ampl = 0.;
				for (size_t i = 0; i < isobars.size(); ++i) {
					ampl += coefficients[i] * isobars[i]->eval(m2);
				}
				retVal.push_back(ampl/normsFree[freeCount]);
				++freeCount;
			}
		} else {
			for (size_t i = 0; i < isobars.size(); ++i){
				retVal.push_back(belleValues[isobarCount]);
				++freeCount;
				++isobarCount;
			}
		}
	}
	return retVal;
}

void checkDerivativesNlopt(std::shared_ptr<logLikelihoodBase> ll, const realVector &params, double delta, bool doSecond) {
	realVector par = params;
	realVector Devl(params.size(), 0.);
	const double evl = ll->nloptCall(par, Devl);
	realVector alibi(params.size(), 0.);
	for (size_t p = 0; p < par.size(); ++p) {
		double dd = delta * par[p];
		par[p] += dd;
		std::cout << Devl[p] << " || " << (ll->nloptCall(par, alibi) - evl)/dd << std::endl;
		par[p] -= dd;
	}
	if (doSecond) {
/*			const std::vector<std::vector<double> > DDevl = ll->DDeval(pa);
		std::cout << "utils::checkDerivatives(...): INFO: Second derivative gotten. Checking numerically now..." << std::endl;
		std::vector<double> D;
		for (size_t a = 0; a < pa.size(); ++a) {
			double dd = delta * pow(std::norm(pa[a]),.5);
			pa[a] += std::complex<double>(dd,0.);
			D = ll->Deval(pa);
			for (size_t b = 0; b < pa.size(); ++b) {
				if (a == b) {
					std::cout << " ! - ! - ! - ! - ! - ! - !" << std::endl;
				}
				std::cout << DDevl[2*a  ][2*b  ] << " ++ " << (D[2*b  ] - Devl[2*b  ])/dd << std::endl;
				std::cout << DDevl[2*a  ][2*b+1] << " ++ " << (D[2*b+1] - Devl[2*b+1])/dd << std::endl;
				if (a == b) {
					std::cout << " ! - ! - ! - ! - ! - ! - !" << std::endl;
				}
			}
			pa[a] += std::complex<double>(-dd,dd);
			D = ll->Deval(pa);
			for (size_t b = 0; b < pa.size(); ++b) {
				if (a == b) {
					std::cout << " ! - ! - ! - ! - ! - ! - !" << std::endl;
				}
				std::cout << DDevl[2*a+1][2*b  ] << " ++ " << (D[2*b  ] - Devl[2*b  ])/dd << std::endl;
				std::cout << DDevl[2*a+1][2*b+1] << " ++ " << (D[2*b+1] - Devl[2*b+1])/dd << std::endl;
				if (a == b) {
					std::cout << " ! - ! - ! - ! - ! - ! - !" << std::endl;
				}
			}
			pa[a] += std::complex<double>(0.,-dd);
		}
*/
	}
}

bool doTheFixingAndCopying(std::shared_ptr<logLikelihood> ll, 
                           const sizeVector& n_waves, 
                           double                     CP_scale_factor, 
                           bool                       copyRightToWrong, 
                           const std::vector<bool>&   fixToZeroList, 
                           bool                       anchorFirst, 
                           bool                       fixBG  = false,
                           complex       BGampl = complex(0.,0.)) {

	const size_t nKright = n_waves[0] + n_waves[1] + n_waves[2];
	const size_t nKwrong = n_waves[3] + n_waves[4] + n_waves[5];
	const size_t nPiPi   = n_waves[6] + n_waves[7] + n_waves[8];
	const size_t n_model  = nKright + nKwrong + nPiPi;

	if (copyRightToWrong) {
		if (fixToZeroList[3] || fixToZeroList[4] || fixToZeroList[5]) {
			std::cout << "Dnine::doTheFixingAndCopying(...): ERROR: Cannot fix Kwrong waves to zero, if they are copied from Kright waves." << std::endl;
			return false;
		}
		size_t countSource   = 0;
		size_t countTarget   = 0;
		for (size_t K = 0; K < 3; ++K) {
			for (size_t a = 0; a < n_waves[K]; ++a) {
// // // // // // // // // // // // // // // //
// // // // // // // // // // // // // // // //
// // // // Copy right to wrong
// // // real to real
				if (!ll->addCopyParameter(2*(nKright + countTarget + a), 2*(countSource +a), 2*K, 1.)) {
					std::cout << "Dnine::doTheFixingAndCopying(...): ERROR: Could not copy parameter " << 2*(nKright + countTarget + a) << " from " << 2*(countSource +a) << "(realRight to realWrong)" << std::endl;
					return false;
				}
// // // // // // // // // // // // // // // //
// // // imag to real
				if (!ll->addCopyParameter(2*(nKright + countTarget + a), 2*(countSource +a)+1, 2*K+1,-1.)) {
					std::cout << "Dnine::doTheFixingAndCopying(...): ERROR: Could not copy parameter " << 2*(nKright + countTarget + a) << " from " << 2*(countSource +a)+1 << "(imagRight to realWrong)" << std::endl;
					return false;
				}
// // // // // // // // // // // // // // // //
// // // imag to imag
				if (!ll->addCopyParameter(2*(nKright + countTarget + a)+1, 2*(countSource +a)+1, 2*K, 1.)) {
					std::cout << "Dnine::doTheFixingAndCopying(...): ERROR: Could not copy parameter " << 2*(nKright + countTarget + a)+1 << " from " << 2*(countSource +a)+1 << "(imagRight to imagWrong)" << std::endl;
					return false;
				}
// // // // // // // // // // // // // // // //
// // // real to imag
				if (!ll->addCopyParameter(2*(nKright + countTarget + a)+1, 2*(countSource +a), 2*K+1, 1.)) {
					std::cout << "Dnine::doTheFixingAndCopying(...): ERROR: Could not copy parameter " << 2*(nKright + countTarget + a)+1 << " from " << 2*(countSource +a) << "(realRight to imagWrong)" << std::endl;
					return false;
				}
// // // // // // // // // // // // // // // //
// // // // // // // // // // // // // // // //
// // // // Copy right to rightCP
// // // real to real
				if (!ll->addCopyParameter(2*(n_model + countSource + a), 2*(countSource + a), -1, CP_scale_factor)) {
					std::cout << "Dnine::doTheFixingAndCopying(...): ERROR: Could not copy parameter " << 2*(n_model + countSource + a) << " from " <<  2*(countSource + a) << "(realRight to realRightCP)" << std::endl;
					return false;
				}
				if (!ll->addCopyParameter(2*(n_model + countSource + a)+1, 2*(countSource + a)+1, -1, CP_scale_factor)) {
					std::cout << "Dnine::doTheFixingAndCopying(...): ERROR: Could not copy parameter " << 2*(n_model + countSource + a)+1 << " from " <<  2*(countSource + a)+1 << "(imagRight to imagRightCP)" << std::endl;
					return false;
				}
// // // // // // // // // // // // // // // //
// // // // // // // // // // // // // // // //
// // // // Copy right to wrongCP
// // // real to real
				if (!ll->addCopyParameter(2*(n_model + nKright + countTarget + a), 2*(countSource +a), 2*K, CP_scale_factor)) {
					std::cout << "Dnine::doTheFixingAndCopying(...): ERROR: Could not copy parameter " << 2*(n_model + nKright + countTarget + a) << " from " << 2*(countSource +a)<< " (realRight to realWrongCP)" << std::endl;
					return false;
				}
// // // // // // // // // // // // // // // //
// // // imag to real
				if (!ll->addCopyParameter(2*(n_model + nKright + countTarget + a), 2*(countSource +a)+1, 2*K+1,-CP_scale_factor)) {
					std::cout << "Dnine::doTheFixingAndCopying(...): ERROR: Could not copy parameter " << 2*(n_model + nKright + countTarget + a) << " from " << 2*(countSource +a)+1 << "(imagRight to realWrongCP)" << std::endl;
					return false;
				}
// // // // // // // // // // // // // // // //
// // // imag to imag
				if (!ll->addCopyParameter(2*(n_model + nKright + countTarget + a)+1, 2*(countSource +a)+1, 2*K, CP_scale_factor)) {
					std::cout << "Dnine::doTheFixingAndCopying(...): ERROR: Could not copy parameter " << 2*(n_model + nKright + countTarget + a)+1 << " from " << 2*(countSource +a)+1 << "(imagRight to imagWrongCP)" << std::endl;
					return false;
				}
// // // // // // // // // // // // // // // //
// // // real to imag
				if (!ll->addCopyParameter(2*(n_model + nKright + countTarget + a)+1, 2*(countSource +a), 2*K+1, CP_scale_factor)) {
					std::cout << "Dnine::doTheFixingAndCopying(...): ERROR: Could not copy parameter " << 2*(n_model + nKright + countTarget + a)+1 << " from " << 2*(countSource +a) << "(realRight to imagWrongCP)" << std::endl;
					return false;
				}
			}
// // // // // // // // // // // // // // // //
// // // // // // // // // // // // // // // //
// // // // Copy wrong to wrongCP
			for (size_t a = n_waves[K]; a < n_waves[K+3]; ++a) { // This covers for the K*(1680), which is only in KpiWrong (but works generally for additional (appended) waves in a sector)
// // // real to real
				if (!ll->addCopyParameter(2*(n_model + nKright + countTarget + a), 2*(nKright + countTarget + a), -1, CP_scale_factor)) {
					std::cout << "Dnine::doTheFixingAndCopying(...): ERROR: Could not copy parameter " << 2*(n_model + countSource + a) << " from " <<  2*(nKright + countTarget + a) << "(realWrong to realWrongCP)" << std::endl;
					return false;
				}
// // // imag to imag
				if (!ll->addCopyParameter(2*(n_model + nKright + countTarget + a)+1, 2*(countSource + a)+1, -1, CP_scale_factor)) {
					std::cout << "Dnine::doTheFixingAndCopying(...): ERROR: Could not copy parameter " << 2*(n_model + nKright + countTarget + a)+1 << " from " <<  2*(nKright + countTarget + a)+1 << "(imagWrong to imagWrongCP)" << std::endl;
					return false;
				}
			}
			countSource += n_waves[K];
			countTarget += n_waves[K+3];
		}	
	} else {
		for (size_t a = 0; a < n_model; ++a) {
			if (!ll->addCopyParameter(2*(n_model + a), 2*(a), -1, CP_scale_factor)) {
				std::cout << "Dnine::doTheFixingAndCopying(...): ERROR: Could not copy parameter " << 2*(n_model + a) << " from " <<  2*a << "(real to realCP)" << std::endl;
				return false;
			}
			if (!ll->addCopyParameter(2*(n_model + a)+1, 2*(a)+1, -1, CP_scale_factor)) {
				std::cout << "Dnine::doTheFixingAndCopying(...): ERROR: Could not copy parameter " << 2*(n_model + a)+1 << " from " <<  2*(a)+1 << "(imag to imagCP)" << std::endl;
				return false;
			}
		}
	}
// // // // // // // // // // // // // // // //
// // // // // // // // // // // // // // // //
// // // // Copy pipi to pipiCP
	for (size_t a = nKright + nKwrong; a < n_model; ++a) {
		if (!ll->addCopyParameter(2*(n_model+a), 2*a, -1, CP_scale_factor)) {
			std::cout << "Dnine::doTheFixingAndCopying(...): ERROR: Could not copy parameter " << 2*(n_model+a) << " from " << 2*a << "(realPiPi to realPiPiCP)" << std::endl;
			return false;
		}
		if (!ll->addCopyParameter(2*(n_model+a)+1, 2*a+1, -1, CP_scale_factor)) {
			std::cout << "Dnine::doTheFixingAndCopying(...): ERROR: Could not copy parameter " << 2*(n_model+a)+1 << " from " << 2*a+1 << "(imagPiPi to imagPiPiCP)" << std::endl;
			return false;
		}

	}
	size_t countAmpl = 0;
	for (size_t w = 0; w < 9; ++w) {
		if (fixToZeroList[w]) {
			for (size_t Ari = 2*countAmpl; Ari < 2*(countAmpl + n_waves[w]); ++Ari) {
				if (!ll->fixParameter(Ari, 0.)) {
					std::cout << "Dnine::doTheFixingAndCopying(...): ERROR: Could not fix parameter " << Ari << " to zero." << std::endl;
					return false;
				}
			}
		}
		countAmpl += n_waves[w];
	}
	if (anchorFirst) {
		if (!ll->fixParameter(1, 0.)) {
			std::cout << "Dnine::doTheFixingAndCopying(...): ERROR: Could not fix the annchor wave's real part to zero" << std::endl;
			return false;
		}
	}
	if (fixBG) {
		if (!ll->fixParameter(4*n_model  , BGampl.real())) { // Fix real part of the BG amplitude
			std::cout << "Dnine::doTheFixingAndCopying(...): ERROR: Could not fix BG amplitude's real part (" <<2*n_model << ") to " << BGampl.real() << std::endl;
			return false;
		}
		if (!ll->fixParameter(4*n_model+1, BGampl.imag())) { // Fix imag part of the BG amplitude
			std::cout << "Dnine::doTheFixingAndCopying(...): ERROR: Could not fix BG amplitude's imag part (" <<2*n_model+1 << ") to " << BGampl.imag() << std::endl;
			return false;
		}
	} else {
		if (!ll->fixParameter(4*n_model+1, 0.)) { // Fix imag part of the BG amplitude to zero
			std::cout << "Dnine::doTheFixingAndCopying(...): ERROR: Could not fix BG amplitude's imag part (" <<2*n_model+1 << ") to zero" << std::endl;
			return false;
		}
	}
	return true;
}
#endif//DNINE__
