#include "getBELLEdata.h"
#include "constants.h"
#include "TTree.h"
#include "TFile.h"
#include "TCut.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <iostream>
#include <fstream>

#include"constants.h"
#include"utils.h"
#include"types.h"

double breakupMomentumSquared(const double M, const double m1, const double m2, const bool allowSubThr = false) {
	const double mSum  = m1 + m2;
	if (not allowSubThr and (M < mSum)) {
		utils::makeError("getBELLEdata::breakupMomentumSquared","Mother mass " + std::to_string(M) + " GeV/c^2 is smaller than sum of daughter masses " + std::to_string(m1) + " + " + std::to_string(m2) +  " GeV/c^2. this should never happen. Aborting...");
		throw;
	}
	const double mDiff = m1 - m2;
	double q2 = (M - mSum) * (M + mSum) * (M - mDiff) * (M + mDiff) / (4 * M * M);
	// check for rounding errors
	if (not allowSubThr and q2 < 0) {
		q2 = 0.;
	}
	return q2;
}

void BELLE_apply_selection_final_Dstar0ToD0piplus_D0ToKspipi(const std::string inFileName, const std::string outFileName, bool massDiffSideBand, bool mDsideBand, int mDsbRegion) {
	TFile *file0 = new TFile(inFileName.c_str(),"READ");
	TTree *tree0 = (TTree*)file0->Get("h1");

	TFile *filenew = new TFile(outFileName.c_str(),"RECREATE");

	std::string cut_string_selection = "p_bestcand_random==1 ";
	if (massDiffSideBand) {
		cut_string_selection += "&& 0.1504<Dstar_massdiff && Dstar_massdiff<0.16 ";
	} else { 
		cut_string_selection += "&& 0.1444<Dstar_massdiff && Dstar_massdiff<0.1464 ";
	}
	if (mDsideBand) {
		if (mDsbRegion == 0) {
			cut_string_selection +=  "&& ( ( D0_mass > 1.815 && D0_mass < 1.835) || ( D0_mass > 1.895 && D0_mass < 1.915) ) ";
		} else if (mDsbRegion == -1) {
			cut_string_selection +=  "&& ( D0_mass > 1.815 && D0_mass < 1.835) ";
		} else if (mDsbRegion ==  1) {
			cut_string_selection +=  "&& ( D0_mass > 1.895 && D0_mass < 1.915) ";
		} else {
			utils::makeError("getBELLEdata::BELLE_apply_selection_final_Dstar0ToD0piplus_D0ToKspipi(...)","Unknown mD sideband region.");
		}
	} else {
		cut_string_selection +=  "&& abs(D0_mass-1.865)<0.015 ";
	}
	cut_string_selection += "&& ( ( (e_expno==43 || e_expno==53 || e_expno==67 || e_expno==67 || e_expno==69 || e_expno==71) && Dstar_cm_ptot>3.1 ) || (Dstar_cm_ptot>2.5) )";

        TCut cut( cut_string_selection.c_str() );

	TTree *treenew = tree0->CopyTree( cut );
	
	treenew->Write();
	filenew->Close();
	file0->Close();
}

void BELLE_apply_selection_final_Dstar0ToD0piplus_D0ToKspipi_massdiffsidebandlarge(std::string infilename, std::string outfilename) {
	TFile *file0 = new TFile(infilename.c_str(),"READ");
	TTree *tree0 = (TTree*)file0->Get("h1");

	TFile *filenew = new TFile(outfilename.c_str(),"RECREATE");

    	std::string cut_string_selection = "p_bestcand_random==1 && 0.1504<Dstar_massdiff && Dstar_massdiff<0.16 && abs(D0_mass-1.865)<0.015 && ( ( (e_expno==43 || e_expno==53 || e_expno==67 || e_expno==67 || e_expno==69 || e_expno==71) && Dstar_cm_ptot>3.1 ) || (Dstar_cm_ptot>2.5) )";

        TCut cut( cut_string_selection.c_str() );

	TTree *treenew = tree0->CopyTree( cut );
	
	treenew->Write();
	filenew->Close();
	file0->Close();
}

void BELLE_apply_selection_final_Dstar0ToD0piplus_D0ToKspipi_massdiffsidebandlarge_andElectronCut(std::string infilename, std::string outfilename) {
	TFile *file0 = new TFile(infilename.c_str(),"READ");
	TTree *tree0 = (TTree*)file0->Get("h1");

	TFile *filenew = new TFile(outfilename.c_str(),"RECREATE");

    	std::string cut_string_selection = "p_bestcand_random==1 && 0.1504<Dstar_massdiff && Dstar_massdiff<0.16 && abs(D0_mass-1.865)<0.015 && ( ( (e_expno==43 || e_expno==53 || e_expno==67 || e_expno==67 || e_expno==69 || e_expno==71) && Dstar_cm_ptot>3.1 ) || (Dstar_cm_ptot>2.5) ) && D0_c1_electron_id > 0.5 && D0_c2_electron_id > 0.5";

        TCut cut( cut_string_selection.c_str() );

	TTree *treenew = tree0->CopyTree( cut );
	
	treenew->Write();
	filenew->Close();
	file0->Close();
}

void BELLE_make_sidebands(std::string infilename, std::string outfilename) {
	TFile *file0 = new TFile(infilename.c_str(),"READ");
	TTree *tree0 = (TTree*)file0->Get("h1");

	TFile *filenew = new TFile(outfilename.c_str(),"RECREATE");

    	std::string cut_string_selection = "p_bestcand_random==1 && 0.1504<Dstar_massdiff && Dstar_massdiff<0.16 && abs(D0_mass-1.865)<0.015 && ( ( (e_expno==43 || e_expno==53 || e_expno==67 || e_expno==67 || e_expno==69 || e_expno==71) && Dstar_cm_ptot>3.1 ) || (Dstar_cm_ptot>2.5) ) && D0_c1_electron_id > 0.5 && D0_c2_electron_id > 0.5";

        TCut cut( cut_string_selection.c_str() );

	TTree *treenew = tree0->CopyTree( cut );
	
	treenew->Write();
	filenew->Close();
	file0->Close();
}

void writeTextDataFile(const std::string outFileName, const realMatrix dataPoints) {
	std::ofstream outFile;
	outFile.open (outFileName.c_str(), std::ios::out);
	for (const realVector& line : dataPoints) {
		for (const double& d : line) {
			outFile << d << " ";
		}
		outFile << std::endl;
	}
	outFile.close();
}

realMatrix readTextDataFile(const std::string inFileName, size_t variablesPerEvent) {
	realMatrix retVal;
	std::ifstream inFile;
	inFile.open (inFileName.c_str(), std::ios::in);
	size_t position = 0;
	double val = 0.;
	realVector line(variablesPerEvent, 0.);
	while (inFile >> val) {
		line[position] = val;
		++position;
		if (position == variablesPerEvent) {
			position = 0;
			retVal.push_back(line);
			line = realVector (variablesPerEvent, 0.);
		}
	}
	if (position != 0) {
		utils::makeError("getBELLEdata::readBinaryDataFile(...)","Could not get full number of events.");
		throw;
	}
	inFile.close();
	return retVal;
}

realMatrix getBELLEevents(std::string inFileName, int SP_sign, bool SIGNSWITCH) {
	if (SP_sign*SP_sign*SP_sign != SP_sign) {
		utils::makeError("getBELLEevents(...)","SP_sign (soft pion) can only be 1, -1, or 0 (for positive, negative, or both).");
		throw;
	}
	TFile *file0   = new TFile(inFileName.c_str(),"READ");
	TTree *treenew = (TTree*)file0->Get("h1");

	size_t nPoints = treenew->GetEntries();
	realMatrix retVal(nPoints, realVector(3,mD0*mD0));

	TBranch* branch = treenew->GetBranch("h1");

	TLeaf* L_m01 = branch-> GetLeaf("D0_massconstrained_mass_squared01");
	TLeaf* L_m02 = branch-> GetLeaf("D0_massconstrained_mass_squared02");
	TLeaf* L_m12 = branch-> GetLeaf("D0_massconstrained_mass_squared12");
	TLeaf* L_spc = branch-> GetLeaf("softpion_charge");

	size_t count = 0;
	for (int i = 0, N = treenew->GetEntries(); i < N; ++i) {
		treenew->GetEntry(i);
		if (L_m01->GetValue(0) == 0. || L_m02->GetValue(0) == 0.|| L_m12->GetValue(0) == 0.) {
			continue;
		}

		if ((L_spc->GetValue(0) - SP_sign)*SP_sign != 0.) {
			continue;
		}
		if (L_spc->GetValue(0) == 1. or not SIGNSWITCH) {
			retVal[count][1] = L_m01->GetValue(0);
			retVal[count][2] = L_m12->GetValue(0);
		} else {
			retVal[count][1] = L_m02->GetValue(0);
			retVal[count][2] = L_m12->GetValue(0);
		}
		++count;	
	}
	utils::makeInfo("getBELLEdata::getBELLEevents(...)",std::to_string(count) + " events have been found for SP_charge of "+ std::to_string(SP_sign) + ".");
	file0->Close();
	retVal.resize(count);
	return retVal;
}

void makeDmassBinnedDalitzs(const std::string inFileName, const std::string outFileName, std::vector<std::pair<double, double> > binning, int SP_sign, bool SIGNSWITCH, int nBinsX, int nBinsY) {
	if (SP_sign*SP_sign*SP_sign != SP_sign) {
		utils::makeError("makeDmassBinedDalitzs(...)","SP_sign (soft pion) can only be 1, -1, or 0 (for positive, negative, or both).");
		throw;
	}

	TFile *file0   = new TFile(inFileName.c_str(),"READ");
	TTree *treenew = (TTree*)file0->Get("h1");

	size_t nPoints = treenew->GetEntries();
	realMatrix retVal(nPoints, realVector(3,mD0*mD0));

	TBranch* branch = treenew->GetBranch("h1");

	TLeaf* L_mD0 = branch-> GetLeaf("D0_mass");
	TLeaf* L_m01 = branch-> GetLeaf("D0_massconstrained_mass_squared01");
	TLeaf* L_m02 = branch-> GetLeaf("D0_massconstrained_mass_squared02");
	TLeaf* L_m12 = branch-> GetLeaf("D0_massconstrained_mass_squared12");
	TLeaf* L_spc = branch-> GetLeaf("softpion_charge");
	std::vector<TH2D> hists;
	std::vector<TH2D> hists_ortho;
	for (std::pair<double, double> bin: binning) {
		if (bin.first > bin.second) {
			utils::makeError("makeDmassBinedDalitzs(...)","Bin borders have to be ordered.");
			throw;
		}
		std::string name = std::string("dalitz_")+std::to_string(bin.first) + std::string("-") + std::to_string(bin.second);
		hists.push_back(TH2D(name.c_str(), name.c_str(), nBinsX, 0.,3., nBinsY ,0. , 3.));
		name += "_ortho";
		hists_ortho.push_back(TH2D(name.c_str(), name.c_str(), nBinsX, 0.,3., nBinsY ,0. , 3.));
	}
	for (int i = 0, N = treenew->GetEntries(); i < N; ++i) {
		treenew->GetEntry(i);
		if (L_m01->GetValue(0) == 0. || L_m02->GetValue(0) == 0.|| L_m12->GetValue(0) == 0.) {
			continue;
		}

		if ((L_spc->GetValue(0) - SP_sign)*SP_sign != 0.) {
			continue;
		}
		double mX;
		double mY;
		double mZ;
		double mD = L_mD0->GetValue(0);
		if (L_spc->GetValue(0) == 1. or not SIGNSWITCH) {
			mX = L_m01->GetValue(0);
			mY = L_m12->GetValue(0);
			mZ = L_m02->GetValue(0);
		} else {
			mX = L_m02->GetValue(0);
			mY = L_m12->GetValue(0);
			mZ = L_m01->GetValue(0);
		}
		for (size_t b = 0; b < binning.size(); ++b) {
			if (mD < binning[b].first || mD > binning[b].second) {
				continue;
			}
			hists[b].Fill(mX, mY);
			hists_ortho[b].Fill(mX, mZ);
		}
	}
	TFile *file = new TFile(outFileName.c_str(), "RECREATE");
	for (TH2D& h: hists) {
		h.Write();
	}
	for (TH2D& h: hists_ortho) {
		h.Write();
	}
	file->Close();
}

void makeIDplot(const std::string inFileName, const std::string outFileName) {
	TFile *file0   = new TFile(inFileName.c_str(),"READ");
	TTree *treenew = (TTree*)file0->Get("h1");
	TBranch* branch = treenew->GetBranch("h1");

	TLeaf* L_mD0 = branch-> GetLeaf("D0_mass");

	TLeaf* L_xK = branch-> GetLeaf("D0_child0_p4_px");
	TLeaf* L_yK = branch-> GetLeaf("D0_child0_p4_py");
	TLeaf* L_zK = branch-> GetLeaf("D0_child0_p4_pz");
	
	TLeaf* L_x1 = branch-> GetLeaf("D0_child1_p4_px");
	TLeaf* L_y1 = branch-> GetLeaf("D0_child1_p4_py");
	TLeaf* L_z1 = branch-> GetLeaf("D0_child1_p4_pz");

	TLeaf* L_x2 = branch-> GetLeaf("D0_child2_p4_px");
	TLeaf* L_y2 = branch-> GetLeaf("D0_child2_p4_py");
	TLeaf* L_z2 = branch-> GetLeaf("D0_child2_p4_pz");

	TLeaf* L_ep1 = branch-> GetLeaf("D0_c1_electron_id");
	TLeaf* L_ep2 = branch-> GetLeaf("D0_c2_electron_id");
	
	double mDmin = 1./0.;
	double mDmax = 0.;

	double m2eeMin = 1./0.;
	double m2eeMax = 0.;

	TH2D hist  = TH2D("mD0_vs_m2ee"   , "mD0 vs. m2(ee)"   , 100, 1.85, 1.88, 100, 0. , 2. );
	TH2D hist2 = TH2D("m2pp_vs_m2ee"  , "m2(pp) vs. m2(ee)", 100, 0.  , 2.1 , 100, 0. , 2. );
	TH2D hist3 = TH2D("cosTe_vs_m2ee" , "cosTe vs. m2(ee)" , 100,-1.  , 1.  , 100, 0. ,pow( 2.,.5) );
	TH2D hist4 = TH2D("mD0_vs_cosTe"  , "mD0 vs. cosTe"    , 100, 1.85, 1.88, 100,-1. , 1. );
	TH2D hist5 = TH2D("cosTp_vs_m2pp" , "cosT vs. m(pp)"  , 100,-1.  , 1.,    100, .26 , pow(2.1,.5) );

	TH2D hist8 = TH2D("wtf" ,"wtf"  , 100,-1.  , 1.,    100,0. , 2.1 );
	TH2D hist6 = TH2D("cosTe_vs_cosTp", "cosTe vs. cosTp"  , 100,-1.  , 1.,    100, -1. , 1. );
	TH2D hist7 = TH2D("m2ee_pee"      , "m2(ee) vs. pee"   ,   5, 0.  , 2.,    5,0. , 1. );

	TH2D ful = TH2D("full_dalitz", "fullDalitz", 100, 0.,2.2, 100 ,0.3 , 3.3);
	TH2D odd = TH2D("odd_dalitz" , "oddDalitz" , 100, 0.,2.2, 100 ,0.3 , 3.3);

	TH3D hist3d = TH3D("m2pipi_m2Kpi_cosT","m2pipi_m2Kpi_cosT",33, 0.,2.1,  33, 0., 2.1, 33, -1.,1.);
	for (int i = 0, N = treenew->GetEntries(); i < N; ++i) {
		treenew->GetEntry(i);

		double mD = L_mD0->GetValue(0);

		double pxK = L_xK->GetValue(0);
		double pyK = L_yK->GetValue(0);
		double pzK = L_zK->GetValue(0);

		TLorentzVector pK;
		pK.SetXYZM(pxK,pyK,pzK,mKs);

		double px1 = L_x1->GetValue(0);
		double py1 = L_y1->GetValue(0);
		double pz1 = L_z1->GetValue(0);

		double pe1 = L_ep1->GetValue(0);
		double pe2 = L_ep2->GetValue(0);

		double pee = pe1 * pe2;

		TLorentzVector p1e;
		TLorentzVector p1p;
		p1e.SetXYZM(px1,py1,pz1, mel);
		p1p.SetXYZM(px1,py1,pz1, mPi);

		double px2 = L_x2->GetValue(0);
		double py2 = L_y2->GetValue(0);
		double pz2 = L_z2->GetValue(0);

		TLorentzVector p2e;
		TLorentzVector p2p;
		p2e.SetXYZM(px2,py2,pz2, mel);
		p2p.SetXYZM(px2,py2,pz2, mPi);

		TLorentzVector P  = p1p + p2p + pK;
		TVector3 restBoost = P.BoostVector();

		TLorentzVector Pe = p1e + p2e + pK;
		TVector3 restBoostE = Pe.BoostVector();

		p1p.Boost(-restBoost);
		p2p.Boost(-restBoost);
		pK.Boost(-restBoost);

		p1e.Boost(-restBoostE);
		p2e.Boost(-restBoostE);

		TLorentzVector p12e = p1e + p2e;
		TLorentzVector p12p = p1p + p2p;
		TVector3 boostE = p12e.BoostVector();
		TVector3 boostP = p12p.BoostVector();
		

		double m2ee = p12e.M2();
		double m2pp = p12p.M2();
		TLorentzVector pKp = pK+p1p;
		double m2Kp = pKp.M2();

		double mpp = pow(m2pp,.5);

		double Q2 = breakupMomentumSquared(mD, mKs, mpp, true);
		if (Q2 < 0.) {
			continue;
		}
		ful.Fill(m2pp, m2Kp);
		double Q = pow(Q2,.5);
		double q = pow(breakupMomentumSquared(mpp, mPi, mPi, true),.5);

		double EK = pow(m2pp + Q*Q, .5)*pow(mKs*mKs + Q*Q, .5)/mpp + Q*Q/mpp;



		p1e.Boost(-boostE);
		p1p.Boost(-boostP);
		pK.Boost(-boostP);
		double pKp1 = (m2Kp - mPi*mPi - mKs*mKs)/2;

		double cosTformula = (pKp1 - pow(mPi*mPi + q*q,.5) *(Q*Q/mpp + pow(m2pp+Q*Q,.5)*pow(mKs*mKs + Q*Q,.5)/mpp))*mpp/q/Q/(pow(m2pp + Q*Q,.5) + pow(mKs*mKs + Q*Q,.5));


		double cosTformula2 =  (pK.E() * p1p.E() - pKp1)/pK.P()/p1p.P();

		TVector3 p12UnitE    = p12e.Vect().Unit();
		TVector3 p1restUnitE = p1e.Vect().Unit();	
		double cosTe = p1restUnitE.Dot(p12UnitE);

		TVector3 p12UnitP    = p12p.Vect().Unit();
		TVector3 p1restUnitP = p1p.Vect().Unit();
		TVector3 pKrestUnit  = pK.Vect().Unit();
		double cosTp = p1restUnitP.Dot(p12UnitP);
		double cosTp_fromK = p1restUnitP.Dot(pKrestUnit);

		mDmin = std::min(mD, mDmin);
		mDmax = std::max(mD, mDmax);

		m2eeMin = std::min(m2eeMin, m2ee);
		m2eeMax = std::max(m2eeMax, m2ee);
		hist.Fill(mD, m2ee);
		hist2.Fill(m2pp, m2ee);
		hist3.Fill(cosTe, pow(m2ee,.5));
		hist5.Fill(cosTp, pow(m2pp,.5));
		hist4.Fill(mD, cosTe);
		hist6.Fill(cosTe, cosTp);
		hist7.Fill(m2ee, pee);
		hist8.Fill(cosTformula2, m2pp);

		if (cosTp*cosTp > .64 && m2pp > 1. && m2pp < 1.5625) {
			odd.Fill(m2pp, m2Kp);
		}

		hist3d.Fill(m2pp, m2Kp, cosTp);

		double unused = EK;
		unused *= unused + cosTp_fromK + cosTformula;

	}
//	std::cout << mDmin << " | " << mDmax << " || || || " << m2eeMin << " | " << m2eeMax << std::endl;

	for (int y =0;y<hist7.GetNbinsY();++y) {
		double val = 0.;
		for (int x = 0; x < hist7.GetNbinsX(); ++x) {
			val += hist7.GetBinContent(x+1, y+1);
		}
		if (val == 0.) {
			continue;
		}
		for (int x = 0; x < hist7.GetNbinsX(); ++x) {
			hist7.SetBinContent(x+1, y+1, hist7.GetBinContent(x+1, y+1)/val);
		}
	}

	TFile* file = new TFile(outFileName.c_str(), "RECREATE");
	hist.Write();
	hist2.Write();
	hist3.Write();
	hist4.Write();
	hist5.Write();
	hist6.Write();
	hist7.Write();
	hist8.Write();

	ful.Write();
	odd.Write();

	hist3d.Write();
	file->Close();

//	std::cout << mDmin << " " << mDmax << " " << " ::: " << m2eeMin << " " << m2eeMax << std::endl;

}

TH2D makeDalitzPlot(std::string name, const realMatrix data, double sMin, double sMax, double nBin) {
	TH2D hist(name.c_str(), name.c_str(), nBin, sMin, sMax, nBin, sMin, sMax);
	for (size_t e = 0; e < data.size(); ++e) {
		hist.Fill(data[e][1], data[e][2]);
	}
	return hist;
}


/*
int main() {

}
	std::string outFileName = "./BELLE_data.root";
//	BELLE_apply_selection_final_Dstar0ToD0piplus_D0ToKspipi_massdiffsidebandlarge(inFileName, outFileName);
	TFile* outFile = new TFile("./histograms.root", "RECREATE");
	TH2D hist = makeDalitzPlot("both_true" ,getEvents(outFileName,  0, true ));
	outFile->cd();
	hist.Write();
	hist = makeDalitzPlot("pos_true"  ,getEvents(outFileName,  1, true ));
	outFile->cd();
	hist.Write();
	hist = makeDalitzPlot("neg_true"  ,getEvents(outFileName, -1, true ));
	outFile->cd();
	hist.Write();
	hist = makeDalitzPlot("both_false",getEvents(outFileName,  0, false));
	outFile->cd();
	hist.Write();
	hist = makeDalitzPlot("pos_false" ,getEvents(outFileName,  1, false));
	outFile->cd();
	hist.Write();
	hist = makeDalitzPlot("neg_false" ,getEvents(outFileName, -1, false));
	outFile->cd();
	hist.Write();
	outFile->Close();
	return 0;
}
*/
