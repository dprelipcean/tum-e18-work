#ifndef GETBELLEDATA
#define GETBELLEDATA
#include<vector>
#include<string>
#include"TH2D.h"
#include"types.h"
void BELLE_apply_selection_final_Dstar0ToD0piplus_D0ToKspipi_massdiffsidebandlarge(std::string infilename, std::string outfilename);
void BELLE_apply_selection_final_Dstar0ToD0piplus_D0ToKspipi_massdiffsidebandlarge_andElectronCut(std::string infilename, std::string outfilename);
void BELLE_apply_selection_final_Dstar0ToD0piplus_D0ToKspipi(const std::string inFileName, const std::string outFileName, bool massDiffSideBand = false, bool mDsideBand = false, int mDsbRegion = 0);

void writeTextDataFile(const std::string outFileName, const realMatrix dataPoints);
realMatrix readTextDataFile(const std::string inFileName, size_t variablesPerEvent);

realMatrix getBELLEevents(std::string inFileName, int SP_sign, bool SIGNSWITCH = true);

void makeDmassBinnedDalitzs(const std::string inFileName, const std::string outFileName, std::vector<std::pair<double, double> > binning, int SP_sign, bool SIGNSWITCH = true, int nBinsX = 100, int nBinsY = 100);

TH2D makeDalitzPlot(std::string name, const realMatrix data, double sMin = 0., double sMax = 4., double nBin = 200);

void makeIDplot(const std::string inFileName, const std::string outFileName);
#endif//GETBELLEDATA
