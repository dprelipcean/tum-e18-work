#include<iostream>
#include<memory>
#include<vector>
#include<chrono>
#include<ctime>
#include<fstream>
#include<iomanip>
#include<limits>
#include <sstream>

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
#include"ROOT_fit_wrapper.h"
#include"branchFileEnding.h"
#include"masterDirectory.h"

#include"TH1D.h"
#include"TH2D.h"
#include"TFile.h"

template<typename T>
TH2D makeDalitzPlot(const std::shared_ptr<T> ampl, const double mMother, const realVector& fsMasses, const size_t nBins = 250, const double delta = 0.1) {
	// Function to produce a two-dimensional histogram
	const double xMin = pow(fsMasses[0] + fsMasses[1],2) - delta;
	const double xMax = pow(mMother     - fsMasses[1],2) + delta;
	const double yMin = pow(fsMasses[0] + fsMasses[2],2) - delta;
	const double yMax = pow(mMother     - fsMasses[1],2) + delta;

	threeParticleMassGenerator gen(mMother, fsMasses,ampl->kinSignature());

	TH2D hist(ampl->name().c_str(), ampl->name().c_str(), nBins, xMin, xMax, nBins, yMin, yMax);

	realVector kin = {mMother * mMother, 0.,0.};
	for (size_t i = 0; i < nBins; ++i) {
		kin[1] = hist.GetXaxis()->GetBinCenter(i+1);
		for (size_t j = 0; j < nBins; ++j) {
			kin[2] = hist.GetYaxis()->GetBinCenter(j+1);
			if (!gen.isValidPoint(kin)) {
				continue;
			}
			complex evl = ampl->eval(kin);
			hist.SetBinContent(i+1, j+1, evl.real()*evl.real() + evl.imag() * evl.imag());
		}
	}
	return hist;
}


int main() {
	const size_t nThreads = 7; // Number of jobs for parallel processing
	omp_set_num_threads(nThreads);

	const size_t seed = size_t( time(NULL) );
	utils::makeInfo("sandbox::main(...)","Seed: " +std::to_string(seed) + ".");
	srand(seed); // Set seed for pseudo random number generator

	const bool writeAmplitudeToROOT = false; // Flags for output
	const bool writeGeneratedEvents = false;

	std::shared_ptr<kinematicSignature> kinSig = std::make_shared<kinematicSignature>(2); // Identifier for the type of process (2 = no Bose symmetrization, since three different fs particles

	// realVector is a simple typedef of std::vector<double>  (also realMatrix, complexVector, ...)
	const realVector fs_masses = {mD0, mPi, mP0}; // Final state masses

	std::shared_ptr<threeParticleMassGenerator> generator  = std::make_shared<threeParticleMassGenerator>(mB0, fs_masses, kinSig); // Generator fir phase space distributed events
	std::shared_ptr<efficiencyFunction>         efficiency = std::make_shared<threeParticlPerfectEfficiency>(kinSig); // Detector efficiency (perfec efficiency for MC studies)

	std::shared_ptr<BELLE_S> S_ang = std::make_shared<BELLE_S>(12, fs_masses); // S-wave angular dependence for the (D0,pi-)  system (Particles 1 and 2)
	std::shared_ptr<BELLE_P> P_ang = std::make_shared<BELLE_P>(13, fs_masses); // P-wave angular dependence for the (D0,pi0)  system (Particles 1 and 3)
	std::shared_ptr<BELLE_D> D_ang = std::make_shared<BELLE_D>(23, fs_masses); // D-wave angular dependence for the (pi-,pi0) system (Particles 2 and 3)

	// The values for masses and widths here are completely arbitrary
	std::shared_ptr<BELLEbreitWigner> DmStar = std::make_shared<BELLEbreitWigner>("DmStar", 10., 5., 0, mB0, fs_masses[2], fs_masses[0], fs_masses[1]);
	std::shared_ptr<BELLEbreitWigner> D0Star = std::make_shared<BELLEbreitWigner>("D0Star", 15., 6., 1, mB0, fs_masses[1], fs_masses[0], fs_masses[2]);
	std::shared_ptr<BELLEbreitWigner> f2     = std::make_shared<BELLEbreitWigner>("f21270", 20., 7., 2, mB0, fs_masses[0], fs_masses[1], fs_masses[2]);

	std::shared_ptr<threeParticleIsobaricAmplitudeNoBose> S_wave = std::make_shared<threeParticleIsobaricAmplitudeNoBose>(12, DmStar->name(), DmStar, S_ang, fs_masses); // Full amplitude
	std::shared_ptr<threeParticleIsobaricAmplitudeNoBose> P_wave = std::make_shared<threeParticleIsobaricAmplitudeNoBose>(13, D0Star->name(), D0Star, P_ang, fs_masses); // Full amplitude
	std::shared_ptr<threeParticleIsobaricAmplitudeNoBose> D_wave = std::make_shared<threeParticleIsobaricAmplitudeNoBose>(23, f2->name(),     f2    , D_ang, fs_masses); // Full amplitude

	std::vector<std::shared_ptr<amplitude> > modelAmplitudes = {S_wave, P_wave, D_wave}; // Model to generate and fit with

	const size_t      integral_points    = 6000*10*10*10*10;
	std::shared_ptr<integrator> integral = std::make_shared<integrator>(integral_points, generator, modelAmplitudes, efficiency); // Class to calculate and hold the normalization integral
	if (!integral->makeIntegrals("ps_integralFile", "ac_integralFile", "ID_integralFile")){ // Tries to read the integral files. If it fails, they are calculated
		utils::makeError("sandbox::main(...)","Could not make the integrals.");
		return 1;
	}
	realVector norms = integral->getNormalizations().second;

	complexVector productionAmplitudes = {complex(1.,0.), complex(.345,.678), complex(0.121,-1.312)}; // Some mock values for the amplitudes

	std::shared_ptr<modelAmplitude> model    = std::make_shared<modelAmplitude>(productionAmplitudes, modelAmplitudes, norms, "modelAmplitude"); // One amplitude of the whole model
	std::shared_ptr<modelGenerator> modelGen = std::make_shared<modelGenerator>(model, generator, efficiency); // Generator for events according to the model

	if (writeAmplitudeToROOT) { // Write the amplitude to a ROOT file
		TH2D dalitz = makeDalitzPlot(model, mB0, fs_masses);
		TFile* outROOT = new TFile("outROOT.root","RECREATE");
		dalitz.Write();
		outROOT->Close();
	}

	realMatrix generatedDataPoints = modelGen->generateDataPoints(100000,100000); // Generate 10000 events according to the model

	if (writeGeneratedEvents) { // Optional output for the data points (e.g. to plot)
		std::ofstream outFile;
		outFile.open("modelGenerated");
		for (const realVector& event : generatedDataPoints) {
			for (const double& v : event) {
				outFile << v << " ";
			}
			outFile << std::endl;
		}
		outFile.close();
	}

	std::shared_ptr<logLikelihood> ll = std::make_shared<logLikelihood>(modelAmplitudes, integral); // logL function

	if (!ll->loadDataPoints(generatedDataPoints, 3)) { // Load data to the logL
		utils::makeError("sandbox::main(...)","Could not load data points.");
		return 1;
	}
	ROOT_fit_wrapper wrapper(ll);

	for (size_t i = 0; i < 100; ++i) { // Make several random attempts, since the logL function is multimodal
		realVector startValues = {utils::random2(), utils::random2(), utils::random2(), utils::random2(), utils::random2(), utils::random2()}; // Random start values

// // // Run the fit
		std::pair<double, complexVector> retVal = wrapper.fit(startValues); // retVal = (bestLogL, vector of resulting parameters)
// // // // // // //

		complex normedOne = retVal.second[1]/retVal.second[0];
		complex normedTwo = retVal.second[2]/retVal.second[0];

		utils::makeInfo("sandbox::main()","Fit results: (" +std::to_string(normedOne.real()) +"," + std::to_string(normedOne.imag()) + "), (" + std::to_string(normedTwo.real()) + "," + std::to_string(normedTwo.imag()) + ") logL = " + std::to_string(retVal.first) + ".");
	}
	return 0;
}
