#include"efficiencyFunction.h"
#include"BELLE_efficiency.h"
#include<iostream>
#include<limits>
#include"constants.h"
#include"math.h"
#include"utils.h"
efficiencyFunction::efficiencyFunction() : _kinSignature(std::make_shared<kinematicSignature>(0)) {}

double efficiencyFunction::eval(const realVector& kin) const {
	(void) kin;
	utils::makeError("efficiencyFunction::operator()","Called base class method. Returning zero.");
	return 0.;
}

threeParticlPerfectEfficiency::threeParticlPerfectEfficiency (std::shared_ptr<kinematicSignature> kinSig) {
	_kinSignature = kinSig;
}

dalitzCut::dalitzCut(const realMatrix& cutLimits) : _cutValues(cutLimits) {
	for (const realVector& lim : _cutValues) {
		if (lim.size() != 4) {
			utils::makeError("dalitzCut(...)","Limit size is not 4: " + std::to_string(lim.size()) + ".");
			throw;
		}
	}
}

double dalitzCut::eval(const realVector& kin) const {
	double retVal = 1.;
	for (const realVector& lim : _cutValues) {
		if (lim[0] < kin[1] && lim[1] > kin[1]) {
			if (lim[2] < kin[2] && lim[3] > kin[2]) {
				retVal = 0.;
			}
		}
	}
	return retVal;
}

double threeParticlPerfectEfficiency::eval(const realVector& kin) const {
	(void) kin;
	return 1.;
}

BELLE_DtoKpipi_efficiency::BELLE_DtoKpipi_efficiency() {
	_kinSignature = std::make_shared<kinematicSignature>(2);
}

double BELLE_DtoKpipi_efficiency::eval(const realVector& kin) const {
	return Efficiency(kin);
}

BELLE_DtoKpipi_efficiency_CP::BELLE_DtoKpipi_efficiency_CP(const realVector& fs_masses) {
	_kinSignature = std::make_shared<kinematicSignature>(2);
	if (fs_masses.size() != 3) {
		std::string errMsg = "Number of final state masses needs to be three. Given:";
		for (const double& mass : fs_masses) {
			errMsg += " " + std::to_string(mass);
		}
		utils::makeError("BELLE_DtoKpipi_efficiency_CP::BELLE_DtoKpipi_efficiency_CP(...)",errMsg);
		throw;
	}
	_fs_masses_square_sum = 0.; 
	for (const double& mass : fs_masses) {
		_fs_masses_square_sum += mass*mass;
	}
}

double BELLE_DtoKpipi_efficiency_CP::eval(const realVector& kin) const {
	realVector CP_kin(3,0.);
	CP_kin[0] = kin[0];
	CP_kin[1] = kin[0] + _fs_masses_square_sum - kin[1] - kin[2];
	CP_kin[2] = kin[2];

	return Efficiency(CP_kin);
}
