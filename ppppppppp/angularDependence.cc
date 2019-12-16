#include"angularDependence.h"
#include"constants.h"
#include"utils.h"
#include<iostream>
angularDependence::angularDependence(std::shared_ptr<kinematicSignature> kinSignature, std::string name):
	_kinSignature(kinSignature), _name(name) {}

complex angularDependence::eval(const realVector& kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		utils::makeError("angularDependence::eval(...)","Number of kinematic variables does not match (" + std::to_string(kin.size()) + " != " + std::to_string(_kinSignature->nKin()) + "). Returning zero.");
		return complex(0.,0.);
	}
	utils::makeError("angularDependence::eval(...)","Calling base-class method. Returning zero.");
	return complex(0.,0.);
}

sameMassZeroS::sameMassZeroS(double fsMass) : angularDependence(std::make_shared<kinematicSignature>(1), "sameMassZeroS"), _fsMass(fsMass) {}

complex sameMassZeroS::eval(const realVector& kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		utils::makeError("sameMassZeroS::eval(...)","Number of kinematic variables does not match (" + std::to_string(kin.size()) + " != " + std::to_string(_kinSignature->nKin()) + "). Returning zero.");
		return complex(0.,0.);
	}
//	std::cout << "Called sameMassZeroS with " << kin[0] << " " << kin[1] << " " << kin[2] << " giving 1." << std::endl;
	return complex(1.,0.);
}

sameMassOneP::sameMassOneP(double fsMass) : angularDependence(std::make_shared<kinematicSignature>(1), "sameMassOneP"), _fsMass(fsMass) {}

complex sameMassOneP::eval(const realVector& kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		utils::makeError("sameMassOneP::eval(...)","Number of kinematic variables does not match (" + std::to_string(kin.size()) + " != " + std::to_string(_kinSignature->nKin()) + "). Returning zero.");
		return complex(0.,0.);
	}
	double retVal = (1 + (kin[1] - _fsMass*_fsMass)/kin[0]) * (kin[0] - kin[1] - 2*kin[2] + 3*_fsMass*_fsMass);
//	std::cout << "Called sameMassOneP with " << kin[0] << " " << kin[1] << " " << kin[2] << " giving " << retVal << std::endl;
	return complex(retVal, 0.);
}

sameMassTwoD::sameMassTwoD(double fsMass) : angularDependence(std::make_shared<kinematicSignature>(1), "sameMassTwoD"), _fsMass(fsMass) {}

complex sameMassTwoD::eval(const realVector& kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		utils::makeError("sameMassTwoD::eval(...)","Number of kinematic variables does not match (" + std::to_string(kin.size()) + " != " + std::to_string(_kinSignature->nKin()) + "). Returning zero.");
		return complex(0.,0.);
	}
//	std::cerr << "sameMassTwoD::eval(...): ERROR: Calculations for this wave are not done yet. Returning zero" << std::endl;
	const double s   = kin[0];
	const double s12 = kin[1];
	const double s13 = kin[2];
	const double m2  = _fsMass*_fsMass;
//	double retVal    = (2*m2*pow(m2-s,2) + (-3* m2*m2 - 6 *m2 * s + s*s)*s12 + 2*(3*m2 - s)*s12*s12 + 7*pow(s12,3))*(m2*m2 - 2*s*s - 2*s*s12 + s12*s12 - 2*m2 * (s+s12))/(648 * s * s * s12);
	double retVal = - (m2*m2 - 2*s*s - 2*s*s12 + s12*s12 - 2*m2*(s+s12));
	retVal *= 2*pow(m2,3) + m2*m2*(-4*s + 9*s12) + s12*(pow(s-s12,2) + 6*(-s+s12)*s13 + 6*s13*s13) + 2*m2*(s*s + 3*s*s12 - 3*s12*(s12 + 3*s13));
	retVal /= 648*s*s*s12;
	return complex (retVal, 0.);
}

sameMassZeroSnonRelativistic::sameMassZeroSnonRelativistic(double fsMass) : angularDependence(std::make_shared<kinematicSignature>(1), "sameMassZeroSnonRelativistic"), _fsMass(fsMass) {}

complex sameMassZeroSnonRelativistic::eval(const realVector& kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		utils::makeError("sameMassZeroSnonRelativistic::eval(...)","Number of kinematic variables does not match (" + std::to_string(kin.size()) + " != " + std::to_string(_kinSignature->nKin() ) +  "). Returning zero.");
		return complex(0.,0.);
	}
//	std::cout << "Called sameMassZeroS with " << kin[0] << " " << kin[1] << " " << kin[2] << " giving 1." << std::endl;
	return complex(1.,0.);
}

sameMassOnePnonRelativistic::sameMassOnePnonRelativistic(double fsMass) : angularDependence(std::make_shared<kinematicSignature>(1), "sameMassOnePnonRelativistic"), _fsMass(fsMass) {}

complex sameMassOnePnonRelativistic::eval(const realVector& kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		utils::makeError("sameMassOnePnonRelativistic::eval(...)","Number of kinematic variables does not match (" + std::to_string(kin.size()) + " != " + std::to_string(_kinSignature->nKin() ) +  "). Returning zero.");
		return complex(0.,0.);
	}
//	std::cerr << "sameMassTwoD::eval(...): ERROR: Calculations for this wave are not done yet. Returning zero" << std::endl;
	const double s   = kin[0];
	const double s12 = kin[1];
	const double s13 = kin[2];
	const double m2  = _fsMass*_fsMass;

	double retVal    = -1./4.; // Use definitions as in the paper
//	double retVal    = - pow(s12/s,.5)/4.; //  Use definitions as by Dima
	retVal *= s - s12 - 2*s13 + 3*m2;
	return complex (retVal, 0.);
}

sameMassTwoDnonRelativistic::sameMassTwoDnonRelativistic(double fsMass) : angularDependence(std::make_shared<kinematicSignature>(1), "sameMassTwoDnonRelativistic"), _fsMass(fsMass) {}

complex sameMassTwoDnonRelativistic::eval(const realVector& kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		utils::makeError("sameMassTwoDnonRelativistic::eval(...)","Number of kinematic variables does not match (" + std::to_string(kin.size()) + " != " + std::to_string(_kinSignature->nKin() ) +  "). Returning zero.");
		return complex(0.,0.);
	}
//	std::cerr << "sameMassTwoD::eval(...): ERROR: Calculations for this wave are not done yet. Returning zero" << std::endl;
	const double s   = kin[0];
	const double s12 = kin[1];
	const double s13 = kin[2];
	const double m2  = _fsMass*_fsMass;
	double retVal = 0.;
	retVal += -(-4*m2 + s12)*(m2*m2 + pow(s-s12,2) - 2* m2 * (s+s12));
	retVal += 3*s12*pow(3*m2 + s - s12 - 2*s13,2);
	retVal /= 48*s;
	return complex (retVal, 0.);
}

ratioOfDependences::ratioOfDependences(std::shared_ptr<angularDependence> numerator, std::shared_ptr<angularDependence> denominator) : 
	angularDependence(numerator->kinSignature(), "ratioOfDependences"), _numerator(numerator), _denominator(denominator) 
{
		if (not (*(numerator->kinSignature()) == *(denominator->kinSignature()))) {
			utils::makeError("ratioOfDependences::ratioOfDependences(...)","Different kinematic signatures found.");
		}
}

complex ratioOfDependences::eval(const realVector& kin) const {
	complex den = _denominator->eval(kin);
	if (den == complex(0.,0.) ) {
		return complex (0.,0.);
	}
	return _numerator->eval(kin)/den;
}

arbitraryMass_S_nonRelativistic::arbitraryMass_S_nonRelativistic(size_t isobarIndex, realVector fsMasses) : angularDependence(std::make_shared<kinematicSignature>(2), "arbitraryMass_S_nonRelativistic"), _isobarIndex(isobarIndex), _fsMasses(fsMasses) {
	if (_fsMasses.size() != 3) {
		utils::makeError("arbitraryMass_S_nonRelativistic(...)","Number of final state paricles masses has to be 3");
		throw;
	}
	if (isobarIndex != 12 and isobarIndex != 13 and isobarIndex != 23) {
		utils::makeError("arbitraryMass_S_nonRelativistic(...)","None of the three possible isobar masses match (12, 13 and 23) the given value:" +std::to_string(isobarIndex ) + ".");
		throw;
	}

};

complex arbitraryMass_S_nonRelativistic::eval(const realVector& kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		utils::makeError("arbitraryMass_S_nonRelativistic::eval(...)","Number of kinematic variables does not match (" + std::to_string(kin.size()) + " != " + std::to_string(_kinSignature->nKin() ) +  "). Returning zero.");
		return complex(0.,0.);
	}
	return complex(1.,0.);
}

arbitraryMass_P_nonRelativistic::arbitraryMass_P_nonRelativistic(size_t isobarIndex, realVector fsMasses) : angularDependence(std::make_shared<kinematicSignature>(2), "arbitraryMass_S_nonRelativistic"), _isobarIndex(isobarIndex), _fsMasses(fsMasses) {
	if (_fsMasses.size() != 3) {
		utils::makeError("arbitraryMass_P_nonRelativistic(...)","Number of final state paricles masses has to be 3.");
		throw;
	}
	if (isobarIndex != 12 and isobarIndex != 13 and isobarIndex != 23) {
		utils::makeError("arbitraryMass_P_nonRelativistic(...)","None of the three possible isobar masses match (12, 13 and 23) the given value:" + std::to_string(isobarIndex) + ".");
		throw;
	}

};

complex arbitraryMass_P_nonRelativistic::eval(const realVector& kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		utils::makeError("arbitraryMass_P_nonRelativistic::eval(...)","Number of kinematic variables does not match (" + std::to_string(kin.size()) + " != " + std::to_string(_kinSignature->nKin() ) +  "). Returning zero.");
		return complex(0.,0.);
	}
	const double& s   = kin[0];
	const double& s12 = kin[1];
	const double& s13 = kin[2];

	const double& m1 = _fsMasses[0];
	const double& m2 = _fsMasses[1];
	const double& m3 = _fsMasses[2];

	const double s23 = s + m1*m1 + m2*m2 + m3*m3 - s12 - s13;
	if (_isobarIndex == 12) {
		return (s23 - s13 - (s - m3*m3) * (m1*m1 - m2*m2)/s12)/4;
	}
	if (_isobarIndex == 13) {
		return (s12 - s23 - (s - m2*m2) * (m3*m3 - m1*m1)/s13)/4;
	}
	if (_isobarIndex == 23) {
		return (s13 - s12 - (s - m1*m1) * (m2*m2 - m3*m3)/s23)/4;
	}
	utils::makeError("arbitraryMass_P_nonRelativistic::eval(...)","None of the isobarCombinations matched. Return 0.");
	return complex(0.,0.);
}

BELLE_S::BELLE_S(size_t isobarIndex, const realVector& fsMasses, std::shared_ptr<kinematicSignature> kinSignature) : angularDependence(kinSignature, "BELLE_S"), _isobarIndex(isobarIndex) {
	if (isobarIndex != 12 and isobarIndex != 13 and isobarIndex != 23) {
		utils::makeError("BELLE_S(...)","None of the three possible isobar masses match (12, 13 and 23) the given value:" + std::to_string(isobarIndex ) + ".");
		throw;
	}
	if (fsMasses.size() != 3) {
		utils::makeError("BELLE_S(...)","Number of final state masses given is not three.");
		throw;
	}
}

complex BELLE_S::eval(const realVector &kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		utils::makeError("BELLE_S::eval(...)","Number of kinematic variables does not match (" + std::to_string(kin.size()) + " != " + std::to_string(_kinSignature->nKin() ) +  "). Returning zero.");
		return complex(0.,0.);
	}
	return complex(1.,0.);
}

BELLE_P::BELLE_P(size_t isobarIndex, const realVector& fsMasses, std::shared_ptr<kinematicSignature> kinSignature) : angularDependence(kinSignature, "BELLE_P"), _isobarIndex(isobarIndex), _fsMasses() {
	if (isobarIndex != 12 and isobarIndex != 13 and isobarIndex != 23) {
		utils::makeError("BELLE_P(...)","None of the three possible isobar masses match (12, 13 and 23) the given value:"  + std::to_string(isobarIndex ) + ".");
		throw;
	}
	if (fsMasses.size() != 3) {
		utils::makeError( "BELLE_P(...)","Number of final state masses given is not three.");
		throw;
	}
	_fsMasses = fsMasses;
}

complex BELLE_P::eval(const realVector &kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		utils::makeError("BELLE_P::eval(...)","Number of kinematic variables does not match (" + std::to_string(kin.size()) + " != " + std::to_string(_kinSignature->nKin() ) +  "). Returning zero.");
		return complex(0.,0.);
	}

	double mD2  = kin[0];
	double mA2  = 0.;
	double mB2  = 0.;
	double mC2  = 0.;
	double mAB2 = 0.;
	double mAC2 = 0.;
	double mBC2 = 0.;

	double sKpiWrong = mD2 + _fsMasses[0]*_fsMasses[0] +  _fsMasses[1]*_fsMasses[1]+ _fsMasses[2]*_fsMasses[2]- kin[1] - kin[2];
	if (_isobarIndex == 12) { // ((piRight, Ks), piWrong)
		mA2 = _fsMasses[0]*_fsMasses[0];
		mB2 = _fsMasses[1]*_fsMasses[1];
		mC2 = _fsMasses[2]*_fsMasses[2];

		mAB2 = kin[1];
		mAC2 = kin[2];
		mBC2 = sKpiWrong;
	} else if (_isobarIndex == 13) { // ((piRight, piWrong), Ks)
		mA2 = _fsMasses[2]*_fsMasses[2];
		mB2 = _fsMasses[0]*_fsMasses[0];
		mC2 = _fsMasses[1]*_fsMasses[1];

		mAB2 = kin[2];
		mAC2 = sKpiWrong;
		mBC2 = kin[1];
	} else if (_isobarIndex == 23) { // ((Ks, piWrong), piRight)
		mA2 = _fsMasses[1]*_fsMasses[1];
		mB2 = _fsMasses[2]*_fsMasses[2];
		mC2 = _fsMasses[0]*_fsMasses[0];

		mAB2 = sKpiWrong;
		mAC2 = kin[1];
		mBC2 = kin[2];
	}
	double retVal = mAC2 - mBC2 + (mD2-mC2)*(mB2-mA2)/(mAB2);
	return complex(retVal, 0.);
}

bool BELLE_P::setFSmasses(const realVector newMasses) {
	if (newMasses.size() != 3) {
		utils::makeError("BELLE_P::setFSmasses(...)","Number of masses given is not three: " + std::to_string( newMasses.size())+ ".");
		return false;
	}
	_fsMasses[0] = newMasses[0];
	_fsMasses[1] = newMasses[1];
	_fsMasses[2] = newMasses[2];
	return true;
}

BELLE_D::BELLE_D(size_t isobarIndex, const realVector& fsMasses, std::shared_ptr<kinematicSignature> kinSignature) : angularDependence(kinSignature, "BELLE_D"), _isobarIndex(isobarIndex) {
	if (isobarIndex != 12 and isobarIndex != 13 and isobarIndex != 23) {
		utils::makeError("BELLE_D(...)","None of the three possible isobar masses match (12, 13 and 23) the given value:"  + std::to_string(isobarIndex )+".");
		throw;
	}
	if (fsMasses.size() != 3) {
		utils::makeError("BELLE_D(...)","Number of final state masses given is not three.");
		throw;
	}
	_fsMasses = fsMasses;
}

complex BELLE_D::eval(const realVector &kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		utils::makeError("BELLE_D::eval(...)","Number of kinematic variables does not match (" + std::to_string(kin.size()) + " != " + std::to_string(_kinSignature->nKin() ) +  "). Returning zero.");
		return complex(0.,0.);
	}

	double mD2  = kin[0];
	double mA2  = 0.;
	double mB2  = 0.;
	double mC2  = 0.;
	double mAB2 = 0.;
	double mAC2 = 0.;
	double mBC2 = 0.;

	double sKpiWrong = mD2 + _fsMasses[0]*_fsMasses[0] +  _fsMasses[1]*_fsMasses[1]+ _fsMasses[2]*_fsMasses[2]- kin[1] - kin[2];
	if (_isobarIndex == 12) { // ((piRight, Ks), piWrong)
		mA2 = _fsMasses[0]*_fsMasses[0];
		mB2 = _fsMasses[1]*_fsMasses[1];
		mC2 = _fsMasses[2]*_fsMasses[2];

		mAB2 = kin[1];
		mAC2 = kin[2];
		mBC2 = sKpiWrong;
	} else if (_isobarIndex == 13) { // ((piWrong, piRight), Ks)
		mA2 = _fsMasses[2]*_fsMasses[2];
		mB2 = _fsMasses[0]*_fsMasses[0];
		mC2 = _fsMasses[1]*_fsMasses[1];

		mAB2 = kin[2];
		mAC2 = sKpiWrong;
		mBC2 = kin[1];
	} else if (_isobarIndex == 23) { // ((Ks, piWrong), piRight)
		mA2 = _fsMasses[1]*_fsMasses[1];
		mB2 = _fsMasses[2]*_fsMasses[2];
		mC2 = _fsMasses[0]*_fsMasses[0];

		mAB2 = sKpiWrong;
		mAC2 = kin[1];
		mBC2 = kin[2];
	}
	double retVal = mAB2 - 2*mD2 - 2*mC2 + pow(mD2 - mC2,2)/mAB2;
	retVal       *= mAB2 - 2*mA2 - 2*mB2 + pow(mA2 - mB2,2)/mAB2;
	retVal       /= -3;
	retVal       += pow(mBC2 - mAC2 + (mD2 - mC2)*(mA2 - mB2)/mAB2,2);
	return complex(retVal, 0.);
}

bool BELLE_D::setFSmasses(const realVector newMasses) {
	if (newMasses.size() != 3) {
		utils::makeError("BELLE_D::setFSmasses(...)","Number of masses given is not three: " + std::to_string( newMasses.size() )+".");
		return false;
	}
	_fsMasses[0] = newMasses[0];
	_fsMasses[1] = newMasses[1];
	_fsMasses[2] = newMasses[2];
	return true;
}

binnedDalitzAmplitude::binnedDalitzAmplitude(std::shared_ptr<angularDependence> angAmp, const realVector& binningXaxis, const realVector& binningYaxis) : 
	angularDependence(angAmp->kinSignature(), angAmp->name() + std::string("_binned")), _angAmp(angAmp), _binningX(binningXaxis), _binningY(binningYaxis) {}

complex binnedDalitzAmplitude::eval(const realVector& kin) const {
	realVector newKin = getBinnedKin(kin);
	if (newKin.size() == 0) {
		utils::makeWarning("binnedDalitzAmplitude::eval(...)","Got empty kinematic vector, returning zero.");
		return complex(0.,0.);
	}
	return _angAmp->eval(newKin);
}

realVector binnedDalitzAmplitude::getBinnedKin(const realVector& kin) const {
	realVector retVal(3, kin[0]); 
	for (size_t i = 0; i < _binningX.size() - 1; ++i) {
		if (_binningX[i] <= kin[1] && _binningX[i+1] > kin[1]) {
			retVal[1] = (_binningX[i] + _binningX[i+1])/2;
			break;
		}
	}
	for (size_t i = 0; i < _binningY.size() - 1; ++i) {
		if (_binningY[i] <= kin[2] && _binningY[i+1] > kin[2]) {
			retVal[2] = (_binningY[i] + _binningY[i+1])/2;
			break;
		}
	}
	if (kin[1] == kin[0]) { // This check works, since kin[1] always has to be smaller for a kinematically valid bin
		utils::makeWarning("binnedDalitzAmplitude::getBinnedKin(...)","No bin for kin[1] found, returning empty vector (kin[1] = " + std::to_string(kin[1]) + ").");
		return realVector();
	}
	if (kin[2] == kin[0]) { // This check works, since kin[1] always has to be smaller for a kinematically valid bin
		utils::makeWarning("binnedDalitzAmplitude::getBinnedKin(...)","No bin for kin[2] found, returning empty vector (kin[2] = " + std::to_string(kin[2]) + ").");
		return realVector();
	}

	return retVal;
}
