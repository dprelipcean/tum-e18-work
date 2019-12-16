#include<iostream>
#include<fstream>
#include <math.h>

#include"phaseSpaceFunctions.h"
realPhaseSpaceBase::realPhaseSpaceBase() {}

twoBodyPhaseSpace::twoBodyPhaseSpace(double m1, double m2) : realPhaseSpaceBase(), _m1(m1), _m2(m2), _sTh((m1+m2)*(m1+m2)) {}

double twoBodyPhaseSpace::eval(double s) const {
	if (s < _sTh) {
		return 0.;
	}	
	return pow(1. -_sTh/s,.5);
}

fourPiPhaseSpace::fourPiPhaseSpace(std::string valueFileName, size_t nStep, double sMin, double sMax, double mPi, double mRho, double gamma, double sanityDelta) : 
	realPhaseSpaceBase(), _nStep(nStep), _sMin(sMin), _sMax(sMax), _step((sMax-sMin)/nStep), _mPi(mPi), _mRho(mRho), _gamma(gamma), 
	_valueFileName(valueFileName), _values(nStep+1,0.)
{
	std::ifstream fin(_valueFileName.c_str(), std::ifstream::in);
	size_t nStepLoad = 0;
	if (!(fin>>nStepLoad)) {
		std::cout << "fourPiPhaseSpace(...): ERROR: Could not load 'nStep' from configuration file" << std::endl; 
		throw;
	}
	if (nStepLoad != _nStep) {
		std::cout << "fourPiPhaseSpace(...): ERROR: _nStep does not match the value in the configuration file" << std::endl;
		throw;
	}
	double loadValue = 0.;
	if (!(fin>>loadValue)) {
		std::cout << "fourPiPhaseSpace(...): ERROR: Could not load 'sMin' from configuration file" << std::endl;
		throw;
	}
	if (loadValue != _sMin) {
		std::cout << "fourPiPhaseSpace(...): ERROR: _sMin does not match the value in the configuration file" << std::endl;
		throw;
	}
	if (!(fin>>loadValue)) {
		std::cout << "fourPiPhaseSpace(...): ERROR: Could not load 'sMax' from configuration file" << std::endl;
		throw;
	}
	if (loadValue != _sMax) {
		std::cout << "fourPiPhaseSpace(...): ERROR: _sMax does not match the value in the configuration file" << std::endl;
		throw;
	}
	if (!(fin>>loadValue)) {
		std::cout << "fourPiPhaseSpace(...): ERROR: Could not load 'mPi' from configuration file" << std::endl;
		throw;
	}
	if (loadValue != _mPi) {
		std::cout << "fourPiPhaseSpace(...): ERROR: _mPi does not match the value in the configuration file" << std::endl;
		throw;
	}
	if (!(fin>>loadValue)) {
		std::cout << "fourPiPhaseSpace(...): ERROR: Could not load 'mRho' from configuration file" << std::endl;
		throw;
	}
	if (loadValue != _mRho) {
		std::cout << "fourPiPhaseSpace(...): ERROR: _mRho does not match the value in the configuration file" << std::endl;
		throw;
	}
	if (!(fin>>loadValue)) {
		std::cout << "fourPiPhaseSpace(...): ERROR: Could not load 'gamma' from configuration file" << std::endl;
		throw;
	}
	if (loadValue != _gamma) {
		std::cout << "fourPiPhaseSpace(...): ERROR: _gamma does not match the value in the configuration file" << std::endl;
		throw;
	}
	size_t count = 0;
	while (fin>>loadValue) {
		if (count == _nStep+1) {
			std::cout << "fourPiPhaseSpace(...): ERROR: _nStep+1 (=" << _nStep+1 << ") range exceeded" << std::endl;
			throw;
		}
		_values[count] = loadValue;
		++count;
	}
	if (count != _nStep+1) {
		std::cout << "fourPiPhaseSpace(...): ERROR: number of loaded points does not match:" << count << "; should be " << _nStep+1 << std::endl;
		throw;
	}
	if (!sanityCheck(sanityDelta)) {
		std::cout << "fourPiPhaseSpace(...): ERROR: Loaded values are insane" << std::endl;
		throw;
	}
}

double fourPiPhaseSpace::eval(double s) const { 
	double retVal = 0.;
	if (s > _sMin) {
		if (s < _sMax) {
			retVal =  interpolate_rho51(s);
		} else {
			retVal = rho52(s);
		}
	}
	return retVal;
}

double fourPiPhaseSpace::rho52(double s) const {
	return pow(1.-_sMin/s, .5);
}

double fourPiPhaseSpace::interpolate_rho51(double s) const {
	size_t i = (size_t)((s-_sMin)/_step);
	double x = (s - _sMin - i*_step)/_step;
	return (1.-x)*_values[i] + x*_values[i+1];
}

bool fourPiPhaseSpace::sanityCheck(double sanityDelta) const {
	double del = rho52(_sMax) - interpolate_rho51(_sMax);
	if (del*del > sanityDelta*sanityDelta) {
		std::cout << "fourPiPhaseSpace::sanityCheck(...): ERROR: Function not continuous at _sMax = " << _sMax << std::endl;
		return false;
	}
	del = 16*_mPi*_mPi - _sMin;
	if (del*del > sanityDelta*sanityDelta) {
		std::cout << "fourPiPhaseSpace::sanityCheck(...): ERROR: _sMin != 16*mPi*mPi: Interpolation does not start at threshold" << std::endl;
		return false;
	}
	for (size_t i = 0; i < _nStep; ++i) {
		if (_values[i+1] < _values[i]) {
			std::cout << "fourPiPhaseSpace::sanityCheck(...): ERROR: Function not increasing at i = " << i << std::endl;
			return false;
		}
		if (_values[i] < 0.) {
			std::cout << "fourPiPhaseSpace::sanityCheck(...): ERROR: Subzero value at i = " << i << std::endl;
			return false;
		}
	}
	return true;
}

