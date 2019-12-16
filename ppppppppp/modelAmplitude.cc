#include"modelAmplitude.h"
#include"utils.h"
#include<iostream>
modelAmplitude::modelAmplitude(complexVector transitionAmplitudes, std::vector<std::shared_ptr<amplitude> > amplitudes, realVector normalizations, std::string name) : 
amplitude(amplitudes.at(0)->kinSignature(), name)
{
	if (amplitudes.size() == 0) {
		utils::makeError("modelAmplitude::modelAmplitude(...)","No amplitudes given.");
		throw;
	}
	_nAmpl = amplitudes.size();
	_kinSignature = amplitudes[0]->kinSignature();
	for (std::shared_ptr<amplitude> a : amplitudes) {
		if (not ( *(a->kinSignature()) == *_kinSignature)) {
			utils::makeError("modelAmplitude::modelAmplitude(...)","Different kinematic signatures encountered.");
			throw;
		}
	}
	_amplitudes = amplitudes;
	if (amplitudes.size() != transitionAmplitudes.size()) {
		utils::makeError("modelAmplitude::modelAmplitude(...)","Size of transitions amplitudes ("+std::to_string(transitionAmplitudes.size())+") does not match (should be " + std::to_string(amplitudes.size()) + ").");
		throw;
	}
	_transitionAmplitudes = transitionAmplitudes;
	if (amplitudes.size() != normalizations.size()) {
		utils::makeError("modelAmplitude::modelAmplitude(...)","Size of normalizations ("+std::to_string(normalizations.size())+") does not match (should be "+std::to_string(amplitudes.size())+".");
		throw;		
	}
	_normalizations = normalizations;
}

complex modelAmplitude::eval(const realVector& kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		utils::makeError("modelAmplitude::ampl(...)","Number of kinematic variables does not match (" +std::to_string(kin.size()) + " != "+ std::to_string(_kinSignature->nKin()) + "). Returning zero.");
		return complex(0.,0.);
	}

	complex retVal(0.,0.);
	for (size_t a = 0; a < _nAmpl; ++a) {
//		std::cout << a << " :: :: :::::: " << std::endl;
//		std::cout << _transitionAmplitudes[a] << " * " <<_amplitudes[a]->eval(kin) << " * "  <<  _normalizations[a] << " = " << _transitionAmplitudes[a] * _amplitudes[a]->eval(kin) * _normalizations[a];
		retVal += _transitionAmplitudes[a] * _amplitudes[a]->eval(kin) * _normalizations[a];
//		std::cout << " -> "  << retVal << std::endl;
	}
//	std::cout << " - E - V - A - L - - - F - I - N - I - S - H - E - D - " << std::endl;
	return retVal;
}

bool modelAmplitude::setTransitionAmplitude(size_t n, complex amp) {
	if (n >= _nAmpl) {
		return false;
	}
	_transitionAmplitudes[n] = amp;
	return true;
}
