#include"modelGenerator.h"
#include<iostream>

modelGenerator::modelGenerator(std::shared_ptr<amplitude> model, std::shared_ptr<generator> gen, std::shared_ptr<efficiencyFunction> efficiency) : 
	modelGenerator::modelGenerator(std::vector<std::shared_ptr<amplitude> >(1, model), gen, efficiency) {};

modelGenerator::modelGenerator(std::vector<std::shared_ptr<amplitude> > model, std::shared_ptr<generator> gen, std::shared_ptr<efficiencyFunction> efficiency):
	_maxFail(1000), _failCount(0), _kinSignature(gen->kinSignature()), _generator(gen), _model(model), _efficiency(efficiency) {
	if (model.size() == 0) {
		std::cout << "modelGenerator::modelGenerator(...): ERROR: No amplitude given" << std::endl;
		throw;
	}
	if (not(*_kinSignature == *(_model[0]->kinSignature()))) {
		std::cout << "modelGenerator::modelGenerator(...): ERROR: Kinematic signatures does not match with generator" << std::endl;
		throw;
	}
	if (_efficiency) {
		if (not(*_kinSignature == *(_efficiency->kinSignature()))) {
			std::cout << "modelGenerator::modelGenerator(...): ERROR: Kinematic signatures does not match with efficiency" << std::endl;
			throw;
		}
	}
}

std::pair<std::pair<size_t, double>, std::vector<std::vector<double > > > modelGenerator::burnIn(size_t nPoints) const {
	double maxWeight = 0.;
	realVector weights(nPoints);
	realMatrix points(nPoints);
	for (size_t p = 0; p < nPoints; ++p) {
		realVector kin = _generator->generate();
		double intens = 0;
		for (const std::shared_ptr<amplitude> ampl: _model) {
			intens += ampl->intens(kin);
		}
		weights[p] = intens;
		points[p]  = kin;
		if (intens > maxWeight) {
			maxWeight = intens;
		}
	}
	realMatrix deWeighted(nPoints);
	size_t remaining = 0;
	for (size_t p = 0; p < nPoints; ++p) {
		if (weights[p] > utils::random() * maxWeight) {
			deWeighted[remaining] = points[p];
			++remaining;
		}
	}
	return std::pair<std::pair<size_t, double>, realMatrix >(std::pair<size_t, double>(remaining, maxWeight), deWeighted);
}

realMatrix modelGenerator::generateDataPoints(size_t nPoints, size_t nBurnIn) const {
	double maxWeight = 0.;
	size_t found     = 0;
	realMatrix retVal(nPoints);
	{
		std::pair<std::pair<size_t, double>, std::vector<std::vector<double > > > burn = burnIn(nBurnIn);
		found     = std::min(burn.first.first, nPoints); 
		maxWeight = burn.first.second;
		for (size_t i = 0; i < found; ++i) {
			retVal[i] = burn.second[i];
		}
	}
	utils::makeInfo("modelGenerator::generateDataPoints(...)","Got " + std::to_string(found) + " of " + std::to_string(nPoints) + " events after burn in.");
	while (found < nPoints) {
		realVector kin = _generator->generate();
		double intens = 0.;
		for (std::shared_ptr<amplitude> ampl : _model) {
			intens += ampl->intens(kin);
		}
		if (intens > utils::random() * maxWeight) {
			retVal[found] = kin;
			++found;
			_failCount = 0;
		} else {
			++_failCount;
			if (_failCount > _maxFail) {
				utils::makeError("modelGenerator::generateDataPoints(...)","Over " + std::to_string(_maxFail) +  " failed attempts. Aborting.");
				throw;
			}
		}
		if (intens > maxWeight) {
			realMatrix newVal(nPoints);
			size_t count = 0;
			for (size_t i = 0; i < found; ++i) {
				if (maxWeight > utils::random() * intens) {
					newVal[count] = retVal[i];
					++count;
				}
			}
//			std::cout << (double)count/found << " vs " << maxWeight/intens << std::endl;
			found     = count;
			maxWeight = intens;
			retVal    = newVal;
			utils::makeInfo("modelGenerator::generateDataPoints(...)","Reweighting: " + std::to_string(count) + " events left.");
		}
	}
	return retVal;
}

std::pair<double, realVector > modelGenerator::getSinglePoint() const { 
	realVector kin = _generator->generate();
	double intens = 0.;
	for (std::shared_ptr<amplitude> ampl : _model) {
		intens += ampl->intens(kin);
	}
	if (_efficiency) {
		intens *= _efficiency->eval(kin);
	}
	return std::pair<double, realVector > (intens, kin);
}
