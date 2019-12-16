#include"logLikelihood.h"

#include<string>
#include<limits>
#include<cmath>        // std::abs

#include<nlopt.hpp>

#include"omp.h"
#include"utils.h"
//#include<fstream>

double ff_nlopt(const realVector &x, realVector &grad, void* f_data) {
	logLikelihood* inst = reinterpret_cast<logLikelihood*>(f_data);
	double retVal   = inst->nloptCall(x,grad);
	return retVal;
}

logLikelihoodBase::logLikelihoodBase (size_t nAmpl, std::shared_ptr<kinematicSignature> kinSignature, std::shared_ptr<integrator> integral) : 
	_extended(true), _kinSignature(kinSignature), _nCalls(0), _nCallsPrint(10000), _nAmpl(nAmpl), _nPoints(0), _nScale(0), _numDelta(1.e-5), _integral(integral), _fixedParameters(), _copyParameters(), _directFitParameters()
 	, _trackFit(false), _storeParams() 
{
	if (_nAmpl == 0) {
		utils::makeError("logLikelihoodBase::logLikelihoodBase(...):","No amplitudes given.");
		throw;
	}
	if (_nAmpl != _integral->nAmpl()) {
		utils::makeError("logLikelihoodBase::logLikelihoodBase(...)","Number of amplitudes does not match.");
		throw;
	}
	if (not (*(_integral->kinSignature()) == *_kinSignature)) {
		utils::makeError("logLikelihoodBase::logLikelihoodBase(...)","Kinematic signatures doe not match.");
		throw;	
	}
	if (not _integral->isIntegrated()) {
		if (not _integral->integrate()) {
			utils::makeError("logLikelihoodBase::logLikelihoodBase(...)","Could not obtain integral.");
			throw;
		}
	}
}

std::pair<double, complexVector > logLikelihoodBase::fitNlopt(const realVector& parameters) {
	if (parameters.size() != getNparFit()) {
		utils::makeError("logLikelihoodBase::fitNlopt(...)","Wrong numer of parameters: parameters.size() = " 
		                                                    + std::to_string(parameters.size()) + " equal getNparFit() = " + std::to_string(getNparFit()));
		return std::pair<double, complexVector >(std::numeric_limits<double>::infinity(), complexVector());
	}
	realVector startPars = parameters;

//	nlopt::opt opt = nlopt::opt(nlopt::LD_LBFGS, getNparFit());
	nlopt::opt opt = nlopt::opt(nlopt::LD_SLSQP, getNparFit());
//	opt.set_ftol_abs(1.E-14);
	opt.set_maxeval(1000000);
	opt.set_min_objective(ff_nlopt, this);

	double bestVal;
	opt.set_vector_storage(100);
	utils::makeInfo("logLikelihoodBase::fitNlopt(...)","nlopt vector storage is: " + std::to_string(opt.get_vector_storage()) + ".");
	utils::makeInfo("logLikelihoodBase::fitNlopt(...)","nlopt ftol rel is: " + std::to_string(opt.get_ftol_rel()) + ".");
	utils::makeInfo("logLikelihoodBase::fitNlopt(...)","nlopt ftol abs is: " + std::to_string(opt.get_ftol_abs()) + ".");

	int optimization_result  = opt.optimize(startPars, bestVal);
	utils::makeInfo("logLikelihoodBase::fitNlopt(...)","Optimization status: " +std::to_string(optimization_result) + ":");
	if (optimization_result == 1) {
		std::cout << " - nlopt status: NLOPT_SUCCESS" << std::endl;
	} else if ( optimization_result == 2 ) {
		std::cout << " - nlopt status: NLOPT_STOPVAL_REACHED" << std::endl;
	} else if ( optimization_result == 3 ) {
		std::cout << " - nlopt status: NLOPT_FTOL_REACHED" << std::endl;
	} else if ( optimization_result == 4 ) {
		std::cout << " - nlopt status: NLOPT_XTOL_REACHED" << std::endl;
	} else if ( optimization_result == 5 ) {
		std::cout << " - nlopt status: NLOPT_MAXEVAL_REACHED" << std::endl;
	} else if ( optimization_result == 6 ) {
		std::cout << " - nlopt status: NLOPT_MAXTIME_REACHED" << std::endl;
	} else if ( optimization_result == -1 ) {
		std::cout << " - nlopt status: NLOPT_FAILURE" << std::endl;
	} else if ( optimization_result == -2 ) {
		std::cout << " - nlopt status: NLOPT_INVALID_ARGS" << std::endl;
	} else if ( optimization_result == -3 ) {
		std::cout << " - nlopt status: NLOPT_OUT_OF_MEMORY" << std::endl;
	} else if ( optimization_result == -4 ) {
		std::cout << " - nlopt status: NLOPT_ROUNDOFF_LIMITED" << std::endl;
	} else if ( optimization_result == -5 ) {
		std::cout << " - nlopt status: NLOPT_FORCED_STOP" << std::endl;
	} else {
		std::cout << " - Unknown nlopt status: " << optimization_result << std::endl;
	}

	_directFitParameters = startPars;

	complexVector retVal = makeProdAmpsFromFitPars(startPars);
	return std::pair<double, complexVector >(bestVal, retVal);
}

std::pair<double, complexVector > logLikelihoodBase::fitNlopt(const complexVector& parameters) {
	utils::makeError("logLikelihoodBase::fitNlopt(const complexVector& prodAmps)","Starting the fit with production amplitudes is not supported anymore.");
	throw;
	realVector startPars = cutGradient(prodAmpsToFullParams(parameters)); // Use this, since no better option now... (will change the start parameters a bit)
	return fitNlopt(startPars);
}

double logLikelihoodBase::nloptCall(const realVector &x, realVector &grad) const {
//	std::cout << "logLikelihoodBase::nloptCall(...) INFO: Called!" << std::endl;
	if (x.size() != getNpar()) {
		utils::makeError("logLikelihood::nloptCall(...)","Number of amplitudes does not match (" + std::to_string(x.size()) + " vs. " + std::to_string(getNpar()) + ").");
		throw;
	}
	

	if (_trackFit) {
		if (_storeParams.size() <= nCalls()) {
			_storeParams.resize(2*nCalls()+1);
		}
		_storeParams[nCalls()] = realVector(x);
	}
	complexVector prodAmps = fullParamsToProdAmps(getFullParameters(x));
//	complexVector prodAmps = makeProdAmpsFromFitPars(x);
	double retVal = eval(prodAmps);
	if (grad.size() > 0) {
		realVector fullGradient = makeFinalGradient(x,Deval(prodAmps));
		for (size_t g = 0; g < grad.size(); ++g) {
			grad[g] = fullGradient[g];
		}		
// //		grad = makeFinalGradient(x,Deval(prodAmps));
// //		if (_ngc) {                                          // DELETE
// //			double delta = 0.00001;                      // DELETE
// //			for (size_t p = 0; p < getNpar(); ++p) {     // DELETE
// //				double dd = delta * _storeParams[p]; // DELETE
// //				_storeParams[p] += dd;               // DELETE
// //				std::cout << "ngc " << p << ": " << (eval(fullParamsToProdAmps(getFullParameters(_storeParams))) - retVal)/dd << " |vs| " << grad[p] << std::endl; // DELETE
// //				_storeParams[p] -= dd;               // DELETE
// //			}                                            // DELETE
// //		}                                                    // DELETE
	}
//	std::cout << "logLikelihoodBase::nloptCall(...) INFO: Getting: " << retVal << std::endl;
	return retVal;
}

realVector logLikelihoodBase::makeFinalGradient(const realVector& params, const realVector& fullGradient) const {

	const size_t nPar         = getNpar();
	const size_t nParTot      = getNparTot();
	const size_t nParNonScale = nParTot - _fixedParameters.size();

	if (params.size() != nPar) {
		utils::makeError("logLikelihoodBase::makeFinalGradient(...)","Wrong parameter size.");
		throw;
	}

	realVector fullParams;
	if (_nScale > 0) {
		fullParams = getFullParameters(params);
	}

	realVector scaleGradients(_nScale, 0.);
	realVector modifiedGradient = fullGradient;
	for (std::tuple<size_t, size_t, int, double> copy : _copyParameters) { // First: modify the entries of the gradient corresponding to the actual parameters & collect the gradient for the scale parameters
		int scl = std::get<2>(copy);
		double scaler = std::get<3>(copy);
		if (scl > -1) {
			scaleGradients[scl] += fullGradient[std::get<0>(copy)] * fullParams[std::get<1>(copy)]*scaler;
			scaler *= params[nParNonScale + scl];
		}
		modifiedGradient[std::get<1>(copy)] += fullGradient[std::get<0>(copy)] * scaler;
	}

	modifiedGradient = cutGradient(modifiedGradient); // Second: Cut away the fixed parameters from the grandient (this leaves the gradient entries for the scale parameters with odd values)
	for (size_t s = 0; s < _nScale; ++s) {
		modifiedGradient[nParNonScale+s] = scaleGradients[s]; // Third: Fix the grandient entries for the scale parameters 
	}
	return modifiedGradient;

}

realMatrix logLikelihoodBase::makeFinalHessian(const realVector& params, const realMatrix& fullHessian) const {

	const size_t nPar = getNpar();
	const size_t nParTot = getNparTot();
	const size_t nParNonScale = nParTot - _fixedParameters.size();

	realVector fullParams;
	if (_nScale > 0) {
		fullParams = getFullParameters(params);
	}
	realMatrix modifiedHessian(fullHessian.size());
// Do eveything for lines first (interprete the hessian as a list of gradients)
	for (size_t p = 0; p < fullHessian.size(); ++p) {
		realVector scaleGradients(_nScale, 0.);
		realVector modifiedGradient = fullHessian[p];
		for (std::tuple<size_t, size_t, int, double> copy : _copyParameters) { // First: modify the entries of the gradient corresponding to the actual parameters & collect the gradient for the scale parameters
			int scl = std::get<2>(copy);
			double scaler = std::get<3>(copy);
			if (scl > -1) {
				scaleGradients[scl] += fullHessian[p][std::get<0>(copy)] * fullParams[std::get<1>(copy)]*scaler;
				scaler *= params[nParNonScale + scl];
			}
			modifiedGradient[std::get<1>(copy)] += fullHessian[p][std::get<0>(copy)] * scaler;
		}
		modifiedGradient = cutGradient(modifiedGradient); // Second: Cut away the fixed parameters from the grandient (this leaves the gradient entries for the scale parameters with odd values)
		for (size_t s = 0; s < _nScale; ++s) {
			modifiedGradient[nParNonScale+s] = scaleGradients[s]; // Third: Fix the grandient entries for the scale parameters 
		}
		modifiedHessian[p] = modifiedGradient;
	}
// The do the same for columns (then turn around this interpretation)
	realMatrix retVal(nPar);
	for (size_t p = 0; p < nPar; ++p) {
		realVector scaleGradients(_nScale, 0.);
		realVector modifiedGradient(nParTot);
		for (size_t i = 0; i < nParTot; ++i) {
		 	modifiedGradient[i] = modifiedHessian[i][p];
		}
		for (std::tuple<size_t, size_t, int, double> copy : _copyParameters) { // First: modify the entries of the gradient corresponding to the actual parameters & collect the gradient for the scale parameters
			int scl = std::get<2>(copy);
			double scaler = std::get<3>(copy);
			if (scl > -1) {
				scaleGradients[scl] += modifiedHessian[std::get<0>(copy)][p] * fullParams[std::get<1>(copy)]*scaler;
				scaler *= params[nParNonScale + scl];
			}
			modifiedGradient[std::get<1>(copy)] += modifiedHessian[std::get<0>(copy)][p] * scaler;
		}
		modifiedGradient = cutGradient(modifiedGradient); // Second: Cut away the fixed parameters from the grandient (this leaves the gradient entries for the scale parameters with odd values)
		for (size_t s = 0; s < _nScale; ++s) {
			modifiedGradient[nParNonScale+s] = scaleGradients[s]; // Third: Fix the grandient entries for the scale parameters 
		}
		retVal[p] = modifiedGradient;
	}
	return retVal;
}


size_t logLikelihoodBase::getNparTot() const {
	size_t retVal = 2*_nAmpl;
	return retVal;
}

size_t logLikelihoodBase::getNpar() const {
	size_t retVal = getNparTot() - _fixedParameters.size() + _nScale;
	return retVal;
}

realVector logLikelihoodBase::getFullParameters(const realVector& params) const {
	const size_t nPar = getNpar();
	const size_t nParTot = getNparTot();
	if (params.size() != nPar) {
		utils::makeError("logLikelihoodBase::getFullParameters(...)","Parameter size mismatch: " + std::to_string(params.size()) + " should be " + std::to_string(nPar) + ".");
		throw;
	}
	realVector retVal(nParTot);
	size_t nextFixed = nParTot;
	if (_fixedParameters.size() > 0) {
		nextFixed = _fixedParameters[0].first;
	}
	size_t fixCount  = 0;
	size_t parCount  = 0;
	for (size_t p = 0; p < nParTot; ++p) {
		if (p == nextFixed) {
			retVal[p] = _fixedParameters[fixCount].second;
			++fixCount;
			if (fixCount == _fixedParameters.size()) {
				nextFixed = nParTot;
			} else {
				nextFixed = _fixedParameters[fixCount].first;
			}
		} else {
			retVal[p] = params[parCount];
			++parCount;
		}
	}
	for (std::tuple<size_t, size_t, int, double> copy : _copyParameters) {
		if (std::isnan(retVal[std::get<0>(copy)])) {
			retVal[std::get<0>(copy)] = 0.;
		}
		double toAdd = retVal[std::get<1>(copy)] * std::get<3>(copy);
		if (std::get<2>(copy) > -1) {
			toAdd *= params[parCount + std::get<2>(copy)];
		}
		retVal[std::get<0>(copy)] += toAdd;
	}
//	parCount += _nScale; // do this, if a third type of parameter should appear
	return retVal;
}

realVector logLikelihoodBase::cutGradient(const realVector& params) const {
	const size_t nPar = getNpar();
	const size_t nParTot = getNparTot();
	if (params.size() != nParTot) {
		utils::makeError("logLikelihoodBase::cutGradient(...)","Parameter size mismatch: " + std::to_string(params.size()) + " should be " + std::to_string(nParTot) + ".");
		throw;
	}
	realVector retVal(nPar);
	size_t nextFixed = nParTot;
	if (_fixedParameters.size() > 0) {
		nextFixed = _fixedParameters[0].first;
	}
	size_t fixCount  = 0;
	size_t parCount  = 0;
	for (size_t p = 0; p < nParTot; ++p) {
		if (p == nextFixed) {
			++fixCount;
			if (fixCount == _fixedParameters.size()) {
				nextFixed = nParTot;
			} else {
				nextFixed = _fixedParameters[fixCount].first;
			}
		} else {
			retVal[parCount] = params[p];
			++parCount;
		}
	}
	std::vector<bool> scalesSet(_nScale, false);
	size_t nScalesSet = 0;
	for (std::tuple<size_t, size_t, int, double> copy : _copyParameters) {
		int scl = std::get<2>(copy);
		if (scl > -1) {
			if (!scalesSet[scl]) {
				double from = params[std::get<1>(copy)];
				if (from == 0.) {
					continue;
				}
				double to   = params[std::get<0>(copy)];
				retVal[parCount+scl] = to/from/std::get<3>(copy);
				scalesSet[scl] = true;
				++nScalesSet;
			}
		}
		if (nScalesSet == _nScale) {
			break;
		}
	}
//	parCount += _nScale; // do this, if a third type of parameter should appear
	return retVal;
}

bool logLikelihoodBase::fixParameter(size_t nPar, double value) {
	if (nPar >= getNparTot()) {
		utils::makeError("logLikelihoodBase::fixParameter(...)","Parameter number too large.");
		return false;
	}
	std::vector<std::pair<size_t, double> > newFixedPars;
	bool appended = false;
	for (const std::pair<size_t, double>& fixed : _fixedParameters) {
		if (fixed.first == nPar) {
			if ((!std::isnan(fixed.second) or !std::isnan(value)) and fixed.second != value) { // NaN are copy parameters. (The same value twice does not couse an error)
				utils::makeError("logLikelihoodBase::fixParameter(...)","Parameter " + std::to_string(nPar) + " already fixed.");
				return false;
			} else {
				return true;
			}
		}
		if (fixed.first > nPar and !appended) {
			newFixedPars.push_back(std::pair<size_t, double>(nPar, value));
			appended = true;
		}
		newFixedPars.push_back(fixed);
	}
	if (!appended)  {
		newFixedPars.push_back(std::pair<size_t, double>(nPar, value));
	}
	_fixedParameters = newFixedPars;
	return true;
}

bool logLikelihoodBase::addCopyParameter(size_t nPar, size_t copyFrom, int nScale, double scaler) {
	if (nScale < -1) {
		utils::makeError("logLikelihoodBase::addCopyParameter(...)","Invalid value for nScale.");
		return false;
	}
	if (scaler == 0.) {
		utils::makeError("logLikelihoodBase::addCopyParameter(...)","Encountered a scaler of zero (fix parameter to zero instead).");
		return false;
	}
	const size_t nParTot = getNparTot();
	if (nPar >= nParTot){
		utils::makeError("logLikelihoodBase::addCopyParameter(...)","nPar > nParTot.");
		return false;
	}
	if (copyFrom >= nParTot) {
		utils::makeError("logLikelihoodBase::addCopyParameter(...)","nPar > nParTot.");
		return false;		
	}
	for (std::tuple<size_t, size_t, int, double> copy : _copyParameters) {
		if (std::get<0>(copy) == copyFrom) {
			utils::makeError("logLikelihoodBase::addCopyParameter(...)","Cannot copy from a copied parameter.");
			return false;
		}
	}
	if (!fixParameter(nPar, std::numeric_limits<double>::quiet_NaN())) { // Set to NaN, since it will have to be overwritten
		utils::makeError("logLikelihoodBase::addCopyParameter(...)","Could not fix the parameter to NaN.");
//		return false; Comment this, since two times fixing should do it
	}
	if (nScale > -1) {
		if (nScale > (int)_nScale) {
			utils::makeError("logLikelihoodBase::addCopyParameter(...)","nScale skipped (Would introduce an unused scale parameter) (nScale = " 
			                                                            + std::to_string(nScale) + ", _nScale = " + std::to_string(_nScale) + ").");
			return false;
		}
		if ((size_t)nScale == _nScale) {
			++_nScale; // Add additional scale parameter
		}
		
	}
	_copyParameters.push_back(std::tuple<size_t, size_t, int, double>(nPar, copyFrom, nScale, scaler));
	return true;
}

realMatrix logLikelihoodBase::DDconstrainedProdAmps(const realVector& params) const {
	complexVector prodAmps = fullParamsToProdAmps(getFullParameters(params));

	utils::makeWarning("logLikelihoodBase::DDconstrainedProdAmps(...)","TODO: Check me!!!");

	size_t nA = 2*prodAmps.size();
	size_t nP = params.size();

	std::vector<std::vector<double > > hessian  = DDeval(prodAmps);
	hessian = makeFinalHessian(params, hessian);
	realMatrix jacobian  = DprodAmpsNumerical(params, _numDelta);
	realMatrix intermed(nP, realVector(nA, 0.));

	for (size_t p = 0; p < nP; ++p) {
		for (size_t q = 0; q < nP; ++q) {
			for (size_t a = 0; a < nA; ++a) {
				intermed[p][a] += jacobian[q][a]*hessian[p][q];
			}
		}
	}
	realMatrix retVal(nA, realVector(nA,0.));
	for (size_t p = 0; p < nP; ++p) {
		for (size_t a = 0; a < nA; ++a) {
			for (size_t b = 0; b < nA; ++b) {
				retVal[a][b] += intermed[p][a] * jacobian[p][b];
			}
		}
	}
	return retVal;
}

realMatrix logLikelihoodBase::DprodAmpsNumerical(const realVector& params, double delta) const {

	utils::makeWarning("logLikelihoodBase::DprodAmpsNumerical(...)","TODO: Check me!!!");

	realVector pars = params; // To be able to change them
	complexVector centerProdAmps = fullParamsToProdAmps(getFullParameters(params));
	realMatrix retVal(params.size(), realVector(2*centerProdAmps.size(), 0.));
	for (size_t p = 0; p < pars.size(); ++p) {
		double dd = pars[p] * delta;
		if (dd == 0.) {
			dd = delta;
		}
		pars[p] += dd;
		complexVector evalProdAmps = fullParamsToProdAmps(getFullParameters(pars));
		pars[p] -= dd;
		for (size_t a = 0; a < 	centerProdAmps.size(); ++a) {
			complex diff = (evalProdAmps[a] - centerProdAmps[a])/dd;
			retVal[p][2*a  ] = diff.real();
			retVal[p][2*a+1] = diff.imag();
		}
	}
	return retVal;
}

complexVector logLikelihoodBase::fullParamsToProdAmps(const realVector& params) const {
	if (params.size() != getNparTot()) {
		utils::makeError("logLikelihoodBase::fullParamsToProdAmps(...)","Parameter size mismatch: " + std::to_string(params.size()) + " should be " + std::to_string(getNparTot()) + ".");
		throw;
	}

	complexVector prodAmps(_nAmpl);
	for(size_t a = 0; a < _nAmpl; ++a) {
		prodAmps[a] = complex(params[2*a], params[2*a+1]);
	}
	return prodAmps;
}


realVector logLikelihoodBase::prodAmpsToFullParams(const complexVector& prodAmps) const {
	if (prodAmps.size() != _nAmpl) {
		utils::makeError("logLikelihoodBase::prodAmpsToFullParams(...)","Parameter size mismatch: " + std::to_string(prodAmps.size()) + " should be " + std::to_string(_nAmpl) + ".");
		throw;
	}
	realVector params(2*_nAmpl);
	for (size_t a = 0; a < _nAmpl; ++a) {
		params[2*a  ] = prodAmps[a].real();
		params[2*a+1] = prodAmps[a].imag();
	}
	return params;
}

logLikelihood::logLikelihood(std::vector<std::shared_ptr<amplitude> > amplitudes, std::shared_ptr<integrator> integral) :
	logLikelihoodBase(amplitudes.size(), integral->kinSignature(), integral), _nSect(1), _maximumIncoherentIntensity(0.), _amplitudes(amplitudes), _points(), _amplitudeCoherenceBorders(0), _contributingWaves(0)
        {
	if (_amplitudes.size() == 0) {
		utils::makeError("logLikelihood::logLikelihood(...)","No amplitudes given.");
		throw;
	}
	_kinSignature = _amplitudes[0]->kinSignature();
	if (not (*_kinSignature == *(integral->kinSignature()))) {
		utils::makeError("logLikelihood::logLikelihood(...)","Kinematic signature of integral differs.");
		throw;
	}
	for (std::shared_ptr<amplitude> a : _amplitudes) {
		if (not(*(a->kinSignature()) == *_kinSignature)) {
			utils::makeError("logLikelihood::logLikelihood(...)","Amplitudes have different kinematic signatures.");
			throw;
		}
	}
	for (size_t a = 0; a < _nAmpl; ++a) {
		std::pair<bool, std::string> waveName =  _integral->getWaveName(a);
		if (not waveName.first) {
			utils::makeError("logLikelihood::logLikelihood(...)","Could not get wave name from integral.");
			throw;
		}
		if (not (_amplitudes[a]->name() == waveName.second)) {
			utils::makeError("logLikelihood::logLikelihood(...)","Wave names do not match.");
			throw;
		}
	}
	if (not setCoherenceBorders(_integral->getCoherenceBorders())) {
		utils::makeError("logLikelihood::logLikelihood(...)","Could not set coherence borders.");
		throw;
	}
}

double logLikelihood::eval(const complexVector& prodAmps) const {
//	std::ofstream outFile;
//	std::string outFileName = "was_geschicht_"+std::to_string(_nAmpl)+".deleteMe";
//	outFile.open(outFileName.c_str());

	double apa = 0.;
	for (complex PA : prodAmps) {
		apa += norm(PA);
	}
	if (prodAmps.size() != _nAmpl) {
		utils::makeError("logLikelihood::eval(...)","Number of production amplitudes does not match: " +std::to_string(prodAmps.size()) + " vs. " + std::to_string(_nAmpl) + ".");
		throw; // Throw here, since a ll = 0. could confuse the program
	}
	double ll = 0.;
#pragma omp parallel for reduction(+:ll)
	for (size_t p = 0; p < _nPoints; ++p) {
		size_t upperSectorBorder = _amplitudeCoherenceBorders[1]; // {0, b1, b2, ..., _nAmpl}.size() = nSector +1
		double intens = 0.;
		complex ampl (0.,0.);
		size_t aCount = 0;
		for (size_t a : _contributingWaves[p]) {
			if (a >= upperSectorBorder) {
				intens += norm(ampl);
				ampl = complex(0.,0.);
				upperSectorBorder = _amplitudeCoherenceBorders[getSector(a)+1];
			}
			ampl += _points[p][aCount] * prodAmps[a];

			++aCount;
		}
		intens += norm(ampl); // Do also for the last sector
//		if (intens == 0.) {
//			continue;
//		}
		ll += log(intens);
	}
//double imll = ll;
//std::cout << " i n t e r m e d i a r y " << imll << std::endl;

	if (_extended) {
		ll -= _integral->totalIntensity(prodAmps, true);
	} else {
		ll -= log(_integral->totalIntensity(prodAmps, true))*(double)_nPoints;
	}
//std::cout << " d e l t a l l " << ll - imll << std::endl;
	if (_maximumIncoherentIntensity >0. ){
		double incoherentIntensity = 0.;
		for (complex pa : prodAmps) {
			incoherentIntensity += norm(pa);
		}
		if (incoherentIntensity > _maximumIncoherentIntensity) {
//			std::cout << " es ist soweit " << incoherentIntensity << std::endl;
			ll -= pow(incoherentIntensity - _maximumIncoherentIntensity,2);
		}
	}
	++_nCalls;
	if (_nCalls%_nCallsPrint == 0) {
		utils::makeInfo("logLikelihood::eval(...)","Call #" + std::to_string(_nCalls) + "; likelihood = " + std::to_string(-ll) + ". Total intensity is: " + std::to_string(_integral->totalIntensity(prodAmps, true)) + ".");
	}
	return -ll;
}

realVector logLikelihood::Deval(const complexVector& prodAmps) const {
	if (prodAmps.size() != _nAmpl) {
		utils::makeError("logLikelihood::Deval(...)","Number of production amplitudes does not match.");;
		throw; // Throw here, since a ll = 0. could confuse the program
	}

	const size_t maxNthread = omp_get_max_threads();
	realMatrix retVal_thread(maxNthread,realVector(2*_nAmpl, 0.));
#pragma omp parallel for
	for (size_t p = 0; p < _nPoints; ++p) {
		const size_t ID = omp_get_thread_num();
		size_t sector = 0;
		size_t upperSectorBorder = _amplitudeCoherenceBorders[1];
		double intens = 0.;
		complexVector sectorAmpls(_nSect, complex(0.,0.));
		size_t aCount = 0;
		for (size_t a : _contributingWaves[p]) {
			if (a >= upperSectorBorder) {
				intens += norm(sectorAmpls[sector]);
				sectorAmpls[sector]= std::conj(sectorAmpls[sector]);
				sector = getSector(a);
				upperSectorBorder = _amplitudeCoherenceBorders[sector+1];
			}
			sectorAmpls[sector] += _points[p][aCount] * prodAmps[a];

			++aCount;
		}
		intens += norm(sectorAmpls[sector]); // Do also for the last sector
		sectorAmpls[sector]= std::conj(sectorAmpls[sector]);
		sector = 0;
		upperSectorBorder = _amplitudeCoherenceBorders[1];
		aCount = 0;
		for (size_t a : _contributingWaves[p]) {
			if (a >= upperSectorBorder) {
				sector = getSector(a);
				upperSectorBorder = _amplitudeCoherenceBorders[sector+1];
			}
			complex factor = _points[p][aCount] * sectorAmpls[sector]/intens * 2.;
			retVal_thread[ID][2*a  ] -= factor.real(); // One factor of -1 since, the NEGATIVE likelihood is used 
			retVal_thread[ID][2*a+1] += factor.imag(); // Two factors of -1: Complex i*i = -1 and since the NEGATIVE likelihood is used

			++aCount;
		}
	}
	realVector retVal(2*_nAmpl, 0.);
	for(size_t i = 0; i < maxNthread; ++i) {
		for (size_t p = 0; p < 2*_nAmpl; ++p) {
			retVal[p] += retVal_thread[i][p];
		}
	}
	realVector Dintegral = _integral->DtotalIntensity(prodAmps, true);
	
	if (_extended) {
		for (size_t a = 0; a < 2*_nAmpl; ++a) {
			retVal[a] += Dintegral[a];
		}
	} else {
		double totalIntens = _integral->totalIntensity(prodAmps, true);
		for (size_t a = 0; a< 2*_nAmpl; ++a) {
			retVal[a] += Dintegral[a]/totalIntens*(double)_nPoints;
		}
	}
	if (_maximumIncoherentIntensity > 0. ){
		double incoherentIntensity = 0.;
		for (complex pa : prodAmps) {
			incoherentIntensity += norm(pa);
		}
		if (incoherentIntensity > _maximumIncoherentIntensity) {
			double fakk = 4*(incoherentIntensity - _maximumIncoherentIntensity);
			for (size_t a = 0; a < _nAmpl; ++a) {
				retVal[2*a  ] += fakk*prodAmps[a].real();
				retVal[2*a+1] += fakk*prodAmps[a].imag();
			}
		}
	}
	return retVal;
}

realMatrix logLikelihood::DDeval(const complexVector& prodAmps) const {
	if (prodAmps.size() != _nAmpl) {
		utils::makeError("logLikelihood::DDeval(...)","Number of production amplitudes does not match.");
		throw; // Throw here, since a ll = 0. could confuse the program
	}
	const size_t maxNthread = omp_get_max_threads();
	std::vector<realMatrix > retVal_thread(maxNthread, realMatrix(2*_nAmpl, realVector (2*_nAmpl, 0.)));
#pragma omp parallel for
	for (size_t p = 0; p < _nPoints; ++p) {
		const size_t ID = omp_get_thread_num();
		size_t sector = 0;
		size_t upperSectorBorder = _amplitudeCoherenceBorders[1];
		double intens = 0.;
		complexVector sectorAmpls(_nSect, complex(0.,0.));
		size_t aCount = 0;
		for (size_t a : _contributingWaves[p]) {
			if (a >= upperSectorBorder) {
				intens += norm(sectorAmpls[sector]);
				sectorAmpls[sector]= std::conj(sectorAmpls[sector]);
				sector = getSector(a);
				upperSectorBorder = _amplitudeCoherenceBorders[sector+1];
			}
			sectorAmpls[sector] += _points[p][aCount] * prodAmps[a];

			++aCount;
		}
		intens += norm(sectorAmpls[sector]); // Do also for the last sector
		sectorAmpls[sector]= std::conj(sectorAmpls[sector]);

		size_t sector_i = 0;
		size_t upperSectorBorder_i = _amplitudeCoherenceBorders[1];
		size_t aiCount = 0;
		for(size_t ai : _contributingWaves[p]) {
			if (ai >= upperSectorBorder_i) {
				sector_i = getSector(ai);
				upperSectorBorder_i = _amplitudeCoherenceBorders[sector_i+1];
			}
			complex factori = _points[p][aiCount] * sectorAmpls[sector_i]/intens * 2.;

			size_t sector_j = 0;
			size_t upperSectorBorder_j = _amplitudeCoherenceBorders[1];
			size_t ajCount = 0;
			for (size_t aj : _contributingWaves[p]) {
				if (aj >= upperSectorBorder_j) {
					sector_j = getSector(aj);
					upperSectorBorder_j = _amplitudeCoherenceBorders[sector_j+1];
				}
				if (sector_i == sector_j) {
					complex factor = _points[p][aiCount] * std::conj(_points[p][ajCount])/intens*2.;
					retVal_thread[ID][2*ai  ][2*aj  ] -= factor.real();
					retVal_thread[ID][2*ai  ][2*aj+1] -= factor.imag();
					retVal_thread[ID][2*ai+1][2*aj  ] += factor.imag();
					retVal_thread[ID][2*ai+1][2*aj+1] -= factor.real();
				}
				complex factorj = -_points[p][ajCount] * sectorAmpls[sector_j]/intens *2.; // -1/intens and then the same factor;

				retVal_thread[ID][2*ai  ][2*aj  ] -= factori.real() * factorj.real();
				retVal_thread[ID][2*ai  ][2*aj+1] += factori.real() * factorj.imag();
				retVal_thread[ID][2*ai+1][2*aj  ] += factori.imag() * factorj.real();
				retVal_thread[ID][2*ai+1][2*aj+1] -= factori.imag() * factorj.imag();

				++ajCount;
			}
			++aiCount;
		}
	}
	realMatrix retVal(2*_nAmpl, realVector (2*_nAmpl, 0.));
	for (size_t i = 0; i < maxNthread; ++i) {
		for (size_t p = 0; p < 2*_nAmpl; ++p) {
			for (size_t q = 0; q < 2*_nAmpl;++q) {
				retVal[p][q] += retVal_thread[i][p][q];
			}
		}
	}

	realMatrix DDintegral = _integral->DDtotalIntensity(prodAmps, true);
	if (_extended) {
		for (size_t i = 0; i < 2*_nAmpl; ++i) {
			for (size_t j = 0; j < 2*_nAmpl; ++j) {
				retVal[i][j] += DDintegral[i][j];
			}
		}
	} else {
		double totalIntens = _integral->totalIntensity(prodAmps, true);
		realVector Dintegral = _integral->DtotalIntensity(prodAmps, true);
		for (size_t i = 0; i < 2*_nAmpl; ++i) {
			for (size_t j = 0; j < 2*_nAmpl; ++j) {
				retVal[i][j] += (DDintegral[i][j]/totalIntens - Dintegral[i]*Dintegral[j]/totalIntens/totalIntens)*(double)_nPoints;
			}
		}
	}
	if (_maximumIncoherentIntensity > 0. ){
		double incoherentIntensity = 0.;
		for (complex pa : prodAmps) {
			incoherentIntensity += norm(pa);
		}
		if (incoherentIntensity > _maximumIncoherentIntensity) {
			double fakk = 4*(incoherentIntensity - _maximumIncoherentIntensity);
			for (size_t a = 0; a < _nAmpl; ++a) {
				retVal[2*a  ][2*a  ] += fakk;
				retVal[2*a+1][2*a+1] += fakk;
				for (size_t b = 0; b < _nAmpl; ++b) {
					retVal[2*a  ][2*b  ] += 8*prodAmps[a].real()*prodAmps[b].real();
					retVal[2*a+1][2*b  ] += 8*prodAmps[a].imag()*prodAmps[b].real();
					retVal[2*a  ][2*b+1] += 8*prodAmps[a].real()*prodAmps[b].imag();
					retVal[2*a+1][2*b+1] += 8*prodAmps[a].imag()*prodAmps[b].imag();
				}
			}
		}
	}
	return retVal;
}

bool logLikelihood::loadDataPoints(const realMatrix& dataPoints, size_t maxNperEvent) {
	if (dataPoints.size() == 0) {
		utils::makeError("logLikelihood::loadDataPoints(...)","Not data points given.");
		return false;
	}
	_nPoints = dataPoints.size();
	_points  = complexMatrix (_nPoints, complexVector (maxNperEvent));
	_contributingWaves= std::vector<sizeVector >(_nPoints, sizeVector(maxNperEvent));
	sizeVector counts(_nPoints, 0);
	bool resized = false;
	for(size_t a = 0; a < _nAmpl; ++a) {
		std::pair<bool, complex > diag = _integral->element(a,a, false);
		if (not diag.first) {
			utils::makeError("logLikelihood::loadDataPoints(...)","Could not get diagonal integral.");
			return false;
		}
		double norm = 0.;
		if (diag.second.real() != 0.) {
			norm = 1./pow(diag.second.real(), .5);
		}
		for (size_t p = 0; p < _nPoints; ++p) {
			complex ampl = _amplitudes[a]->eval(dataPoints[p]);
			if (std::isnan(ampl.real()) || std::isnan(ampl.imag())) {
				utils::makeError("logLikelihood::loadDataPoints(...)","NaN amplitude encountered for wave '" + _amplitudes[a]->name() + "': " + utils::to_string(ampl) + ".");
				utils::makeError("logLikelihood::loadDataPoints(...)","Kinematics of data point " + std::to_string(p) + " are: ");
				for (const double& point: dataPoints[p]) {
					std::cout << point << " ";
				}
				std::cout << std::endl;
				return false;
			}

			if (norm == 0. && ampl != 0.) {
				utils::makeError("logLikelihood::loadDataPoints(...)","Zero-norm wave has non-vanishing amplitude.");
				return false;
			}
			ampl *= norm;
			if (ampl != complex(0.,0.)) {
				if (counts[p] == maxNperEvent) {
					utils::makeWarning("logLikelihood::loadDataPoints(...)","Found more non-zero amplitudes for an event, than set by maxNperEvent = " 
					                                                        + std::to_string(maxNperEvent) + " (Set this value higher to avoid dynamic resizing).");
					maxNperEvent *= 2;
					for (size_t p2 = 0; p2 < _nPoints; ++p2) {
						_contributingWaves[p2].resize(maxNperEvent);
						_points[p2].resize(maxNperEvent);
					}
					resized = true;
				}
				_contributingWaves[p][counts[p]] = a;
				_points[p][counts[p]] = ampl;
				++counts[p];
			}
		}
	}
	if (resized) {
		size_t maxSize = 0;
		for (size_t size : counts) {
			maxSize = std::max(maxSize, size);
		}
		utils::makeInfo("logLikelihood::loadDataPoints(...)","Found a maximum size of " + std::to_string(maxSize) + ".");
	}
	for (size_t p = 0; p < _nPoints; ++p) {
		_contributingWaves[p].resize(counts[p]);
		_points[p].resize(counts[p]);
	}
	return true;
}

bool logLikelihood::setCoherenceBorders(sizeVector borders) {
	if (borders[0] != 0) {
		utils::makeError("logLikelihood::setCoherenceBorders(...)","First border has to be zero.");
		return false;
	}
	_amplitudeCoherenceBorders = sizeVector(borders.size(),0);
	for (size_t i = 1; i < borders.size(); ++i) {
		if (borders[i] <= _amplitudeCoherenceBorders[i-1]) {
			utils::makeError("logLikelihood::setCoherenceBorders(...)","Borders are nor ordered.");
			return false;
		}
		 _amplitudeCoherenceBorders[i] = borders[i];
	}
	if (_amplitudeCoherenceBorders[_amplitudeCoherenceBorders.size()-1] != _nAmpl) {
		utils::makeError("logLikelihood::setCoherenceBorders(...)","Last border has to be _nAmpl ( = " + std::to_string(_nAmpl) + " != " 
		                                                           +std::to_string(_amplitudeCoherenceBorders[_amplitudeCoherenceBorders.size()-1]) + ").");
		return false;
	}
	_nSect = _amplitudeCoherenceBorders.size()-1;
	return true;
}

size_t logLikelihood::getSector(size_t a) const {
	for (size_t s = 0; s < _nSect; ++s) {
		if (_amplitudeCoherenceBorders[s] <= a && _amplitudeCoherenceBorders[s+1] > a) {
			return s;
		}
	}
	utils::makeError("logLikelihood::getSector(...)","No sector found for amplitude index " + std::to_string(a) + " (_nAmpl = " + std::to_string(_nAmpl) + ").");
	throw;
	return 0;
}

linearLogLikelihood::linearLogLikelihood(std::vector<std::shared_ptr<amplitude> > amplitudes, std::shared_ptr<integrator> integral, realMatrix linearTransformation) :
	logLikelihood(amplitudes, integral), _transformation(linearTransformation) {

	if (_transformation.size() == 0) {
		utils::makeError("linearLogLikelihood::linearLogLikelihood","_transformation has size 0.");
		throw;
	}
	if (_transformation[0].size() == 0) {
		utils::makeError("linearLogLikelihood::linearLogLikelihood","_transformation[0] has size 0.");
		throw;
	}
	_dimFull      = _transformation.size();
	_dimTransform = _transformation[0].size();
	if (_dimFull < _dimTransform) {
		utils::makeWarning("linearLogLikelihood::linearLogLikelihood","Input dimension larger than fit dimension, leading to artificial parameters and linear dependences. (Ignore, if you know what you are doing).");
	}
	for (const realVector& line : _transformation) {
		if (line.size() != _dimTransform) {
			utils::makeError("linearLogLikelihood::linearLogLikelihood","Lines of the transformation are not equal.");
			throw;
		}
	}
// DO NOT CHECK NPAR HERE, SINCE THIS CAN AN WILL CHANGE (E.G. COPY PARAMETERS). / // /// //// ///// ////// /////// ////// ///// //// /// // /
};

double linearLogLikelihood::nloptCall(const realVector& x, realVector& grad) const {
	if (x.size() != getNparFit()) {
		utils::makeError("linearLogLikelihood::nloptCall(...)","Number of parameters does not match (" + std::to_string(x.size()) + " vs. " + std::to_string(getNparFit()) + ").");
	}
	realVector fullParameters = getFullRealParametersLinear(x);
	realVector gradFull;
	if (grad.size() > 0) {
		gradFull = realVector(_dimFull, 0.);
	}
	if (_trackFit) {
		if (_storeParams.size() <= nCalls()) {
			_storeParams.resize(2*nCalls()+1);
		}
		_storeParams[nCalls()] = realVector(x);
	}
	double retVal = logLikelihoodBase::nloptCall(fullParameters, gradFull);
	if (grad.size() > 0) {
		for (size_t i = 0; i < _dimFull; ++i) {
			for (size_t j = 0; j < _dimTransform; ++j) {
				grad[j] += _transformation[i][j] * gradFull[i];
			}
		}
	}
	return retVal;
};

double linearLogLikelihood::linearEval(const realVector& params) const {
	if (params.size() != getNparFit()) {
		utils::makeError("linearLogLikelihood::linearEval(...)","Number of parameters does not match (" + std::to_string(params.size()) + " vs. " + std::to_string(getNparFit()) + ").");
	}
	complexVector prodAmps = makeProdAmpsFromFitPars(params);
	return logLikelihood::eval(prodAmps);
}

realVector linearLogLikelihood::linearDeval(const realVector& params) const {
	if (params.size() != getNparFit()) {
		utils::makeError("linearLogLikelihood::linearDeval(...)","Number of parameters does not match (" + std::to_string(params.size()) + " vs. " + std::to_string(getNparFit()) + ").");
	}
	complexVector prodAmps = makeProdAmpsFromFitPars(params);
	realVector fullRealParams          = getFullRealParametersLinear(params);
	realVector gradFull                = makeFinalGradient(fullRealParams,Deval(prodAmps));
	realVector retVal(getNparFit(), 0.);
	
	for (size_t i = 0; i < _dimFull; ++i) {
		for (size_t j = 0; j < _dimTransform; ++j) {
			retVal[j] += _transformation[i][j] * gradFull[i];
		}
	}

	return retVal;	
}

realMatrix linearLogLikelihood::linearDDeval(const realVector& params) const {
	if (params.size() != getNparFit()) {
		utils::makeError("linearLogLikelihood::linearDDeval(...)","Number of parameters does not match (" + std::to_string(params.size()) + " vs. " + std::to_string(getNparFit()) + ").");
	}
	complexVector prodAmps = makeProdAmpsFromFitPars(params);
	realVector fullRealParams          = getFullRealParametersLinear(params);
	realMatrix hessFull  = makeFinalHessian(fullRealParams, logLikelihood::DDeval(prodAmps));
	realMatrix intermedHessian(_dimFull, realVector(_dimTransform,0.));
	
	for (size_t i = 0; i < _dimFull; ++i) {
		for (size_t j = 0; j < _dimFull; ++j) {
			for (size_t k = 0; k < _dimTransform; ++k) {
				intermedHessian[i][k] += _transformation[j][k] * hessFull[i][j];
			}
		}
	}
	realMatrix  retVal(_dimTransform, realVector(_dimTransform,0.));
	for (size_t i = 0; i < _dimTransform; ++i) {
		for (size_t j = 0; j < _dimFull; ++j) {
			for (size_t k = 0; k < _dimTransform; ++k) {
				retVal[k][i] += _transformation[j][k] * intermedHessian[j][i];
			}
		}
	}

	return retVal;	
}

realVector linearLogLikelihood::getFullRealParametersLinear(const realVector& fitPars) const {
	if (fitPars.size() != _dimTransform) {
		utils::makeError("linearLogLikelihood::getFullRealParametersLinear(...)","Number of fit parameters (" + std::to_string(fitPars.size()) + ") does not match the matrix dimension (" 
		                                                                     + std::to_string(_dimTransform) + ").");
		throw;
	}
	realVector fullRealParameters(_dimFull, 0.);
	for (size_t i = 0; i < _dimFull; ++i) {
		for (size_t j = 0; j < _dimTransform; ++j) {
			fullRealParameters[i] += _transformation[i][j] * fitPars[j];
		}
	}
	return fullRealParameters;
}

complexVector linearLogLikelihood::makeProdAmpsFromFitPars(const realVector& fitPars) const {
	return fullParamsToProdAmps(getFullParameters(getFullRealParametersLinear(fitPars)));
}

logLikelihood_withPrior::logLikelihood_withPrior(std::vector<std::shared_ptr<amplitude> > amplitudes, std::shared_ptr<integrator> integral) :
	logLikelihood(amplitudes, integral), _noNormWarn(false), _interferencePriorStrength(0.), _nPrior(0), _priorStrengths(), _priorDirections() {}

double logLikelihood_withPrior::eval(const complexVector& prodAmps) const {
	double retVal = logLikelihood::eval(prodAmps);
	for (size_t p = 0; p < _nPrior; ++p) {
		complex scalarProduct(0.,0.);
		for (size_t a = 0; a < _nAmpl; ++a) {
			scalarProduct += prodAmps[a] * _priorDirections[p][a];
		}
		retVal += _priorStrengths[p] * std::norm(scalarProduct);
	}
	if (_interferencePriorStrength != 0.) {
		double coherent = _integral->totalIntensity(prodAmps, false);
		double incoherent = 0.;
		for (const complex& a : prodAmps) {
			incoherent += std::norm(a);
		}
		retVal += interferencePriorFunc(coherent, incoherent);
	}
	if (_nCalls%_nCallsPrint == 0) {
		utils::makeInfo("logLikelihood_withPrior::eval(...)", "Call #" + std::to_string(_nCalls) + " likeWithPrior: " + std::to_string(retVal) + ".");
	}
	return retVal;
}

realVector logLikelihood_withPrior::Deval(const complexVector& prodAmps) const {
	realVector retVal = logLikelihood::Deval(prodAmps);
	for (size_t p = 0; p < _nPrior; ++p) {
		complex scalarProduct(0.,0.);
		for (size_t a = 0; a < _nAmpl; ++a) {
			scalarProduct += prodAmps[a] * _priorDirections[p][a];
		}
		for (size_t a = 0; a < _nAmpl; ++a) {
			complex deriv = 2.*scalarProduct*std::conj(_priorDirections[p][a]);
			retVal[2*a  ] += _priorStrengths[p] * deriv.real();
			retVal[2*a+1] += _priorStrengths[p] * deriv.imag();
		}
	}
        if (_interferencePriorStrength != 0.) {
                double coherent = _integral->totalIntensity(prodAmps, false);
                double incoherent = 0.;
                for (const complex& a : prodAmps) {
                        incoherent += std::norm(a);
                }
                std::pair<double,double> Dipf = DinterferencePriorFunc(coherent, incoherent);
		realVector Dcoherent = _integral->DtotalIntensity(prodAmps, false);
		for (size_t a = 0; a < _nAmpl; ++a) {
			retVal[2*a  ] += Dipf.first * Dcoherent[2*a  ] + Dipf.second * 2. * prodAmps[a].real();
			retVal[2*a+1] += Dipf.first * Dcoherent[2*a+1] + Dipf.second * 2. * prodAmps[a].imag();
		}
        }
	return retVal;
}

realMatrix logLikelihood_withPrior::DDeval(const complexVector& prodAmps) const {
	realMatrix retVal = logLikelihood::DDeval(prodAmps);
	for (size_t p = 0; p < _nPrior; ++p) {
		for (size_t a_i = 0; a_i < _nAmpl; ++a_i) {
			for (size_t a_j = 0; a_j < _nAmpl; ++a_j) {
				complex product = 2.*_priorDirections[p][a_i] * std::conj(_priorDirections[p][a_j]);
				retVal[2*a_i  ][2*a_j  ] += _priorStrengths[p] * product.real();
				retVal[2*a_i  ][2*a_j+1] += _priorStrengths[p] *-product.imag();
				retVal[2*a_i+1][2*a_j  ] += _priorStrengths[p] * product.imag();
				retVal[2*a_i+1][2*a_j+1] += _priorStrengths[p] * product.real();
			}
		}
	}
        if (_interferencePriorStrength != 0.) {

                double coherent = _integral->totalIntensity(prodAmps, false);
                double incoherent = 0.;
                for (const complex& a : prodAmps) {
                        incoherent += std::norm(a);
                }
                std::pair<double,double> Dipf = DinterferencePriorFunc(coherent, incoherent);
		Dipf.first *= 1.;
		realVector DDipf     = DDinterferencePriorFunc(coherent, incoherent);
                realVector Dcoherent = _integral->DtotalIntensity(prodAmps, false);
		realMatrix DDcoherent = _integral->DDtotalIntensity(prodAmps, false);

                for (size_t a = 0; a < _nAmpl; ++a) {
			retVal[2*a  ][2*a  ] += 2.*Dipf.second;
			retVal[2*a+1][2*a+1] += 2.*Dipf.second;

			for (size_t b = 0; b < _nAmpl; ++b) {

				retVal[2*a  ][2*b  ] += Dipf.first * DDcoherent[2*a  ][2*b  ];	
				retVal[2*a+1][2*b  ] += Dipf.first * DDcoherent[2*a+1][2*b  ];	
				retVal[2*a  ][2*b+1] += Dipf.first * DDcoherent[2*a  ][2*b+1];	
				retVal[2*a+1][2*b+1] += Dipf.first * DDcoherent[2*a+1][2*b+1];	

				retVal[2*a  ][2*b  ] += DDipf[0] * Dcoherent[2*a  ] * Dcoherent[2*b  ];
				retVal[2*a+1][2*b  ] += DDipf[0] * Dcoherent[2*a+1] * Dcoherent[2*b  ];
				retVal[2*a  ][2*b+1] += DDipf[0] * Dcoherent[2*a  ] * Dcoherent[2*b+1];
				retVal[2*a+1][2*b+1] += DDipf[0] * Dcoherent[2*a+1] * Dcoherent[2*b+1];

				retVal[2*a  ][2*b  ] += DDipf[1] * Dcoherent[2*a  ] * 2. * prodAmps[b].real();
				retVal[2*a+1][2*b  ] += DDipf[1] * Dcoherent[2*a+1] * 2. * prodAmps[b].real();
				retVal[2*a  ][2*b+1] += DDipf[1] * Dcoherent[2*a  ] * 2. * prodAmps[b].imag();
				retVal[2*a+1][2*b+1] += DDipf[1] * Dcoherent[2*a+1] * 2. * prodAmps[b].imag();

				retVal[2*a  ][2*b  ] += DDipf[1] * 2. * prodAmps[a].real() * Dcoherent[2*b  ];
				retVal[2*a+1][2*b  ] += DDipf[1] * 2. * prodAmps[a].imag() * Dcoherent[2*b  ];
				retVal[2*a  ][2*b+1] += DDipf[1] * 2. * prodAmps[a].real() * Dcoherent[2*b+1];
				retVal[2*a+1][2*b+1] += DDipf[1] * 2. * prodAmps[a].imag() * Dcoherent[2*b+1];

				retVal[2*a  ][2*b  ] += DDipf[2] * 4. * prodAmps[a].real() * prodAmps[b].real();
				retVal[2*a+1][2*b  ] += DDipf[2] * 4. * prodAmps[a].imag() * prodAmps[b].real();
				retVal[2*a  ][2*b+1] += DDipf[2] * 4. * prodAmps[a].real() * prodAmps[b].imag();
				retVal[2*a+1][2*b+1] += DDipf[2] * 4. * prodAmps[a].imag() * prodAmps[b].imag();

			}

                }
        }

	return retVal;
}

bool logLikelihood_withPrior::setInterferencePriorStrength(double strength) {
	if (strength < 0.) {
		utils::makeWarning("logLikelihood_withPrior::setInterferencePriorStrenght(...)","Strength < 0. leads to instable likelihood.");
	}
	_interferencePriorStrength = strength;
	return true;
}

double logLikelihood_withPrior::interferencePriorFunc(double coherent, double incoherent) const {
	double retVal = _interferencePriorStrength *  pow(coherent/incoherent - 1., 2);
	return retVal;
}

std::pair<double, double> logLikelihood_withPrior::DinterferencePriorFunc(double coherent, double incoherent) const {
	double fakk = 2*_interferencePriorStrength*(coherent/incoherent-1.);
	double Dcoh = fakk/incoherent;
	double Dinc = -fakk*coherent/pow(incoherent,2);
	return std::pair<double,double>(Dcoh, Dinc);
}

realVector logLikelihood_withPrior::DDinterferencePriorFunc(double coherent, double incoherent) const {
	realVector retVal(3,0.);
	retVal[0] = 2*_interferencePriorStrength/pow(incoherent,2); // DDcoherent
	retVal[1] = 2*_interferencePriorStrength*(1./pow(incoherent,2) - 2.*coherent/pow(incoherent,3)); //DcoherentDincoherent
	retVal[2] = 2*_interferencePriorStrength*(3.*pow(coherent,2)/pow(incoherent,4) - 2.*coherent/pow(incoherent,3)); //DDincoherent
	return retVal;
}

bool logLikelihood_withPrior::addPriorDirection(double strength, const complexVector direction) {
	if (strength < 0.) {
		utils::makeError("logLikelihood_withPrior::addPriorDirection(...)","Prior strength should not be below zero: " + std::to_string(strength) + ".");
		return false;
	}
	if (direction.size() != _nAmpl) {
		utils::makeError("logLikelihood_withPrior::addPriorDirection(...)","Prior direction has the wrong size: "+std::to_string(direction.size())+" (should be "+std::to_string(_nAmpl)+").");
		return false;
	}
	double norm = 0.;
	for (const complex d : direction) {
		norm += std::norm(d);
	}
	if (norm < .99 || norm > 1.01) {
		if (!_noNormWarn) {
			utils::makeWarning("logLikelihood_withPrior::addPriorDirection(...)","Prior direction is not normalized. Rescale. Use the 'stregth' parameter to control the prior strength.");
			_noNormWarn = true;
		}
	}
	_priorDirections.push_back(direction);
	const double scale = pow(std::norm(norm), .5);
	for (size_t a = 0; a < _nAmpl; ++a) {
		_priorDirections[_nPrior][a] /= scale;
	}
	_priorStrengths.push_back(strength);
	++_nPrior;
	return true;	
}

logLikelihoodAllFree::logLikelihoodAllFree (realVector binning, std::vector<std::shared_ptr<angularDependence> > freedAmplitudes, std::shared_ptr<integrator> integral) :
	logLikelihoodBase ((binning.size()-1)*freedAmplitudes.size(), integral->kinSignature(), integral), 
	_nBins(binning.size()-1), 
	_binning(binning), 
	_eventsPerBin(std::vector<sizeVector >(binning.size()-1, sizeVector(binning.size() - 1, 0))),
	_amplitudesInBin(std::vector<complexMatrix >(binning.size()-1, 
	                 complexMatrix(binning.size()-1, 
	                 complexVector(freedAmplitudes.size(), 
	                 complex(0.,0.))))) {
	for (size_t b = 0; b < _nBins; ++b) {
		if (_binning[b+1] <= _binning[b]) {
			utils::makeError("logLikelihoodAllFree::logLikelihoodAllFree(...)","Binning not ordered.");
			throw;
		}
	}
	if (_kinSignature->nKin() != 3) {
		utils::makeError("logLikelihoodAllFree::logLikelihoodAllFree(...)","Number of kinemtaic variables is not three. Unknown process. Abortring!");
		throw;
	}
}

double logLikelihoodAllFree::eval(const complexVector& prodAmps) const {
	utils::makeError("logLikelihoodAllFree::eval(...)","NOT IMPLEMENTED (prodAmps.size() = " + std::to_string(prodAmps.size()) + ").");
	throw;
	return 0;
}

realVector logLikelihoodAllFree::Deval(const complexVector& prodAmps) const {
	utils::makeError("logLikelihoodAllFree::Deval(...)","NOT IMPLEMENTED (prodAmps.size() = " + std::to_string(prodAmps.size()) + ").");
	throw;
	realVector();
}

realMatrix  logLikelihoodAllFree::DDeval(const complexVector& prodAmps) const {
	utils::makeError("logLikelihoodAllFree::DDeval(...)","NOT IMPLEMENTED (prodAmps.size() = " + std::to_string(prodAmps.size()) + ").");
	throw;
	realMatrix();
}

bool logLikelihoodAllFree::loadDataPoints(const realMatrix& dataPoints, size_t maxNperEvent) {
	if (dataPoints.size() == 0) {
		utils::makeError("logLikelihood::loadDataPoints(...)","Not data points given (maxNperEvent = " + std::to_string(maxNperEvent) + ").");
		return false;
	}
	const double s = dataPoints[0][0];
	_nPoints = dataPoints.size();
	std::vector<std::vector<std::pair<double, double> > > averageMasses(_nBins, std::vector<std::pair<double, double> >(_nBins, std::pair<double, double>(0.,0.)));
	for (const realVector& point : dataPoints) {
		std::pair<bool, std::pair<size_t, size_t> > bin = findBin(point);
		if (not bin.first) {
			utils::makeError("logLikelihood::loadDataPoints(...)","Bin not found.");
			return false;
		}
		_eventsPerBin[bin.second.first][bin.second.second]        += 1;
		averageMasses[bin.second.first][bin.second.second].first  += point[1];
		averageMasses[bin.second.first][bin.second.second].second += point[2];
	}
	for (size_t b1 = 0; b1 < _nBins; ++b1) {
		for (size_t b2 = 0; b2 < _nBins; ++b2) {
			size_t nBin = _eventsPerBin[b1][b2];
			if (nBin == 0) {
				continue;
			}
			realVector binPoint = {s, averageMasses[b1][b2].first/nBin, averageMasses[b1][b2].second/nBin};
			utils::makeError("logLikelihood::loadDataPoints(...)","NOT IMPLEMENTED: Not fully implemented.");
			throw;
		}
	}
	return true;
}

std::pair<bool, std::pair<size_t, size_t> > logLikelihoodAllFree::findBin(const realVector& point) const {
	bool found1 = false;
	size_t bin1 = 0;
	bool found2 = false;
	size_t bin2 = 0;
	if (point.size() != _kinSignature->nKin()) {
		utils::makeError("logLikelihoodAllFree::findBin(...)","Number of kinematic variables does not match.");
		return std::pair<bool, std::pair<size_t, size_t> >(false, std::pair<size_t, size_t>(0,0));
	}
	const double s12 = point[1];
	const double s13 = point[2];
	for (size_t b = 0; b < _nBins; ++b) {
		const double sMin = _binning[b];
		const double sMax = _binning[b+1];
		if (sMin < s12 and s12 <= sMax) {
			if (found1) {
				utils::makeError("logLikelihoodAllFree::findBin(...)","Two bins found for s_{12}... Returning false.");;
				return  std::pair<bool, std::pair<size_t, size_t> >(false, std::pair<size_t, size_t>(0,0));
			}
			found1 = true;
			bin1   = b;
		}
		if (sMin < s13 and s13 <= sMax) {
			if (found2) {
				utils::makeError("logLikelihoodAllFree::findBin(...)","Two bins found for s_{13}... Returning false.");
				return  std::pair<bool, std::pair<size_t, size_t> >(false, std::pair<size_t, size_t>(0,0));
			}
			found2 = true;
			bin2   = b;
		}
		if (found1 and found2) { // This breaks the total check of the binnig, i.e. if a second bin would be found, but ordering is checked in the constuctor anyway...
			break;
		}
	}
	return std::pair<bool, std::pair<size_t, size_t> >(found1 and found2, std::pair<size_t, size_t>(bin1, bin2));
}
