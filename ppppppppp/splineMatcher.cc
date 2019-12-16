#include "splineMatcher.h"

splineMatcher::splineMatcher(const std::vector<std::shared_ptr<spline_deg> >& functionsFirstBin){
	_degree = functionsFirstBin.size() - 2;
	_nBin   = 1;
	_transformationMatrix = getM(functionsFirstBin);
	_functions = functionsFirstBin;
}

realMatrix splineMatcher::getM(const std::vector<std::shared_ptr<spline_deg> >& functions) {
	// Creates tand inverts the matrix, that makes a single-bin matching to the values: 
	// {f^(n')(flow), f^(n-1')(low), ..., f"(low), f'(low), f(low), f(high)}
	// M: values -> coeffs
	if (functions.size() == 0) {
		utils::makeError("splineMatcher::getM(...)","No functions received.");
		return realMatrix();
	}
	const size_t degree = functions.size() - 2;
	realMatrix retVal(degree+2, realVector(degree+2, 0.));
	const double mMin = functions[0]->getParameter(0).second;
	const double mMax = functions[0]->getParameter(1).second;
	for (size_t i = 0; i < degree+2; ++i) {
		for (size_t d = 0; d < degree+1; ++d) {
			retVal[i][d] = functions[i]->Deval(mMin, degree-d);
		}
		retVal[i][degree+1] = functions[i]->Deval(mMax,0.);
	}
	return utils::inv(retVal);
}

realMatrix splineMatcher::getLastN() const {
	// Creadtes the values f^(n')(high), f^(n-1')(high), ..., f"(high), f'(high)
	// For a single bin, to match the next one to
	// N coeffs -> values
	const double m = _functions[(_nBin-1)*(_degree+2)]->getParameter(1).second;
	realMatrix N(_degree+2, realVector(_degree,0.));
	for (size_t i = 0 ; i < _degree+2; ++i) {
		for (size_t j = 0; j < _degree; ++j) {
			N[i][j] = _functions[(_nBin-1)*(_degree+2)+i]->Deval(m,_degree-j);
		}
	}
	return N;
}

bool splineMatcher::addBin(const std::vector<std::shared_ptr<spline_deg> >& newFunctions) {
	if (newFunctions[0]->getParameter(0).second != _functions[_functions.size()-1]->getParameter(1).second) {
		utils::makeError("splineMatcher::addBin(...)","Bin borders not contiguous.");
		return false;
	}

	if (newFunctions.size() != _degree+2) {
		utils::makeError("splineMatcher::addBin(...)", "Wrong number of fiunctions: " + std::to_string(newFunctions.size()) + " instead of " + std::to_string(_degree+2) + ".");
		return false;
	}
	realMatrix M = getM(newFunctions); // Make relevent function evlauations
	realMatrix N = getLastN();
	realMatrix trans(_degree+2+_nBin, realVector(_degree+2, 0.)); // Construct transformation matrix from params to values
	for (size_t i = 0; i < _degree+1+_nBin;++i){for(size_t j = 0; j<_degree+2;++j){for(size_t k=0;k<_degree;++k){
		trans[i][k] += _transformationMatrix[i][(_degree+2)*(_nBin-1)+j] * N[j][k];
	}}}
	trans[_degree+  _nBin][_degree  ] = 1.; // Lower bin border value directly taken
	trans[_degree+1+_nBin][_degree+1] = 1.; // Upper bin border value directly taken
	// trans: params -> coeffs
	for (realVector& line : _transformationMatrix) {
		line.resize((_degree+2)*(_nBin+1), 0.); // Create the new entries
	}
	_transformationMatrix.push_back(realVector((_degree+2)*(_nBin+1), 0.)); // Add line for the new parameter
	for (size_t i =0; i < _degree+2; ++i){
		size_t I = i + (_degree+2)*_nBin; // Only the back part of the matrox is changed
		for(size_t j = 0; j < _degree+2; ++j){for(size_t k=0; k<_degree+_nBin+2;++k){
			_transformationMatrix[k][I] +=  trans[k][j] * M[j][i]; // The index ordeign here is the forst place to look at for errors
	}}}
	++_nBin;
	for (std::shared_ptr<spline_deg> f : newFunctions) {
		_functions.push_back(f);
	}
	return true;
}
