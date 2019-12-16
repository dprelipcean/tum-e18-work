#include"integrator.h"
#include"utils.h"
#include"types.h"
#include<iostream>
#include<fstream>
#include<iomanip>
#include<limits>

#include"omp.h"
integrator::integrator(size_t integralPoints, std::shared_ptr<generator> pointGenerator, const std::vector<std::shared_ptr<amplitude> >& amplitudes, std::shared_ptr<efficiencyFunction>& efficiency):
	_isIntegrated(false), _numLim(1.e-15), _nAmpl(amplitudes.size()), _nPoints(integralPoints), _amplitudes(amplitudes), _generator(pointGenerator), _efficiency(efficiency), _integralMatrix(), _accCorrIntegralMatrix(), _realIntegralMatrix(), _realAccCorrMatrix() {
	if (_amplitudes.size() == 0) {
		utils::makeError("integrator::integrator(...)","No amplitudes given.");
		throw;
	}
	_kinSignature = _amplitudes[0]->kinSignature();
	for (std::shared_ptr<amplitude> a : _amplitudes) {
		if (not (*(a->kinSignature()) == *_kinSignature)) {
			utils::makeError("integrator::integrator(...)","Amplitudes with different signatures found.");
			throw;
		}
	}
	if (not (*(_generator->kinSignature()) == *(_amplitudes[0]->kinSignature()))) {
		utils::makeError("integrator::integrator(...)","Kinematic signatures in amplitudes and generator differ.");
		throw;
	}
	if (not (*(_efficiency->kinSignature()) == *(_amplitudes[0]->kinSignature()))) {
		utils::makeError("integrator::integrator(...)","Kinematic signatures in amplitudes and efficiency differ.");
		throw;
	}
	if (!utils::checkComplexDouble()) {
		utils::makeError("integrator::integrator(...)","Complex double has wrong structure.");
		throw;
	}
	_amplitudeCoherenceBorders = {0, _nAmpl};
}

bool integrator::integrate() {
	const size_t maxNthread = omp_get_max_threads();
	std::vector<complexMatrix> integralMatrix_omp 
	    (maxNthread, complexMatrix(nAmpl(), complexVector(nAmpl(), complex(0.,0.))));
	std::vector<complexMatrix> accCorrIntegralMatrix_omp
	    (maxNthread, complexMatrix(nAmpl(), complexVector(nAmpl(), complex(0.,0.))));
#pragma omp parallel for
	for (size_t point = 0; point < _nPoints; ++point) {
		const size_t ID = omp_get_thread_num();
		realVector kin = _generator->generate();
		complexVector ampl(nAmpl(), complex(0.,0.));
		double eff = _efficiency->eval(kin);
		size_t countAmp = 0;
		for (std::shared_ptr<amplitude> a : _amplitudes) {
			ampl[countAmp] = a->eval(kin);
//			if (std::isnan(ampl[countAmp].real()) || std::isnan(ampl[countAmp].imag())) {
//				std::cout << "integrator::loadIntegrals(...): ERROR: NaN amplitude encountered for " << a->name() << std::endl;
//				return false;
//			}
			++countAmp;
		}
		for (size_t i = 0; i < nAmpl(); ++i) {
			if (ampl[i] == complex(0.,0.)) {
				continue;
			}
			for (size_t j = 0; j < nAmpl(); ++j) {
				if (ampl[j] == complex(0.,0.)) {
					continue;
				}
				integralMatrix_omp[ID][i][j] += std::conj(ampl[i])*ampl[j];
				accCorrIntegralMatrix_omp[ID][i][j] += std::conj(ampl[i])*ampl[j]*eff;
			}
		}
		if (point % 1000 == 0) { 
			std::cout << "#" << point << std::endl;
		}
	}
	complexMatrix integralMatrix(nAmpl(), complexVector(nAmpl(), complex(0.,0.)));
	complexMatrix accCorrIntegralMatrix(nAmpl(), complexVector (nAmpl(), complex(0.,0.)));
	for (size_t ID = 0; ID < maxNthread; ++ID) {
		for (size_t i = 0; i < nAmpl(); ++i) {
			for (size_t j = 0; j < nAmpl(); ++j) {
				integralMatrix[i][j] += integralMatrix_omp[ID][i][j];
				accCorrIntegralMatrix[i][j] += accCorrIntegralMatrix_omp[ID][i][j];
			}
		}
	}
	for (size_t i = 0; i < nAmpl(); ++i) {
		for (size_t j = 0; j < nAmpl(); ++j) {
			integralMatrix[i][j]        /= _nPoints;
			accCorrIntegralMatrix[i][j] /= _nPoints;
		}
	}
	if (not setIntegrals(integralMatrix, accCorrIntegralMatrix)) {
		utils::makeError("integrator::loadIntegrals(...)","Could not set integrated matrices.");
		return false;
	}
	return true;
}

bool integrator::loadIntegrals(const std::string& psFileName, const std::string& accFileName) {
	// Trust the content of the files, as long as the dimensions match
	std::pair<bool, complexMatrix > ps_load = utils::readMatrixFromTextFile<complex>(psFileName, _amplitudes.size());
	if (!ps_load.first) {
		utils::makeError("integrator::loadIntegrals(...)","Could not load ps matrix from '" + psFileName + "'.");
		return false;
	}
	std::pair<bool, complexMatrix> ac_load = utils::readMatrixFromTextFile<complex>(accFileName, _amplitudes.size());
	if (!ac_load.first) {
		utils::makeError("integrator::loadIntegrals(...)","Could not load ac matrix from '" + accFileName + "'.");
		return false;
	}
	if (not setIntegrals(ps_load.second, ac_load.second)) {
		utils::makeError("integrator::loadIntegrals(...)","Could not set loaded matrices.");
		return false;
	}
	return true;
}

bool integrator::setIntegrals(const complexMatrix& ps_integral, const complexMatrix& ac_integral) {
	if (ps_integral.size() != _nAmpl) {
		utils::makeError("integrator::setIntegrals(...)","ps_integral has the wrong size.");
		return false;
	}
	if (ac_integral.size() != _nAmpl) {
		utils::makeError("integrator::setIntegrals(...)","ac_integral has the wrong size.");
		return false;
	}
	for (size_t a = 0; a < _nAmpl; ++a) {
		if (ps_integral[a].size() != _nAmpl) {
			utils::makeError("integrator::setIntegrals(...)","ps_integral has the wrong size (line " +std::to_string(a) + ").");
			return false;
		}	
		if (ac_integral[a].size() != _nAmpl) {
			utils::makeError("integrator::setIntegrals(...)","ac_integral has the wrong size (line " +std::to_string(a) + ").");
			return false;
		}	
	}
	for (size_t a_i = 0; a_i < _nAmpl; ++a_i) {
		if (ps_integral[a_i][a_i].imag() != 0.) {
			utils::makeError("integrator::setIntegrals(...)","ps_integral has non-vanisihng imaginary part on the diagonal: " + utils::to_string(ps_integral[a_i][a_i]) + 
			                                                 "(wave '" + _amplitudes[a_i]->name() + "').");
			return false;
		}
		if (ac_integral[a_i][a_i].imag() != 0.) {
			utils::makeError("integrator::setIntegrals(...)","ac_integral has non-vanisihng imaginary part on the diagonal: " + utils::to_string(ac_integral[a_i][a_i]) + 
			                                                 "(wave '" + _amplitudes[a_i]->name() + "').");
			return false;
		}
		for (size_t a_j = 0; a_j < a_i; ++a_j) {
			if (norm(ps_integral[a_i][a_j] - std::conj(ps_integral[a_j][a_i])) > _numLim ) {
				utils::makeError("integrator::setIntegrals(...)","ps_integral is not hermitian: " + utils::to_string(ps_integral[a_i][a_j]) + " vs. " 
				                                                 +utils::to_string(ps_integral[a_j][a_i]) + ".");
				return false;
			}
			if (norm(ac_integral[a_i][a_j] - std::conj(ac_integral[a_j][a_i])) > _numLim ) {
				utils::makeError("integrator::setIntegrals(...)","ac_integral is not hermitian: " + utils::to_string(ac_integral[a_i][a_j]) + " vs. " 
				                                                 +utils::to_string(ac_integral[a_j][a_i]) + ".");
				return false;
			}
		}
	}
	bool nanAmpl = false;
	for (size_t a = 0; a < _nAmpl; ++a) {
		if (std::isnan(ac_integral[a][a].real()) || std::isnan(ac_integral[a][a].imag())) {
			nanAmpl = true;
			utils::makeError("integrator::setIntegrals(...)","NaN in ac integral for wave '" + _amplitudes[a]->name() + "': " + utils::to_string(ac_integral[a][a])+ ".");
		}
	}
       for (size_t a = 0; a < _nAmpl; ++a) {
		if (std::isnan(ps_integral[a][a].real()) || std::isnan(ps_integral[a][a].imag())) {
			nanAmpl = true;
			utils::makeError("integrator::setIntegrals(...)","NaN in ps integral for wave '" + _amplitudes[a]->name() + "': " + utils::to_string(ps_integral[a][a]) + ".");
		}
	}

	if (nanAmpl) {
		utils::makeError("integrator::setIntegrals(...)","NaN amplitudes encountered.");
		return false;
	}
	_integralMatrix = ps_integral;
	_accCorrIntegralMatrix =  ac_integral;
	if (not makeRealMatrices()) {
		utils::makeError("integrator::setIntegrals(...)","Could bit create real versions of the matrices.");
		return false;
	}
	_isIntegrated = true;
	return true;
}

bool integrator::makeRealMatrices() {
	_realIntegralMatrix = realMatrix(2*_nAmpl, realVector(2*_nAmpl, 0.));
	_realAccCorrMatrix  = realMatrix(2*_nAmpl, realVector(2*_nAmpl, 0.));
	for (size_t a_i = 0; a_i < _nAmpl; ++a_i) {
		for (size_t a_j = 0; a_j < _nAmpl; ++a_j) {
			double norm = pow((_integralMatrix[a_i][a_i] * _integralMatrix[a_j][a_j]).real(),.5);
			if (norm == 0.) {
				continue;
			}
			_realIntegralMatrix[2*a_i  ][2*a_j  ] = _integralMatrix[a_i][a_j].real()/norm;
			_realIntegralMatrix[2*a_i  ][2*a_j+1] =-_integralMatrix[a_i][a_j].imag()/norm;
			_realIntegralMatrix[2*a_i+1][2*a_j  ] = _integralMatrix[a_i][a_j].imag()/norm;
			_realIntegralMatrix[2*a_i+1][2*a_j+1] = _integralMatrix[a_i][a_j].real()/norm;

			_realAccCorrMatrix[2*a_i  ][2*a_j  ] = _accCorrIntegralMatrix[a_i][a_j].real()/norm;
			_realAccCorrMatrix[2*a_i  ][2*a_j+1] =-_accCorrIntegralMatrix[a_i][a_j].imag()/norm;
			_realAccCorrMatrix[2*a_i+1][2*a_j  ] = _accCorrIntegralMatrix[a_i][a_j].imag()/norm;
			_realAccCorrMatrix[2*a_i+1][2*a_j+1] = _accCorrIntegralMatrix[a_i][a_j].real()/norm;
		}
	}
	return true;
}

complexMatrix integrator::getIntegralMatrix(bool accCorr) const {
	if (not _isIntegrated) {
		utils::makeError("integrator::getIntegralMatrix()","Not integrated yet.");
	}
	if (accCorr){
		return _accCorrIntegralMatrix;
	}
	return _integralMatrix;
}

std::pair<bool, complex> integrator::element(size_t i, size_t j, bool accCorr) const {
	if (not _isIntegrated) {
		utils::makeError("integrator::element()","Not integrated yet.");
		return std::pair<bool, complex>(false, complex(0.,0.));
	}
	if (i >= nAmpl()) {
		utils::makeError("integrator::element()","First index too big.");
		return std::pair<bool, complex>(false, complex(0.,0.));
	}
	if (j >= nAmpl()) {
		utils::makeError("integrator::element()","Second index too big.");
		return std::pair<bool, complex>(false, complex(0.,0.));
	}
	if (accCorr){
		return std::pair<bool, complex>(true, _accCorrIntegralMatrix[i][j]);	
	}
	return std::pair<bool, complex>(true, _integralMatrix[i][j]);	
}

std::pair<bool, realVector> integrator::getNormalizations(bool accCorr) const {
	if (!_isIntegrated) {
		utils::makeError("integrator::getNormalizations(...)","Cannot get normalizations before integration.");
		return std::pair<bool, realVector>(false, realVector());
	}
	realVector retVal(_nAmpl);
	complex val;
	for (size_t a = 0; a < _nAmpl; ++a) {
		if (accCorr) {
			val = _accCorrIntegralMatrix[a][a];
		} else {
			val = _integralMatrix[a][a];
		}
		if (val.imag() != 0.) {
			utils::makeError("integrator::getNormalizations(...)","non-vanishing imaginary part on the diagonal (wave: '" + _amplitudes[a]->name() + "').");
			return std::pair<bool, realVector>(false, realVector());
		}
		if (val.real() == 0.) {
			retVal[a] = 0.;	
		} else {
			retVal[a] = 1./pow(val.real(), .5);
		}
	}
	return std::pair<bool, realVector>(true, retVal);
}

double integrator::totalIntensity(const complexVector& prodAmpl, bool accCorr) const {
	if (not _isIntegrated) {
		utils::makeError("integrator::totalIntensity(...)","Not integrated yet. Returning zero.");
		return 0.;
	}
	if (prodAmpl.size() != nAmpl()) {
		utils::makeError("integrator::totalIntensity(...)","Number of production amplitudes ("+std::to_string(prodAmpl.size())+") does not match (should be "+std::to_string(nAmpl())+"). Returning zero.");
		return 0.;
	}
//	const std::complexVector &integralMatrix = accCorr ? _accCorrIntegralMatrix : _integralMatrix;
//	complex retVal(0.,0.);
//	std::cout << "-------------------------------------" << std::endl;
//	for (size_t i = 0; i < nAmpl(); ++i) {
//		for (size_t j = 0; j < nAmpl(); ++j) {
//			double norm = pow(_integralMatrix[i][i].real() * _integralMatrix[j][j].real(), .5); // Always normalized to phaseSpace
//			if (norm == 0.) {
//				continue;
//			}
//			if (std::conj(prodAmpl[i]) * integralMatrix[i][j] * prodAmpl[j]/norm !=  0.) {
//				std::cout << i << " " << j << "resse "  <<  std::conj(prodAmpl[i]) << " " << integralMatrix[i][j]/norm << " " << prodAmpl[j] << std::endl;
//			}
//			retVal += std::conj(prodAmpl[i]) * integralMatrix[i][j] * prodAmpl[j]/norm;			
//		}
//	}
	double alternative = 0.;
	const double *reIm = (double*)&prodAmpl[0]; // Trust this at the moment
	const realMatrix &real_matrix = accCorr ? _realAccCorrMatrix : _realIntegralMatrix;
	for (size_t a_i = 0; a_i < 2*_nAmpl; ++a_i) {
		for (size_t a_j = 0; a_j < 2*_nAmpl; ++a_j) {
//			if ( reIm[a_i]*realMatrix[a_i][a_j]*reIm[a_j] != 0.) {
//				std::cout << a_i << " " << a_j << "innse " << reIm[a_i]<< " " <<rea_matrix[a_i][a_j] << " " << reIm[a_j] << std::endl;
//			}
			alternative += reIm[a_i]*real_matrix[a_i][a_j]*reIm[a_j];
		}
	}
	return alternative;	
}

realVector integrator::DtotalIntensity(const complexVector& prodAmpl, bool accCorr) const {
// Return a vector of length 2*nAmpl() for the derivative w.r.t. re, and im : {dRe0, dIm0, dRe1, ..., dImnAmpl()}
	if (not _isIntegrated) {
		utils::makeError("integrator::DtotalIntensity(...)","Not integrated yet. Returning empty vector.");
		return realVector();
	}
	if (prodAmpl.size() != nAmpl()) {
		utils::makeError("integrator::DtotalIntensity(...)","Number of production amplitudes does not match. Returning zero.");
		return realVector();
	}
	realVector retVal(2*nAmpl(), 0.);
	const realMatrix &matrix = accCorr ? _realAccCorrMatrix : _realIntegralMatrix;
	const double* reIm=  (double*)&prodAmpl[0];
	for (size_t a_i = 0; a_i < 2*_nAmpl; ++a_i) {
		for (size_t a_j = 0; a_j < 2*_nAmpl; ++ a_j) {
			retVal[a_i] += (matrix[a_i][a_j] + matrix[a_j][a_i])*reIm[a_j];
		}
	}

//	const complexMatrix &integralMatrix = accCorr ? _accCorrIntegralMatrix : _integralMatrix;
//	for (size_t j = 0; j < nAmpl(); ++j) {
//		for (size_t i = 0; i < nAmpl(); ++i) {
//			double norm = pow(_integralMatrix[i][i].real() * _integralMatrix[j][j].real(), .5); // Always normalized to phaseSpace
//			if (norm == 0.) {
//				continue;
//			}
//			complex  factor = integralMatrix[i][j] * prodAmpl[j]/norm*2.;
//			retVal[2*i  ] += factor.real();
//			retVal[2*i+1] += factor.imag();
//
//		}
//	}
	return retVal;
}

realMatrix integrator::DDtotalIntensity(const complexVector& prodAmpl, bool accCorr) const {
	if (not _isIntegrated) {
		utils::makeError("integrator::DDtotalIntensity(...)","Not integrated yet. Returning empty vector.");
		return realMatrix();
	}
	if (prodAmpl.size() != nAmpl()) {
		utils::makeError("integrator::DDtotalIntensity(...)","Number of production amplitudes does not match. Returning zero.");
		return realMatrix();
	}
	realMatrix retVal(2*nAmpl(), realVector(2*nAmpl(),0.));
	const complexMatrix &integralMatrix = accCorr ? _accCorrIntegralMatrix : _integralMatrix;
	for (size_t ai = 0; ai < nAmpl(); ++ai) {
		for (size_t aj = 0; aj < nAmpl(); ++aj) {
			double norm = pow(_integralMatrix[ai][ai].real() * _integralMatrix[aj][aj].real(), .5); // Always normalized to phaseSpace
			if (norm == 0.) {
				continue;
			}
			complex factor = integralMatrix[ai][aj]/norm*2.;
			retVal[2*ai  ][2*aj  ] += factor.real();
			retVal[2*ai  ][2*aj+1] -= factor.imag();
			retVal[2*ai+1][2*aj  ] += factor.imag();
			retVal[2*ai+1][2*aj+1] += factor.real();
		}
	}
	return retVal;
}

bool integrator::writeToFile(const std::string& fileName, bool accCorr) const {
	if (not _isIntegrated) {
		utils::makeError("integrator::writeToFile(...)","Not integrated yet. Returning empty vector.");
		return false;
	}
	std::ofstream outFile;
	outFile.open (fileName.c_str());
	outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1);
	for (size_t i = 0; i < nAmpl(); ++i) {
		for (size_t j = 0; j < nAmpl(); ++j) {
			if (accCorr) {
				outFile << _accCorrIntegralMatrix[i][j] << " ";
			} else {
				outFile << _integralMatrix[i][j] << " ";
			}
		}
		outFile << std::endl;
	}
	outFile.close();
	return true;
}

bool integrator::makeIntegrals(const std::string& ps_fileName, const std::string& ac_fileName, const std::string& ID_fileName) {
	if (_isIntegrated) {
		utils::makeWarning("integrator::makeIntegrals(...)","Already integrated. Do nothing.");
		return true;
	}
	if (loadIntegrals(ps_fileName, ac_fileName)) {
		utils::makeInfo( "integrator::makeIntegrals(...)","Integrals successfully loaded from '" + ps_fileName + "' and '" + ac_fileName + "'.");
		if (!checkIDfile(ID_fileName)) {
			utils::makeError("intgrator::makeIntegrals(...)","Amplitudes do not seem to match the ID file '" + ID_fileName + "'. Integral files might be out of date.");
			return false;
		}
		return true;
	}
	if (!integrate()) {
		utils::makeError("integrator::makeIntegrals(...)","Integration failed.");
		return false;
	} else {
		utils::makeInfo("integrator::makeIntegrals(...)","Integration successful.");
	}
	if (!writeToFile(ps_fileName, false)) {
		utils::makeError("integrator::makeIntegrals(...)","Could not write phase-space integrals to '" + ps_fileName + "'.");
		return false;
	}
	if (!writeToFile(ac_fileName, true)) {
		utils::makeError("integrator::makeIntegrals(...)","Could not write acceptance corrected integrals to '"+ac_fileName+"'.");
		return false;
	}
	if (!writeIDfile(ID_fileName)) {
		utils::makeError("integrator::makeIntegrals(...)","Could not write the ID file to: '"+ID_fileName+"'.");
		return false;
	}
	utils::makeInfo("integrator::makeIntegrals(...)","Integrals successfully written to '"+ps_fileName+"' and '"+ac_fileName+"'.");
	return true;
}

bool integrator::setNpoints(size_t n) {
	_nPoints = n;
	return true;
}

std::pair<bool, std::string> integrator::getWaveName(size_t i) const {
	if (i >= nAmpl()) {
		return std::pair<bool, std::string>(false, "");
	}
	std::string waveName = _amplitudes[i]->name();
	return std::pair<bool, std::string>(true, waveName);
}

std::pair<bool, size_t> integrator::getNpoints() const {
	return std::pair<bool, size_t>(true, _nPoints);
}

bool integrator::writeIDfile(const std::string& outFileName) const {
	const realVector checkKin = _generator->generate();

	std::ofstream outFile;
	outFile.open (outFileName.c_str());
	outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1);
	for (double var : checkKin) {
		outFile << var << std::endl;
	}
	outFile << _efficiency->eval(checkKin) << std::endl;

	for (const std::shared_ptr<amplitude> a : _amplitudes) {
		outFile << a->eval(checkKin) << std::endl;
	}
	for (const std::shared_ptr<amplitude> a : _amplitudes) {
		outFile << a->name() << std::endl;
	}
	outFile.close();
	return true;
}

bool integrator::checkIDfile(const std::string& inFileName) const {
	const size_t nKin = _generator->nKin();

	realVector checkKin(nKin, 0.);
	std::ifstream fin(inFileName.c_str(), std::ifstream::in);
	
	double val;
	for (size_t k = 0; k < nKin; ++k) {
		if (!(fin >> val)) {
			utils::makeError("integrator::checkIDfile(...)","Could not load kinematics "+std::to_string(k)+" from '"+inFileName+"'.");
			return false;
		}
		checkKin[k] = val;
	}
	if (!(fin >> val )) {
		utils::makeError("integrator::checkIDfile(...)","Could not efficiency value from '"+inFileName+"'.");	
	}
	if (abs(_efficiency->eval(checkKin) - val) > _numLim) {
			utils::makeError("integrator::checkIDfile(...)","Efficiency value for does not match (" +std::to_string(_efficiency->eval(checkKin))+" vs. "+std::to_string(val)+ ").");
			return false;
	}

	complex amplIn;
	for (size_t a = 0; a < _nAmpl; ++a) {
		if (!(fin >> amplIn)) {
			utils::makeError("integrator::checkIDfile(...)","Could not amplitude value for "+_amplitudes[a]->name()+" (index "+std::to_string(a)+") from '"+inFileName+"'.");
			return false;
		}
		if (std::norm(_amplitudes[a]->eval(checkKin) - amplIn) > _numLim) {
			utils::makeError("integrator::checkIDfile(...)","Amplitude value for "+_amplitudes[a]->name()+" do not match ("
			                                                +utils::to_string(_amplitudes[a]->eval(checkKin))+" vs. "+utils::to_string(amplIn)+").");
			return false;
		}
	}
	std::string nameIn;
	for (size_t a = 0; a < _nAmpl; ++a) {
		if ( fin >> nameIn ) {
			if (nameIn != _amplitudes[a]->name()) {
				utils::makeError("integrator::checkIDfile(...)","Amplitude names do not match: read '" + nameIn + "' should be '" + _amplitudes[a]->name() + "' (wave index " +std::to_string(a) + ").");
				return false;
			}
		} else {
			if (a == 0) {
				utils::makeWarning("integrator::checkIDfile(...)","Could not read any ampliutude names. Probably a legacy file. First amplitude should be: '" + _amplitudes[a]->name() + "'.");
				break;
			} else {
				utils::makeError("integrator::checkIDfile(...)","Could not read a name for amplitude '" + _amplitudes[a]->name() + "' (wave index " + std::to_string(a) + ").");
				return false;
			}
		}
	}
	fin.close();
	return true;
}

bool integrator::setCoherenceBorders(sizeVector& borders) {
	if (borders[0] != 0) {
		utils::makeError("integrator::setCoherenceBorders(...)","First border has to be zero.");
		return false;
	}
	_amplitudeCoherenceBorders = sizeVector(borders.size(),0);
	for (size_t i = 1; i < borders.size(); ++i) {
		if (borders[i] <= _amplitudeCoherenceBorders[i-1]) {
			utils::makeError("integrator::setCoherenceBorders(...)","Borders are nor ordered.");
			return false;
		}
		 _amplitudeCoherenceBorders[i] = borders[i];
	}
	if (_amplitudeCoherenceBorders[_amplitudeCoherenceBorders.size()] != _nAmpl) {
		utils::makeError("integrator::setCoherenceBorders(...)","Last border has to be _nAmpl.");
		return false;
	}
	return true;
}

bool integrator::addIncoherentSector(std::shared_ptr<integrator> sector) {
	if (not isIntegrated()) {
		utils::makeError("integrator::addIncoherentSector(...)","Cannot add second sector before integration.");
		return false;
	}
	if (not sector->isIntegrated()) {
		utils::makeError("integrator::addIncoherentSector(...)","Cannot add non-integrated sector.");
		return false;
	}
	if (not (*_kinSignature == *(sector->kinSignature()))) {
		utils::makeError("integrator::addIncoherentSector(...)","Kinematic signatures have to match.");
		return false;
	}
	if (_nPoints != sector->getNpoints().second) {
		utils::makeError("integrator::addIncoherentSector(...)","Number of points has to match.");
		return false;
	}
	
	size_t nNew   = sector->nAmpl();
	size_t nTotal = _nAmpl + nNew;
	complexMatrix new_ps_integral(nTotal, complexVector(nTotal, complex(0.,0.)));
	complexMatrix new_ac_integral(nTotal, complexVector(nTotal, complex(0.,0.)));
	for (size_t i = 0; i < _nAmpl; ++i) {
		for (size_t j = 0; j < _nAmpl; ++j) {
			new_ps_integral[i][j] = _integralMatrix[i][j];
			new_ac_integral[i][j] = _accCorrIntegralMatrix[i][j];
		}	
	}
	for (size_t i = 0; i < nNew; ++i) {
		for (size_t j = 0; j < nNew; ++j) {
			new_ps_integral[_nAmpl+i][_nAmpl+j] = sector->element(i,j,false).second;
			new_ac_integral[_nAmpl+i][_nAmpl+j] = sector->element(i,j,true).second;
		}
	}
	for (size_t i : sector->getCoherenceBorders()) {
		if (i > 0) { // Do not add the leading 0
			_amplitudeCoherenceBorders.push_back(_nAmpl + i);
		}
	}
	_nAmpl = nTotal;
	if (not setIntegrals(new_ps_integral, new_ac_integral)) {
		utils::makeError("integrator::addIncoherentSector(...)","Could not set new integral matrices.");
		return false;
	}
	std::vector<std::shared_ptr<amplitude> > newAmplitudes;
	for (std::shared_ptr<amplitude> a : _amplitudes) {
		newAmplitudes.push_back(a);
	}
	for (std::shared_ptr<amplitude> a : (*sector)._amplitudes) {
		newAmplitudes.push_back(a);
	}
	_amplitudes = newAmplitudes;
	return true;	
}
