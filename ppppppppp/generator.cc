#include"generator.h"
#include<iostream>
#include"utils.h"
generator::generator() : _kinSignature(std::make_shared<kinematicSignature>(0)), _maxFail(1000), _failCount(0) {}

realVector generator::generate() const {
	std::cerr << "generator::generate(): ERROR: Called base class method. Returning empty vector." << std::endl;
	return realVector();
}

threeParticleMassGenerator::threeParticleMassGenerator(double initialMass, const realVector& fsMasses, std::shared_ptr<kinematicSignature> kinSig):
	generator(), _initialMass(initialMass), _fsMasses(fsMasses) {
		if (_fsMasses.size() != 3) {
			std::cerr << "threeParticleMassGenerator::threeParticleMassGenerator(...): Not three fsMasses given: " << _fsMasses.size() << std::endl;
			throw;
		}
		_kinSignature = kinSig;
}

double s3v(const realVector& v) {
	double retVal = v[0]*v[0];
	retVal -= v[1]*v[1];
	retVal -= v[2]*v[2];
	return retVal;
}

bool threeParticleMassGenerator::isValidPoint(const realVector& kin) const {
	const double s   = kin[0];
	const double s12 = kin[1];
	const double s13 = kin[2];
	const double m12 = pow(_fsMasses[0],2);
	const double m22 = pow(_fsMasses[1],2);
	const double m32 = pow(_fsMasses[2],2);
	const double s23 = s + m12 + m22 + m32 - s12 - s13;
	if (s23 < pow(_fsMasses[1]+_fsMasses[2],2)) {
		return false;
	}
	if (s23 > pow(_initialMass - _fsMasses[0], 2)) {
		return false;
	}
	const double p12 = utils::breakupMomentumSquared(s, s23, m12);
	if (p12 < 0.) {
		return false;
	}
	const double p22 = utils::breakupMomentumSquared(s, s13, m22);
	if (p22 < 0.) {
		return false;
	}
	const double p32 = utils::breakupMomentumSquared(s, s12, m32);
	if (p32 < 0.) {
		return false;
	}
	const double q12 = utils::breakupMomentumSquared(s12, m12, m22);
	if (q12 < 0.) {
		return false;
	}
	const double q13 = utils::breakupMomentumSquared(s13, m12, m32);
	if (q13 < 0.) {
		return false;
	}
	const double q23 = utils::breakupMomentumSquared(s23, m22, m32);
	if (q23 < 0.) {
		return false;
	}

	const double p1 = pow(p12, .5);
	const double p2 = pow(p22, .5);
	const double p3 = pow(p32, .5);

	if (p1 + p2 < p3) {
		return false;
	}
	if (p1 + p3 < p2) {
		return false;
	}
	if (p2 + p3 < p1) {
		return false;
	}
	return true;
}

realVector threeParticleMassGenerator::generate() const {
	realVector retVal(3, _initialMass*_initialMass);
	double sMin12 = std::pow(_fsMasses[0] + _fsMasses[1],2);
	double sMin13 = std::pow(_fsMasses[0] + _fsMasses[2],2);

	double sMax12 = std::pow(_initialMass - _fsMasses[2],2);
	double sMax13 = std::pow(_initialMass - _fsMasses[1],2);
	while (true) {
		double s12  = sMax12 - sMin12;
		s12 *= utils::random();
		s12 += sMin12;
		retVal[1] = s12;
	
		double s13 = sMax13 - sMin13;
		s13 *= utils::random();
		s13 += sMin13;
		retVal[2] = s13;

		if (isValidPoint(retVal)) {
			_failCount = 0;
			break;
		} else {
			++_failCount;
			if (_failCount > _maxFail) {
				throw;
			}
		}
	}
	return retVal;
}

dalitzScanGenerator::dalitzScanGenerator(double            initialMass, 
		                         const realVector& fsMasses, 
		                         const realVector& massesX,
		                         const realVector& massesY,
		                         std::shared_ptr<kinematicSignature> kinSig):
	threeParticleMassGenerator(initialMass, fsMasses, kinSig),
	indexX(0), indexY(0), massesX(massesX), massesY(massesY) {}

bool dalitzScanGenerator::ended() const {
	if (indexX == 0) {
		if (indexY == massesY.size()) {
			return true;
		}
	}
	return false;
}

realVector dalitzScanGenerator::getNext() const {
	if (ended()) {
		return realVector();
	}

	double massX = massesX[indexX];
	double massY = massesY[indexY];

	++indexX;
	if (indexX == massesX.size()) {
		indexX = 0;
		++indexY;
	}
	return {_initialMass * _initialMass, massX * massX, massY * massY};
}

realVector dalitzScanGenerator::generate() const {
	realVector retVal;
	bool hasPoint = false;
	while (true) {
		retVal = getNext();
		if (retVal.size() == 0) {
			break;
		}
		if (isValidPoint(retVal)) {
			hasPoint = true;
			break;
		}
	}
	if (!hasPoint) {
		utils::makeWarning("dalitzScanGenerator::generate()","No more point to be generated. Returning empty vector.");
	}
	return retVal;
}
