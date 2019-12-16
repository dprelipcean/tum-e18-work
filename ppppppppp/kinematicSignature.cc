#include<iostream>
#include"kinematicSignature.h"

kinematicSignature::kinematicSignature(size_t identifier) : _identifier(identifier) {};

bool kinematicSignature::operator==(const kinematicSignature& other) const {
	return _identifier == other._identifier;
}

realMatrix kinematicSignature::getBoseSymmetrizedKinematics(const realVector& kin) const {
	if (kin.size() != nKin()) {
		std::cerr << "kinematicSignature::getBoseSymmetrizedKinematics(...): ERROR: Number of kinematic variables does not match. Returning empty vetor." << std::endl;
		return realMatrix();
	}
	switch(_identifier) {
		case 0: {return realMatrix();}
		case 1: {realMatrix retVal(2); retVal[0] = {kin[0], kin[1], kin[2]}; retVal[1] = {kin[0], kin[2], kin[1]}; return retVal;}
		case 2: {realMatrix retVal(1); retVal[0] = {kin[0], kin[1], kin[2]}; return retVal;}
	}
	std::cerr << "kinematicSignature::getBoseSymmetrizedKinematics(): ERROR: Unknown identifier: " << _identifier << ". Returning empty vector" <<  std::endl;
	return realMatrix();
}

size_t kinematicSignature::nKin() const {
	switch (_identifier) {
		case 0: return 0;
		case 1: return 3;
		case 2: return 3;
	}
	std::cerr << "kinematicSignature::nKin(): ERROR: Unknown identifier: " << _identifier << std::endl;
	return 0;
}

sizeVector kinematicSignature::isobarMassIndices() const { 
	switch (_identifier) {
		case 0: return sizeVector();
		case 1: return sizeVector(1,1);
		case 2: std::cout << "kinematicSignature::isobarMassIndices(): ERROR: Not properly implemented for _identifier = 2." << std::endl; return sizeVector();
	}
	std::cerr << "kinematicSignature::nKin(): ERROR: Unknown identifier: " << _identifier << std::endl;
	return sizeVector();
}

void kinematicSignature::print() const {
	switch (_identifier) {
		case 0: std::cout << "Empty kinematic signature" << std::endl; return;
		case 1: std::cout << "Signature for one initial state going into three final state particles.\nThe kinematic variables are {s, s_{12}, s_{13}\nParticles 2 and 3 are assumed to be identical" <<  std::endl; return;
		case 2: std::cout << "Signature for one initial state going into three final state particles.\nThe kinematic variables are {s, s_{12}, s_{13}\nNo particles are identical" << std::endl; return;
	}
	std::cerr << "kinematicSignature::print(): ERROR: Unknown identifier: " << _identifier << std::endl;
}
