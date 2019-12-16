#include"massShape.h"
#include"utils.h"
#include<iostream>
#include<math.h>

massShape::massShape(std::string name, realVector pars, std::vector<std::string> parNames):
	_nPar(pars.size()), _name(name), _parameters(pars) {
	if (parNames.size() == _nPar) {
		_parameterNames = parNames;
	} else {
		if (parNames.size() > 0) {
			std::cerr << "massShape::massShape(...): ERROR: Number of parameter names does not match and is not zero." << std::endl;
		}
		_parameterNames = std::vector<std::string>();
		for (size_t i = 0; i < _nPar; ++i) {
			_parameterNames.push_back(std::string("unname_parameter_")+std::to_string(i));
		}
	}
}

complex massShape::eval(double m) const {
	std::cerr << "massShape::eval(double): ERROR: Evaluating the bas class at mass m = " << m << ", returning zero." << std::endl;
	return complex(0.,0.);
}

bool massShape::setParameter(size_t n, double val) {
	if (n < _nPar) {
		_parameters[n] = val;
		return true;
	}
	return false;
}
	
std::pair<bool, double> massShape::getParameter(size_t n) const {
	if (n<_nPar){
		return std::pair<bool, double>(true, _parameters[n]);	
	}
	return std::pair<bool, double>(false, 0.);
}

bool massShape::setParName(size_t n, std::string name) {
	if(n<_nPar) {
		_parameterNames[n] = name;
		return true;
	}
	return false;
}

std::pair<bool, std::string> massShape::getParName(size_t n) const {
	if (n<_nPar) {
		return std::pair<bool, std::string>(true, _parameterNames[n]);
	}
	return std::pair<bool, std::string>(false, "");
}
//------------------------------------------------------------------------------
simpleBW::simpleBW(double mass, double width):
	massShape("simpleBW", {mass, width}, {"m0","G0"}) {}

complex simpleBW::eval(double s) const {
	double num = _parameters[0] * _parameters[1];
	complex den = complex(_parameters[0]*_parameters[0] - s, -_parameters[0]*_parameters[1]);
	return num/den;
}
//------------------------------------------------------------------------------
stepLike::stepLike(double sMin, double sMax) :
	massShape("stepLike|" + std::to_string(sMin) + "|" + std::to_string(sMax) + "|", {sMin, sMax}, {"sMin","sMax"}) {}

complex stepLike::eval(double s) const {
	if (s > _parameters[0] and s <= _parameters[1]) {
		return complex(1.,0.);
	}
	return complex(0.,0.);
}
//------------------------------------------------------------------------------
sawtooth::sawtooth(double sMin, double sMax) :
	massShape("sawtooth|" + std::to_string(sMin) + "|" + std::to_string(sMax) + "|", {sMin, sMax}, {"sMin","sMax"}) {}

complex sawtooth::eval(double s) const {
	if (s > _parameters[0] and s <= _parameters[1]) {
		double x = (s - _parameters[0])/(_parameters[1] - _parameters[0]);
		return complex(x - .5,0.);
	}
	return complex(0.,0.);
}
//------------------------------------------------------------------------------
spike::spike(double sMin, double sEval, double sMax) :
	massShape("spike|"+std::to_string(sMin)+"|"+std::to_string(sEval)+"|"+std::to_string(sMax)+"|",
	          {sMin,sEval,sMax},{"sMin","sEval","sMax"}) {}

complex spike::eval(double s) const {
	if (s <= _parameters[0] or s > _parameters[2]) {
		return complex(0.,0.);
	}
	if (s <= _parameters[1]) {
		double x = (s-_parameters[0])/(_parameters[1]-_parameters[0]);
		return complex(x,0.);
	}
	double x = (s-_parameters[1])/(_parameters[2] - _parameters[1]);
	return complex(1.-x,0.);
}
//------------------------------------------------------------------------------
spline_deg::spline_deg(double sMin, double sMax, size_t deg) :
	massShape("spline_deg"+std::to_string(deg)+"|" + std::to_string(sMin) + "|" + std::to_string(sMax) + "|", {sMin, sMax}, {"sMin","sMax"}), _deg(deg) {
	if (deg > 5) {
		utils::makeError("spline_deg::spline_deg(...)","Polynominal degree > 5 not implemented. (Maximum spline degree is 4)");
		throw;
	}
}

complex spline_deg::eval(double s) const {
	if (s < _parameters[0] or s > _parameters[1]) { // Be aware, that the lower limit is included here (needed in the matching).
		return complex(0.,0.);
	}
	if (_deg == 0) {
		return complex(1.,0.);
	}
	double x = (s-_parameters[0])/(_parameters[1]-_parameters[0]);
	double retVal = 0.;
	if (_deg == 1) {
		retVal = x -.5;
	} else if (_deg == 2) {
		retVal = x*x - x + 1./6.;
	} else if (_deg == 3) {
		retVal = x*x*x - 3./2.*x*x + 3./5. *x -1./20;
	} else if (_deg == 4) {
		retVal = x*x*x*x -2. *x*x*x + 9./7. * x*x - 2./7.*x +1./70.;
	} else if (_deg ==5) {
			retVal = x*x*x*x*x - 5./2. *x*x*x*x + 20./9. *x*x*x - 5./6.*x*x + 5./42.*x - 1./252.;
	}
	return complex(retVal, 0.);
}

double spline_deg::Deval(double s, size_t nth) const {
	if (s < _parameters[0] or s > _parameters[1]) { // Be aware, that the lower limit is included here (needed for the spline matching).
		return 0.;
	} 
	if (nth > _deg) {
		return 0.;
	}

	double retVal = 0.;
	double dxds = 1./(_parameters[1]-_parameters[0]); // dxds = binWidth
	double x    = (s-_parameters[0])*dxds;
	if (nth == 0) {
		retVal = eval(s).real();
	} else if (nth == 1) {
		if (_deg == 1) {
			retVal = dxds;
		} else if (_deg == 2) {
			retVal = dxds * (2*x - 1.);
		} else if (_deg == 3) {
			retVal = dxds * (3*x*x - 3.*x + 3./5.);
		} else if (_deg == 4) {
			retVal = dxds * (4*x*x*x - 6. * x*x + 18./7.* x - 2./7.);
		} else if (_deg == 5) {
			retVal = dxds * (5 *x*x*x*x - 10.*x*x*x + 20./3. *x*x - 5./3. *x + 5./42);	
		}
	} else if (nth == 2) {
		if (_deg == 2) {
			retVal = 2*dxds*dxds;
		} else if (_deg == 3) {
			retVal = dxds*dxds*(6.*x - 3.);
		} else if (_deg == 4) {
			retVal = dxds*dxds*(12.*x*x - 12.*x + 18./7.);
		} else if (_deg == 5) {
			retVal = dxds*dxds*(20.*x*x*x - 30.*x*x + 40./3. * x - 5./3.);
		}
	} else if (nth == 3) {
		if (_deg == 3) {
			retVal = dxds*dxds*dxds*6.;
		} else if (_deg == 4) {
			retVal = dxds*dxds*dxds*(24.*x - 12.);
		} else if (_deg == 5) {
			retVal = dxds*dxds*dxds*(60.*x*x - 60.*x + 40./3.);
		}
	} else if (nth == 4) {
		if (_deg == 4) {
			retVal = pow(dxds,4) * 24.;
		} else if (_deg == 5) {
			retVal = pow(dxds,4) * (120.*x - 60.);
		}
	} else if (nth == 5) {
		retVal = pow(dxds,5) * 120.;
	} else {
		utils::makeError("spline_deg::Deval(...)",std::to_string(nth) + "th derivative not defined for a _deg = " + std::to_string(_deg) + " polynomial. Returning 0.");
	}
	if (std::isnan(retVal)) {
		std::cout << "ich wurde nan bei: s = " << s << "; x = " << x << " [" << _parameters[0] << "," << _parameters[1] << ") _deg = " << _deg << " D = " << nth << std::endl;
	}
	return retVal;
}
//------------------------------------------------------------------------------
constant::constant() :
	massShape("constant", {}, {}) {};

complex constant::eval(double s) const {
	return complex(s/s,0.);
}
//------------------------------------------------------------------------------
zeroMode0pp::zeroMode0pp(double s, double m2) :
	massShape("zeroMode0pp", {s, m2}, {"s", "m2"}) {};

complex zeroMode0pp::eval(double s12) const {
	double retVal = (_parameters[0] - 3*s12 + 3*_parameters[1])/8;
//	std::cout << "Called zeroMode0pp with " << s12 << " giving " << retVal << std::endl;
	return complex(retVal, 0.);
}
//------------------------------------------------------------------------------
zeroMode1mm::zeroMode1mm(double s, double m2) :
	massShape("zeroMode0pp", {s, m2}, {"s", "m2"}) {};

complex zeroMode1mm::eval(double s12) const {
	double retVal = _parameters[0]/(_parameters[0] + s12 - _parameters[1]);
//	std::cout << "Called zeroMode1mm with " << s12 << " giving " << retVal << std::endl;
	return complex(retVal, 0.);
}
//------------------------------------------------------------------------------
polynomialMassShape::polynomialMassShape(complexVector coefficients, double baseExponent) :
	massShape(std::string("polynomialMassShape_deg") + std::to_string(coefficients.size()-1), {}, {}), _polDeg(0), _baseExponent(baseExponent) {

	if (!utils::checkComplexDouble()) {
		std::cout << "polynomialMassShape::polynomialMassShape(...): ERROR: complex* is not double real, double imag array" << std::endl;
		throw;
	}
	for (size_t c = 0; c < coefficients.size(); ++c) {
		_parameters.push_back(coefficients[c].real());
		_parameterNames.push_back(std::string("c_") + std::to_string(c) + std::string("_r"));
		_parameters.push_back(coefficients[c].imag());
		_parameterNames.push_back(std::string("c_") + std::to_string(c) + std::string("_i"));
	}
	_nPar   = 2*coefficients.size();
	_polDeg = coefficients.size();
}

complex polynomialMassShape::eval(double s12) const {
	complex retVal(0.,0.);
	complex* coeffs = (complex*)&_parameters[0];
	for (size_t c = 0; c < _polDeg; ++c) {
		retVal += coeffs[c] * pow(s12, _baseExponent*c);
	}
	return retVal;
}
//------------------------------------------------------------------------------
BELLEbreitWigner::BELLEbreitWigner(std::string name, double mass, double width, size_t spin, double motherMass, double bachelorMass, double daughterMass1, double daughterMass2) :
	massShape("BELLEbw_"+name, {mass, width}, {name+"_mass", name+"_width"}), _spin(spin), _motherMass(motherMass), _bachelorMass(bachelorMass), _daughterMass1(daughterMass1), _daughterMass2(daughterMass2), _Rr(1.5), _RD(5.)
{
//	if (_motherMass - _bachelorMass < _parameters[0]) {
//		std::cout  << "BELLEbreitWigner::BELLEbreitWigner(...): ERROR: On shell resonance mass of '" << _name << "' too heavy for decay of mother particle: " << _motherMass << " -> " << _parameters[0] << " + " << _bachelorMass << std::endl;
//		throw;
//	}
	if (_daughterMass1 + _daughterMass2 > _parameters[0]) {
		std::cout << "BELLEbreitWigner::BELLEbreitWigner(...): ERROR: On shell resonance mass of '" << _name << "' too light for decay into daughter particles: " << _parameters[0] << " -> " << _daughterMass1 << " + " << _daughterMass2 << std::endl;
		throw;
	}	

	if (_spin > 2) {
		std::cout << "BELLEbreitWigner::BELLEbreitWigner(...): ERROR: Spin > 2 not supportet yet" << std::endl;
		throw;
	}
}

complex BELLEbreitWigner::eval(double s12) const {

	const double m12   = pow(s12, .5);

	if ((m12 < _daughterMass1 + _daughterMass2) || (m12 > _motherMass - _bachelorMass)) {
		return complex(0.,0.);
	}

	const double S          = _motherMass*_motherMass;
	const double sr         = _parameters[0]*_parameters[0];
	const double sDaughter1 = _daughterMass1*_daughterMass1;
	const double sDaughter2 = _daughterMass2*_daughterMass2;
	const double sBatch     = _bachelorMass*_bachelorMass;

	const double pr  = pow(pow(sr - sDaughter1 - sDaughter2, 2) - 4*sDaughter1*sDaughter2, .5)/2/_parameters[0];
	const double pAB = pow(pow(s12 - sDaughter1 - sDaughter2, 2) - 4*sDaughter1*sDaughter2, .5)/2/m12;

//	const double pD   = pow(pow(S - sr - sBatch, 2) - 4*sr*sBatch, .5)/2/_motherMass;
	const double pABC = pow(pow(S - s12 - sBatch, 2) - 4*s12*sBatch, .5)/2/_motherMass;

	double Fr = 1.;
	double FD = 1.;
	if (_spin == 1) {
		Fr = pow((pow(_Rr*pr,2)+1)/(pow(_Rr*pAB,2)+1),.5);
//		FD = pow((pow(_RD*pD,2)+1)/(pow(_RD*pABC,2)+1),.5);
		FD = pow(1./(pow(_RD*pABC,2)+1),.5);
	} else if (_spin == 2) {
		const double xr =  _Rr*_Rr*pr *pr;
		const double xAB = _Rr*_Rr*pAB*pAB;
		Fr = pow((pow(xr-3.,2) + 9*xr)/(pow(xAB-3.,2) + 9*xAB),.5);

//		const double xD   = _RD*_RD*pD  *pD;
		const double xABC = _RD*_RD*pABC*pABC;
//		FD = pow((pow(xD-3.,2) + 9*xD)/(pow(xABC-3.,2) + 9*xABC),.5);
		FD = pow(1./(pow(xABC-3.,2) + 9*xABC),.5);
	}
	
	const double Gamma = _parameters[1]* _parameters[0]/m12 * Fr*Fr * pow(pAB/pr, 2*_spin+1);
	
	complex retVal = complex(Fr*FD,0.)/complex(_parameters[0]*_parameters[0] - s12, - _parameters[0]*Gamma);
//	if (std::isnan(retVal.real()) || std::isnan(retVal.imag())) {
//		std::cout << "s12 = " << s12 << "; sDaughter1 = " << sDaughter1 << "; sDaughter2 = " << sDaughter2 << std::endl;
//		std::cout << "pr = " << pr << "; pAB = " << pAB << std::endl;
//		std::cout << "Fr = " << Fr << "; FD = " << FD << "; Gamma = " << Gamma << std::endl;
//	}
	return retVal;
}
//------------------------------------------------------------------------------
BELLE_LASS_KpiS::BELLE_LASS_KpiS(const realVector& parameters, double mPi, double mKs, bool usePS) : massShape("BELLE_LASS_KpiS", parameters, {"a", "r", "M0", "G0", "phiF", "phiR", "phiRsin", "F", "R", "MMax"}), _usePS(usePS), _mPi(mPi), _mKs(mKs) {}

complex BELLE_LASS_KpiS::eval(double s12) const {
	double M  = pow(s12, .5);
	complex retVal(0.,0.);
	if (M > _mPi + _mKs && M < _parameters[9]){
		const double q  = pow(utils::breakupMomentumSquared(s12, _mKs*_mKs, _mPi*_mPi), .5);
		const double q0 = pow(utils::breakupMomentumSquared(_parameters[2]*_parameters[2], _mKs*_mKs, _mPi*_mPi), .5);

		const double G = _parameters[3] * _parameters[2]/M * q/q0;
		const double deltaR = atan(_parameters[2] * G/(_parameters[2]*_parameters[2] - s12));
		const double deltaF = atan(2. * _parameters[0] * q / (2. + _parameters[0] * _parameters[1] * q * q));

		retVal += _parameters[7] * sin(deltaF + _parameters[4]) * exp(complex(0., deltaF + _parameters[4])) + 
		          _parameters[8] * sin(deltaR + _parameters[6]) * exp(complex(0., deltaR + _parameters[5])) * exp(complex(0., 2* ( deltaF + _parameters[4])));
		if (_usePS) {
			retVal *= .5*M/q;
		}

/* Copied the code form LAURA++, but this version has much less parameters... check again when Stefan is back
		const double rVal     = _parameters[1];
		const double aVal     = _parameters[0];
		const double resMass  = _parameters[2];
		const double resWidth = _parameters[3];
		const double qRatio = q/q0;
		const double totWidth = resWidth*qRatio*(resMass/M);
		const double massSqTerm = resMass*resMass - s12;
		complex altern(massSqTerm, resMass*totWidth);
		altern *= (resMass*resMass*resWidth/q0)/(massSqTerm*massSqTerm + resMass*resMass*totWidth*totWidth);

		const double tandeltaB = (2.0*aVal*q)/(2.0 + aVal*rVal*q*q);
		const double tanSq = tandeltaB*tandeltaB;
		const double cos2PhaseShift = (1.0 - tanSq) / (1.0 + tanSq);
		const double sin2PhaseShift = 2.0*tandeltaB / (1.0 + tanSq);
		complex phaseShift(cos2PhaseShift, sin2PhaseShift);

		altern *= phaseShift;
		const double qcotdeltaB = 1.0/aVal + (rVal*q*q)/2.0;
		complex bkgAmplitude(qcotdeltaB, q);
		bkgAmplitude *= M/(qcotdeltaB*qcotdeltaB + q*q);

		altern += bkgAmplitude;

		std::cout << retVal << " " << altern << " but check me again, when Stafen is here, since LAURA++ does not use phiFR, FRm" << std::endl;
*/
	}
	return retVal;
}
//------------------------------------------------------------------------------
BELLE_LASS::BELLE_LASS(const realVector& parameters, double mPi, double mK) : massShape("BELLE_LASS", parameters, {"a", "r", "M0", "G0", "phiF", "phiR", "F", "R"}), _mPi(mPi), _mK(mK) {}

complex BELLE_LASS::eval(double s12) const {

//	double pi180inv = 1.0/EvtConst::radToDegrees;
//	double ampl =_amp;
//	double theta=_phase;

//	double s = mab2;

	const double _a    = _parameters[0]; //_LASS_a;
	const double _r    = _parameters[1]; //_LASS_r;
	const double _R    = _parameters[7]; //_LASS_R;
	const double _phiR = _parameters[5]; //_LASS_phi_R;
	const double _F    = _parameters[6]; //_LASS_F;
	const double _phiF = _parameters[4]; //_LASS_phi_F;

	// T = R sin(deltaR) exp[i (deltaR+phiR)] exp[i2 (deltaB+phiB)] + B sin(deltaB+phiB) exp[i (deltaB+phiB)]
	//   = R exp[i (phiR+2phiB)] sin(deltaR) exp[i deltaR] exp[i2 deltaB] + B exp[i phiB] { cos(phiB) + cot(deltaB) sin(phiB) } sin(deltaB) exp[i deltaB]
	//   = R exp[i (phiR+2phiB) m0*Gamma(m)/(m0*m0-m*m-i*m0*Gamma(m)) exp[i2 deltaB] + B exp[i phiB] { cos(phiB) + cot(deltaB) sin(phiB) m/[q*cot(deltaB)-i*q]
	// The propagator is T/rho, where rho = 2 q/sqrt(s) is the two-body phase space

//	const double gamma = _parameters[3]; //_gammaR;
//	const double bwm   = _parameters[2]; //_massR;

	const double mR     = _parameters[2]; //_massR;
	const double gammaR = _parameters[3]; //_gammaR;

	double fR=1.0; // K*0(1430) has spin zero
	int power=1; // power is 1 for spin zero

	const double mAB=pow(s12,.5); //   (_p4_d1+_p4_d2).mass();

	const double mA=_mPi; //_p4_d1.mass();
	const double mB=_mK; //_p4_d2.mass();

	const double pAB=pow( (((mAB*mAB-mA*mA-mB*mB)*(mAB*mAB-mA*mA-mB*mB)/4.0) -
			  mA*mA*mB*mB)/(mAB*mAB),.5);
	const double q=pAB;

	const double pR=pow( (((mR*mR-mA*mA-mB*mB)*(mR*mR-mA*mA-mB*mB)/4.0) -
			  mA*mA*mB*mB)/(mR*mR),.5);

	if (std::isnan(pAB)) {
		return complex(0.,0.);
	}

	// compute running width g
	const double g = gammaR*pow(pAB/pR,power)*(mR/mAB)*fR*fR;

	const complex propagator_relativistic_BreitWigner = 1./complex(mR*mR - mAB*mAB,-mR*g);

	// non-resonant phase shift
	const double cot_deltaF = 1.0/(_a*q) + 0.5*_r*q;
	const double qcot_deltaF = 1.0/_a + 0.5*_r*q*q;

	// calculate resonant part
	const complex expi2deltaF = complex(qcot_deltaF, q)/ complex(qcot_deltaF, -q);

	const complex resonant_term_T = _R * complex(cos(_phiR + 2 * _phiF), sin(_phiR + 2 * _phiF)) * propagator_relativistic_BreitWigner * mR * gammaR * mR / pR * expi2deltaF;

	// calculate non-resonant part
	const complex  non_resonant_term_F = _F * complex(cos(_phiF), sin(_phiF)) * (cos(_phiF) + cot_deltaF * sin(_phiF)) * sqrt(s12) / complex(qcot_deltaF, -q);

	// sum up non-resonant and resonant terms
	const complex LASS_contribution = non_resonant_term_F + resonant_term_T;

//	// convert complex to TComplex
//	TComplex LASS_contribution_TComplex (LASS_contribution._rpart, LASS_contribution._ipart);

//	TComplex matrixEl = ampl * TComplex(cos(theta*pi180inv), sin(theta*pi180inv)) * LASS_contribution_TComplex;
	return LASS_contribution;
}
//------------------------------------------------------------------------------
BELLE_KMatrix::BELLE_KMatrix(const realVector& parameters) : massShape("BELLE_KMatrix", parameters, 
	{"Kmatrix_beta1_Amplitude","Kmatrix_beta1_Phase",
	 "Kmatrix_beta2_Amplitude","Kmatrix_beta2_Phase",
	 "Kmatrix_beta3_Amplitude","Kmatrix_beta3_Phase",
	 "Kmatrix_beta4_Amplitude","Kmatrix_beta4_Phase",
	 "Kmatrix_beta5_Amplitude","Kmatrix_beta5_Phase",
	 "Kmatrix_f_prod_11_Amplitude", "Kmatrix_f_prod_11_Phase",
	 "Kmatrix_f_prod_12_Amplitude", "Kmatrix_f_prod_12_Phase",
	 "Kmatrix_f_prod_13_Amplitude", "Kmatrix_f_prod_13_Phase",
	 "Kmatrix_f_prod_14_Amplitude", "Kmatrix_f_prod_14_Phase",
	 "Kmatrix_f_prod_15_Amplitude", "Kmatrix_f_prod_15_Phase",
	 "Kmatrix_s_prod_0"}) {}

complex BELLE_KMatrix::eval(double s12) const {
	double s = s12;

	//Define the complex coupling constants
	//The convention is as follow
	//i=0 --> pi+ pi-
	//i=1 --> KK
	//i=2 --> 4pi
	//i=3 --> eta eta
	//i=4 --> eta eta'
	realMatrix g(5,realVector(5,0.)); // Coupling constants. The first index is the pole index. The second index is the decay channel

	//pi+pi- channel
	g[0][0]=0.22889;
	g[1][0]=0.94128;
	g[2][0]=0.36856;
	g[3][0]=0.33650;
	g[4][0]=0.18171;
	//K+K- channel
	g[0][1]=-0.55377;
	g[1][1]=0.55095;
	g[2][1]=0.23888;
	g[3][1]=0.40907;
	g[4][1]=-0.17558;
	//4pi channel
	g[0][2]=0;
	g[1][2]=0;
	g[2][2]=0.55639;
	g[3][2]=0.85679;
	g[4][2]=-0.79658;
	//eta eta channel
	g[0][3]=-0.39899;
	g[1][3]=0.39065;
	g[2][3]=0.18340;
	g[3][3]=0.19906;
	g[4][3]=-0.00355;
	//eta eta' channel
	g[0][4]=-0.34639;
	g[1][4]=0.31503;
	g[2][4]=0.18681;
	g[3][4]=-0.00984;
	g[4][4]=0.22358;

	// Pole masses
	realVector ma(5,0.);   // Pole masses. The unit is in GeV

	ma[0]=0.651;
	ma[1]=1.20360;
	ma[2]=1.55817;
	ma[3]=1.21000;
	ma[4]=1.82206;

	// scattering data
	// double s_0_scatt=-3.92637;
	// double s_A=1.0;
	// double s_A0=-0.15;

	//Now define the K-matrix pole
	complex n11,n12,n13,n14,n15,n21,n22,n23,n24,n25,n31,n32,n33,n34,n35,n41,n42,n43,n44,n45,n51,n52,n53,n54,n55;
	double  rho1sq,rho2sq,rho4sq,rho5sq;//,rho3sq
	complex rho1,rho2,rho3,rho4,rho5;
	complexVector rho(5,complex(0.,0.));
	complex pole,SVT,Alder;
	complex det;
	complexMatrix i(5,complexVector(5, complex(0.,0.)));  //inverse of the (I-iKp) matrix
	realMatrix f(5, realVector(5,0.));

	//Initalize the mass of the resonance
	double mpi   = 0.13957;
	double mK    = 0.493677;     //using charged K value
	double meta  = 0.54775;    //using PDG value
	double metap = 0.95778;   //using PDG value

	//Initialize the matrix to value zero
	complexMatrix K(5, complexVector(5, complex(0.,0.)));
//	for(Int_t k=0;k<5;k++) {
//		for(Int_t l=0;l<5;l++) {
//			i[k][l]=complex(0,0);
//			K[k][l]=complex(0,0);
//			f[k][l]=0;
////			f_scatt[k][l]=0;
//		}
//		rho[k]=0;
//	}

	//Input the _f[i][j] scattering data
	double s_scatt=-3.92637;
	double sa=1.0;
	double sa_0=-0.15;

	f[0][0]=0.23399;  // f^scatt
	f[0][1]=0.15044;
	f[0][2]=-0.20545;
	f[0][3]=0.32825;
	f[0][4]=0.35412;

	f[1][0]=f[0][1];
	f[2][0]=f[0][2];
	f[3][0]=f[0][3];
	f[4][0]=f[0][4];

	//Construct the phase-space factor
	//For eta-eta' there is no difference term
	rho1sq=(1.0-(pow((mpi+mpi),2)/s));   //pi+ pi- phase factor
	if(rho1sq >=0) {
		rho1=complex(sqrt(rho1sq),0);
	} else {
		rho1=complex(0,sqrt(-rho1sq));
	}
	rho[0]=rho1;

	rho2sq=(1.0-(pow((mK+mK),2)/s));
	if(rho2sq >=0) {
		rho2=complex(sqrt(rho2sq),0);
	} else {
		rho2=complex(0,sqrt(-rho2sq));
	}
	rho[1]=rho2;
	//using the A&S 4pi phase space Factor:
	rho3=complex(0,0);
	//Shit, not continue
	if(s<=1) {
		double real = 1.2274+0.00370909/(s*s) - (0.111203)/(s) - 6.39017*s +16.8358*s*s - 21.8845*s*s*s + 11.3153*s*s*s*s;
		double cont32=sqrt(1.0-(16.0*mpi*mpi));
		rho3=complex(cont32*real,0);
	} else {
		rho3=complex(sqrt(1.0-(16.0*mpi*mpi/s)),0);
	}
	rho[2]=rho3;
	//
	rho4sq=(1.0-(pow((meta+meta),2)/s));
	if(rho4sq>=0) {
		rho4=complex(sqrt(rho4sq),0);
	} else {
		rho4=complex(0,sqrt(-rho4sq));
	}
	rho[3]=rho4;
	//
	rho5sq=(1.0-(pow((meta+metap),2)/s));
	if(rho5sq >=0) {
		rho5=complex(sqrt(rho5sq),0);
	} else {
		rho5=complex(0,sqrt(-rho5sq));
	}
	rho[4]=rho5;

	//sum the pole
	//equation (3) in the E791 K-matrix paper
	for(size_t k = 0; k < 5 ; ++k) {
		for(size_t l = 0; l < 5; ++l) {
			for (size_t pole_index = 0; pole_index < 5; ++pole_index) {
				double A = g[pole_index][k]*g[pole_index][l];
				double B = ma[pole_index]*ma[pole_index]-s;
				K[k][l]  = K[k][l]+complex(A/B,0);
			}
		}
	}

	//add the SVT part
	for(size_t k = 0; k < 5; ++k) {
		for(size_t l = 0; l < 5; ++l) {
			double C=f[k][l]*(1.0-s_scatt);
			double D=(s-s_scatt);
			K[k][l]=K[k][l]+complex(C/D,0);
		}
	}

	//Include the Alder zero term:
	for(size_t k = 0; k < 5; ++k) {
		for(size_t l = 0; l < 5; ++l) {
			double E=(s-(sa*mpi*mpi*0.5))*(1.0-sa_0);
			double F=(s-sa_0);
			K[k][l]=K[k][l]*complex(E/F,0);
		}
	}

	n11=complex(1,0)-complex(0,1)*K[0][0]*rho[0];
	n12=complex(0,0)-complex(0,1)*K[0][1]*rho[1];
	n13=complex(0,0)-complex(0,1)*K[0][2]*rho[2];
	n14=complex(0,0)-complex(0,1)*K[0][3]*rho[3];
	n15=complex(0,0)-complex(0,1)*K[0][4]*rho[4];

	n21=complex(0,0)-complex(0,1)*K[1][0]*rho[0];
	n22=complex(1,0)-complex(0,1)*K[1][1]*rho[1];
	n23=complex(0,0)-complex(0,1)*K[1][2]*rho[2];
	n24=complex(0,0)-complex(0,1)*K[1][3]*rho[3];
	n25=complex(0,0)-complex(0,1)*K[1][4]*rho[4];

	n31=complex(0,0)-complex(0,1)*K[2][0]*rho[0];
	n32=complex(0,0)-complex(0,1)*K[2][1]*rho[1];
	n33=complex(1,0)-complex(0,1)*K[2][2]*rho[2];
	n34=complex(0,0)-complex(0,1)*K[2][3]*rho[3];
	n35=complex(0,0)-complex(0,1)*K[2][4]*rho[4];

	n41=complex(0,0)-complex(0,1)*K[3][0]*rho[0];
	n42=complex(0,0)-complex(0,1)*K[3][1]*rho[1];
	n43=complex(0,0)-complex(0,1)*K[3][2]*rho[2];
	n44=complex(1,0)-complex(0,1)*K[3][3]*rho[3];
	n45=complex(0,0)-complex(0,1)*K[3][4]*rho[4];

	n51=complex(0,0)-complex(0,1)*K[4][0]*rho[0];
	n52=complex(0,0)-complex(0,1)*K[4][1]*rho[1];
	n53=complex(0,0)-complex(0,1)*K[4][2]*rho[2];
	n54=complex(0,0)-complex(0,1)*K[4][3]*rho[3];
	n55=complex(1,0)-complex(0,1)*K[4][4]*rho[4];

	  //Compute determinant
	det = (n15*n24*n33*n42*n51 - n14*n25*n33*n42*n51 - n15*n23*n34*n42*n51 +
		 n13*n25*n34*n42*n51 + n14*n23*n35*n42*n51 - n13*n24*n35*n42*n51 -
		 n15*n24*n32*n43*n51 + n14*n25*n32*n43*n51 + n15*n22*n34*n43*n51 -
		 n12*n25*n34*n43*n51 - n14*n22*n35*n43*n51 + n12*n24*n35*n43*n51 +
		 n15*n23*n32*n44*n51 - n13*n25*n32*n44*n51 - n15*n22*n33*n44*n51 +
		 n12*n25*n33*n44*n51 + n13*n22*n35*n44*n51 - n12*n23*n35*n44*n51 -
		 n14*n23*n32*n45*n51 + n13*n24*n32*n45*n51 + n14*n22*n33*n45*n51 -
		 n12*n24*n33*n45*n51 - n13*n22*n34*n45*n51 + n12*n23*n34*n45*n51 -
		 n15*n24*n33*n41*n52 + n14*n25*n33*n41*n52 + n15*n23*n34*n41*n52 -
		 n13*n25*n34*n41*n52 - n14*n23*n35*n41*n52 + n13*n24*n35*n41*n52 +
		 n15*n24*n31*n43*n52 - n14*n25*n31*n43*n52 - n15*n21*n34*n43*n52 +
		 n11*n25*n34*n43*n52 + n14*n21*n35*n43*n52 - n11*n24*n35*n43*n52 -
		 n15*n23*n31*n44*n52 + n13*n25*n31*n44*n52 + n15*n21*n33*n44*n52 -
		 n11*n25*n33*n44*n52 - n13*n21*n35*n44*n52 + n11*n23*n35*n44*n52 +
		 n14*n23*n31*n45*n52 - n13*n24*n31*n45*n52 - n14*n21*n33*n45*n52 +
		 n11*n24*n33*n45*n52 + n13*n21*n34*n45*n52 - n11*n23*n34*n45*n52 +
		 n15*n24*n32*n41*n53 - n14*n25*n32*n41*n53 - n15*n22*n34*n41*n53 +
		 n12*n25*n34*n41*n53 + n14*n22*n35*n41*n53 - n12*n24*n35*n41*n53 -
		 n15*n24*n31*n42*n53 + n14*n25*n31*n42*n53 + n15*n21*n34*n42*n53 -
		 n11*n25*n34*n42*n53 - n14*n21*n35*n42*n53 + n11*n24*n35*n42*n53 +
		 n15*n22*n31*n44*n53 - n12*n25*n31*n44*n53 - n15*n21*n32*n44*n53 +
		 n11*n25*n32*n44*n53 + n12*n21*n35*n44*n53 - n11*n22*n35*n44*n53 -
		 n14*n22*n31*n45*n53 + n12*n24*n31*n45*n53 + n14*n21*n32*n45*n53 -
		 n11*n24*n32*n45*n53 - n12*n21*n34*n45*n53 + n11*n22*n34*n45*n53 -
		 n15*n23*n32*n41*n54 + n13*n25*n32*n41*n54 + n15*n22*n33*n41*n54 -
		 n12*n25*n33*n41*n54 - n13*n22*n35*n41*n54 + n12*n23*n35*n41*n54 +
		 n15*n23*n31*n42*n54 - n13*n25*n31*n42*n54 - n15*n21*n33*n42*n54 +
		 n11*n25*n33*n42*n54 + n13*n21*n35*n42*n54 - n11*n23*n35*n42*n54 -
		 n15*n22*n31*n43*n54 + n12*n25*n31*n43*n54 + n15*n21*n32*n43*n54 -
		 n11*n25*n32*n43*n54 - n12*n21*n35*n43*n54 + n11*n22*n35*n43*n54 +
		 n13*n22*n31*n45*n54 - n12*n23*n31*n45*n54 - n13*n21*n32*n45*n54 +
		 n11*n23*n32*n45*n54 + n12*n21*n33*n45*n54 - n11*n22*n33*n45*n54 +
		 n14*n23*n32*n41*n55 - n13*n24*n32*n41*n55 - n14*n22*n33*n41*n55 +
		 n12*n24*n33*n41*n55 + n13*n22*n34*n41*n55 - n12*n23*n34*n41*n55 -
		 n14*n23*n31*n42*n55 + n13*n24*n31*n42*n55 + n14*n21*n33*n42*n55 -
		 n11*n24*n33*n42*n55 - n13*n21*n34*n42*n55 + n11*n23*n34*n42*n55 +
		 n14*n22*n31*n43*n55 - n12*n24*n31*n43*n55 - n14*n21*n32*n43*n55 +
		 n11*n24*n32*n43*n55 + n12*n21*n34*n43*n55 - n11*n22*n34*n43*n55 -
		 n13*n22*n31*n44*n55 + n12*n23*n31*n44*n55 + n13*n21*n32*n44*n55 -
		 n11*n23*n32*n44*n55 - n12*n21*n33*n44*n55 + n11*n22*n33*n44*n55);

	  //The 1st row of the inverse matrix. This matrix is {(I-iKp)^-1}_0j
	i[0][0]   = (n25*n34*n43*n52 -
		     n24*n35*n43*n52 - n25*n33*n44*n52 + n23*n35*n44*n52 +
		     n24*n33*n45*n52 - n23*n34*n45*n52 - n25*n34*n42*n53 +
		     n24*n35*n42*n53 + n25*n32*n44*n53 - n22*n35*n44*n53 -
		     n24*n32*n45*n53 + n22*n34*n45*n53 + n25*n33*n42*n54 -
		     n23*n35*n42*n54 - n25*n32*n43*n54 + n22*n35*n43*n54 +
		     n23*n32*n45*n54 - n22*n33*n45*n54 - n24*n33*n42*n55 +
		     n23*n34*n42*n55 + n24*n32*n43*n55 - n22*n34*n43*n55 -
		     n23*n32*n44*n55 + n22*n33*n44*n55)/det;

	i[0][1]   = (-n15*n34*n43*n52 +
		     n14*n35*n43*n52 + n15*n33*n44*n52 - n13*n35*n44*n52 -
		     n14*n33*n45*n52 + n13*n34*n45*n52 + n15*n34*n42*n53 -
		     n14*n35*n42*n53 - n15*n32*n44*n53 + n12*n35*n44*n53 +
		     n14*n32*n45*n53 - n12*n34*n45*n53 - n15*n33*n42*n54 +
		     n13*n35*n42*n54 + n15*n32*n43*n54 - n12*n35*n43*n54 -
		     n13*n32*n45*n54 + n12*n33*n45*n54 + n14*n33*n42*n55 -
		     n13*n34*n42*n55 - n14*n32*n43*n55 + n12*n34*n43*n55 +
		     n13*n32*n44*n55 - n12*n33*n44*n55)/det;

	i[0][2]   = (n15*n24*n43*n52 -
		     n14*n25*n43*n52 - n15*n23*n44*n52 + n13*n25*n44*n52 +
		     n14*n23*n45*n52 - n13*n24*n45*n52 - n15*n24*n42*n53 +
		     n14*n25*n42*n53 + n15*n22*n44*n53 - n12*n25*n44*n53 -
		     n14*n22*n45*n53 + n12*n24*n45*n53 + n15*n23*n42*n54 -
		     n13*n25*n42*n54 - n15*n22*n43*n54 + n12*n25*n43*n54 +
		     n13*n22*n45*n54 - n12*n23*n45*n54 - n14*n23*n42*n55 +
		     n13*n24*n42*n55 + n14*n22*n43*n55 - n12*n24*n43*n55 -
		     n13*n22*n44*n55 + n12*n23*n44*n55)/det;

	i[0][3]   = (-n15*n24*n33*n52 +
		     n14*n25*n33*n52 + n15*n23*n34*n52 - n13*n25*n34*n52 -
		     n14*n23*n35*n52 + n13*n24*n35*n52 + n15*n24*n32*n53 -
		     n14*n25*n32*n53 - n15*n22*n34*n53 + n12*n25*n34*n53 +
		     n14*n22*n35*n53 - n12*n24*n35*n53 - n15*n23*n32*n54 +
		     n13*n25*n32*n54 + n15*n22*n33*n54 - n12*n25*n33*n54 -
		     n13*n22*n35*n54 + n12*n23*n35*n54 + n14*n23*n32*n55 -
		     n13*n24*n32*n55 - n14*n22*n33*n55 + n12*n24*n33*n55 +
		     n13*n22*n34*n55 - n12*n23*n34*n55)/det;

	i[0][4]   = (n15*n24*n33*n42 -
		     n14*n25*n33*n42 - n15*n23*n34*n42 + n13*n25*n34*n42 +
		     n14*n23*n35*n42 - n13*n24*n35*n42 - n15*n24*n32*n43 +
		     n14*n25*n32*n43 + n15*n22*n34*n43 - n12*n25*n34*n43 -
		     n14*n22*n35*n43 + n12*n24*n35*n43 + n15*n23*n32*n44 -
		     n13*n25*n32*n44 - n15*n22*n33*n44 + n12*n25*n33*n44 +
		     n13*n22*n35*n44 - n12*n23*n35*n44 - n14*n23*n32*n45 +
		     n13*n24*n32*n45 + n14*n22*n33*n45 - n12*n24*n33*n45 -
		     n13*n22*n34*n45 + n12*n23*n34*n45)/det;

	// convert betas and f_prods to TComplex

	complexVector _beta(5, complex(0.,0.));

	_beta[0] = _parameters[0] * complex(cos(_parameters[1]),sin(_parameters[1]));
	_beta[1] = _parameters[2] * complex(cos(_parameters[3]),sin(_parameters[3]));
	_beta[2] = _parameters[4] * complex(cos(_parameters[5]),sin(_parameters[5]));
	_beta[3] = _parameters[6] * complex(cos(_parameters[7]),sin(_parameters[7]));
	_beta[4] = _parameters[8] * complex(cos(_parameters[9]),sin(_parameters[9]));

	complex _fr11prod = _parameters[10] * complex(cos(_parameters[11]),sin(_parameters[11]));
	complex _fr12prod = _parameters[12] * complex(cos(_parameters[13]),sin(_parameters[13]));
	complex _fr13prod = _parameters[14] * complex(cos(_parameters[15]),sin(_parameters[15]));
	complex _fr14prod = _parameters[16] * complex(cos(_parameters[17]),sin(_parameters[17]));
	complex _fr15prod = _parameters[18] * complex(cos(_parameters[19]),sin(_parameters[19]));

	double _s0prod = _parameters[20];

	complexVector U1j(5, complex(0.,0.));
	for(size_t j = 0; j < 5; ++j) {
		U1j[j] = i[0][j];
	}

	// compute product of inverse matrix times production vector, split production vector into two pieces
	complex value0(0.,0.);
	complex value1(0.,0.);

	  // compute inverse_matrix times first part of production vector, sum all the poles
	for(size_t l = 0; l < 5; ++l) {
		for (size_t pole_index = 0; pole_index < 5; ++pole_index) {
			complex A = _beta[pole_index]*g[pole_index][l];
			double B               = ma[pole_index]*ma[pole_index]-s;
			value0                += U1j[l] * A / B;
		}
	}

	// compute inverse_matrix times second part of production vector
	value1 += U1j[0]*_fr11prod;
	value1 += U1j[1]*_fr12prod;
	value1 += U1j[2]*_fr13prod;
	value1 += U1j[3]*_fr14prod;
	value1 += U1j[4]*_fr15prod;
	value1 *= (1-_s0prod)/(s-_s0prod);

	// Computes the final F0 vector by adding the two parts
	complex value = value0 + value1;

//	TComplex F_0_TComplex (value._rpart, value._ipart);

	return value;

}
//------------------------------------------------------------------------------
CLEO_I2::CLEO_I2(const realVector& parameters, double mPi) : massShape("CLEO_I2", parameters, 
	{"a","b","c","d","m_min", "m_max", "DeltaEta"}), _mPi(mPi) {}

complex CLEO_I2::eval(double s12) const {
	double delta_0_2 = -_parameters[0] * pow(s12/4 - _mPi*_mPi,.5)/(1. + _parameters[1] * s12 + _parameters[2]*s12*s12 + _parameters[3] * s12*s12*s12);
	complex phaseFactor(cos(2*delta_0_2), sin(2*delta_0_2));

	double eta_0_2 = 1.;
	double m12 = pow(s12, .5);
	if (m12 > _parameters[4]) {
		if (m12 < _parameters[5]) {
			eta_0_2 = 1. - _parameters[6]/2 * (1. - cos(M_PI*(m12 - _parameters[4])/(_parameters[4] - _parameters[5])));
		} else {
			eta_0_2 = 1. - _parameters[6];
		}
	}

	return (eta_0_2 * phaseFactor - 1.)/complex(0.,2.);
}
//------------------------------------------------------------------------------
CLEO_BW::CLEO_BW(const std::string& name, const realVector parameters, double mMother, double mBachelor, double mDaughter1, double mDaughter2, size_t L, double R, double r, bool exponentialFF) :
   massShape("CLEO_BW_" + name, parameters, {"m0_" + name, "G0_"+name}),
   _mMother(mMother), _mBachelor(mBachelor), _mDaughter1(mDaughter1), _mDaughter2(mDaughter2), _L(L), _R(R), _r(r), _expFF(exponentialFF) {
	if (_L > 2) {
		std::cout << "CLEO_BW::CLEO_BW: ERROR: L > 2 not implemented yet." << std::endl;
		throw;	
	}
}

complex CLEO_BW::eval(double s12) const { 

	double P = pow(utils::breakupMomentumSquared(_mMother*_mMother, s12, _mBachelor*_mBachelor), .5);
	double p = pow(utils::breakupMomentumSquared(s12, _mDaughter1*_mDaughter1, _mDaughter2*_mDaughter2), .5);

	double Pv = pow(utils::breakupMomentumSquared(_mMother*_mMother, _parameters[0] * _parameters[0], _mBachelor*_mBachelor), .5);
	double pv = pow(utils::breakupMomentumSquared(_parameters[0]*_parameters[0], _mDaughter1*_mDaughter1, _mDaughter2*_mDaughter2), .5);

	double Q  = _R*P;
	double Qv = _R*Pv;

	double q  = _r*p;
	double qv = _r*pv;

	double FF = 1.;
	double ff = 1.;
	if (_L == 0 && _expFF) {
		FF = exp(-(Q*Q - Qv *Qv)/12);
		ff = exp(-(q*q - qv *qv)/12);
	} else if (_L == 1) {
		FF = pow((1. + Qv*Qv)/(1. + Q*Q), .5);
		ff = pow((1. + qv*qv)/(1. + q*q), .5);
	} else if (_L == 2) {
		FF = pow((9. + 3*Qv*Qv + Qv*Qv*Qv*Qv)/(9. + 3*Q*Q + Q*Q*Q*Q),.5);
		ff = pow((9. + 3*qv*qv + qv*qv*qv*qv)/(9. + 3*q*q + q*q*q*q),.5);
	}
	double Gamma = _parameters[1];
	Gamma *= _parameters[0]/pow(s12,.5);
	Gamma *= pow(p/pv, 2*_L+1);
	Gamma *= ff*ff;

	complex BW_denominator(_parameters[0] * _parameters[0] - s12, - _parameters[0] * Gamma);
	return complex(FF,0.)/BW_denominator;
};

