#ifndef PHASESPACEFUNCTIONS__
#define PHASESPACEFUNCTIONS__
#include<iostream>
#include<string>
#include<vector>

#include"types.h"

class realPhaseSpaceBase {
	public:
		realPhaseSpaceBase();
		virtual double eval(double s) const {std::cout << "realPhaseSpaceBase::eval(...): WARNING: Base class method called with s = " << s << std::endl; return 0.;}
};

class twoBodyPhaseSpace : public realPhaseSpaceBase {
	public:
		twoBodyPhaseSpace(double m1, double m2);

		double eval(double s) const override;
	protected:
		double _m1;
		double _m2;
		double _sTh;
};

class fourPiPhaseSpace : public realPhaseSpaceBase{
	public:
		fourPiPhaseSpace(const std::string valueFileName, size_t nStep, double sMin, double sMax, double mPi, double mRho, double gamma, double sanityDelta = 1.e-6);

		double eval(double s) const override;
		double interpolate_rho51(double s) const;
		double rho52(double s) const;

		bool sanityCheck(double sanityDelta) const;
	protected:
		size_t _nStep;
		double _sMin;
		double _sMax;
		double _step;
		double _mPi;
		double _mRho;
		double _gamma;

		std::string _valueFileName;
		realVector  _values;
};
#endif//PHASESPACEFUNCTIONS__
