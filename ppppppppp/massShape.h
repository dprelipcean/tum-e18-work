#ifndef MASSSHAPE__
#define MASSSHAPE__
#include<string>
#include<complex>
#include<vector>
#include"types.h"
class massShape {
	public:
		massShape (std::string name, realVector pars, std::vector<std::string> parNames = std::vector<std::string>());

		virtual complex eval (double s) const;

		bool                         setParameter (size_t n, double val);
		std::pair<bool, double>      getParameter (size_t n)                   const;
		bool                         setParName   (size_t n, std::string name);
		std::pair<bool, std::string> getParName   (size_t n)                   const;

		std::string                  name         () const {return _name;}

	protected:
		size_t                   _nPar;
		std::string              _name;
		realVector      _parameters;
		std::vector<std::string> _parameterNames;
};

class simpleBW : public massShape {
	public:
		simpleBW (double mass, double width);

		complex eval (double s) const override;
};

class stepLike : public massShape {
	public:
		stepLike (double sMin, double sMax);

		complex eval (double s) const override;
};

class sawtooth : public massShape {
	public:
		sawtooth (double sMin, double sMax);

		complex eval (double s) const override;
};

class spike : public massShape {
	public:
		spike ( double sMin, double sEval, double sMax);

		complex eval (double s) const override;
};

class spline_deg : public massShape {
	public:
		spline_deg (double sMin, double sMax, size_t deg);

		complex eval(double s) const override;

		double Deval(double s, size_t nth) const;

	private:
		size_t _deg;
};

class constant : public massShape {
	public:
		constant ();

		complex eval (double s) const override;
};

class zeroMode0pp : public massShape {
	public:
		zeroMode0pp (double s, double m2);

		complex eval (double s) const override;
};

class zeroMode1mm : public massShape {
	public:
		zeroMode1mm (double s, double m2);

		complex eval (double s) const override;
};

class polynomialMassShape : public massShape {
	public:
		polynomialMassShape(complexVector coefficients, double baseExponent = 1.);

		complex eval (double s) const override;
	protected:
		size_t _polDeg;
		double _baseExponent;	
};

class BELLEbreitWigner : public massShape {
	public:
		BELLEbreitWigner(std::string name, double mass, double width, size_t spin, double motherMass, double bachelorMass, double daughterMass1, double daughterMass2);

		complex eval (double s) const override;
	protected:
		size_t _spin;

		double _motherMass;
		double _bachelorMass;
		double _daughterMass1;
		double _daughterMass2;

		double _Rr;		
		double _RD;
};

class BELLE_LASS_KpiS : public massShape {
	public:
		BELLE_LASS_KpiS(const realVector& parameters, double mPi, double mKs, bool usePS = false);
		
		complex eval(double s12) const override;
	protected:
		bool   _usePS;
		double _mPi;
		double _mKs;
	
};

class BELLE_LASS : public massShape {
	public:
		BELLE_LASS(const realVector& parameters, double mPi, double mK);

		complex eval(double s12) const override;

	protected:
		double _mPi;
		double _mK;

};

class BELLE_KMatrix : public massShape {
	public:
		BELLE_KMatrix(const realVector& parameters);

		complex eval(double s12) const override;
};

class CLEO_I2 : public massShape { // (arXiv:0802.4214v2 [hep-ex] 12 Sep 2008)
	public:
		CLEO_I2(const realVector& parameters, double mPi);

		complex eval(double s12) const override;
	protected:
		double _mPi;

};

class CLEO_BW : public massShape { // (arXiv:0802.4214v2 [hep-ex] 12 Sep 2008)
	public:
		CLEO_BW(const std::string& name, const realVector parameters, double mMother, double mBachelor, double mDaughter1, double mDaughter2, size_t L, double R, double r, bool exponentialFF);
		complex eval(double s12) const override;

	protected:
		double _mMother;
		double _mBachelor;
		double _mDaughter1;
		double _mDaughter2;
		size_t _L;

		double _R;
		double _r;

		bool _expFF;
};
#endif//MASSSHAPE__
