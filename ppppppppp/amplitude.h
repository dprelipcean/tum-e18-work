#ifndef AMPLITUDE__
#define AMPLITUDE__
#include<string>
#include<complex>
#include<vector>
#include<memory>
#include"massShape.h"
#include"kinematicSignature.h"
#include"angularDependence.h"
#include"efficiencyFunction.h"

class amplitude {
	public:
		amplitude (std::shared_ptr<kinematicSignature> kinSignature, std::string name);

		virtual complex eval (const realVector& kin) const;
		double  intens (const realVector& kin) const {return std::norm(eval(kin));}

		virtual const std::shared_ptr<massShape> getMassShape () const {return nullptr;}

		std::string                         name         () const {return _name;}
		size_t                              nKin         () const {return _kinSignature->nKin();}
		std::shared_ptr<kinematicSignature> kinSignature () const {return _kinSignature;}


	protected:
		std::shared_ptr<kinematicSignature> _kinSignature;
		std::string                         _name;

};

class constantAmplitude : public amplitude {
	public:
		constantAmplitude(std::shared_ptr<kinematicSignature> kinSignature);
		complex eval(const realVector& kin) const override;

};

class threeParticleIsobaricAmplitude : public amplitude {
	public:
		threeParticleIsobaricAmplitude(bool boseSymmetrize, std::string name, std::shared_ptr<massShape> shape, std::shared_ptr<angularDependence> angDep);

		const std::shared_ptr<massShape> getMassShape () const override {return _massShape;}
		
		complex eval(const realVector& kin) const override;
		complex evalSingleSymmTerm(const realVector& kin) const;

	private:
		bool                               _bose;
		std::shared_ptr<massShape>         _massShape;
		std::shared_ptr<angularDependence> _angDep;
};

class threeParticleIsobaricAmplitudeNoBose : public amplitude {
	public:
		threeParticleIsobaricAmplitudeNoBose(size_t isobarIndex, std::string name, std::shared_ptr<massShape> shape, std::shared_ptr<angularDependence> angDep, realVector fsMasses);

		const std::shared_ptr<massShape> getMassShape () const override {return _massShape;}

		complex eval(const realVector& kin) const override;
	private:
		size_t                             _isobarIndex;
		double                             _sumFSmasses;
		std::shared_ptr<massShape>         _massShape;
		std::shared_ptr<angularDependence> _angDep;
	
};

class dalitzMonomialAmplitude : public amplitude {
	public:
		dalitzMonomialAmplitude(std::shared_ptr<kinematicSignature> kinSignature, double exponent1, double exponent2);

		complex eval(const realVector& kin) const override;
	private:
		double _exponent1;
		double _exponent2;

};

class dalitzPolynomialAmplitude : public amplitude {
	public:
		dalitzPolynomialAmplitude(std::shared_ptr<kinematicSignature> kinSignature, const std::string configurationFile, double xMin, double xMax, double yMin, double yMax);

		complex eval(const realVector& kin) const override;
		std::pair<double, double> getXY(const realVector& kin) const;

	private:
		size_t              _nTerms;
		double              _xMin;
		double              _xMax;
		double              _xWidth;
		double              _yMin;
		double              _yMax;
		double              _yWidth;
		sizeVector _xExponents;
		sizeVector _yExponents;
		realVector _coefficients;
};

class mCosTintensPolynomial : public amplitude {
	public:
		mCosTintensPolynomial(std::shared_ptr<kinematicSignature> kinSignature, const std::string configurationFile, double motherMass, realVector fsMasses, size_t isobarCombination = 12);

		complex eval(const realVector& kin) const override;
		std::pair<double, double> getMcosT(const realVector& kin) const;
		std::pair<double, double> getXY(const std::pair<double, double> mCosT) const;

		bool setCosTLimits(const std::pair<double,double> newLimits);
		bool setMlimits(const std::pair<double,double> newLimits);
	private:
		size_t                   _isobarCombination;
		size_t                   _nTerms;
		double                   _mWidth;
		double                   _cosTwidth;
		std::pair<double,double> _mLimits;
		std::pair<double,double> _cosTlimits;
		sizeVector      _xExponents;
		sizeVector      _yExponents;
		realVector      _coefficients;
		realVector      _fsMasses;
};

class lookupAmplitudeIntens : public amplitude {
	public:
		lookupAmplitudeIntens(std::shared_ptr<kinematicSignature> kinSignature, const std::string& name, double sMinX, double widthX, double sMinY, double widthY, const realMatrix & intensities);

		virtual complex eval(const realVector& kin) const override;
	protected:
		size_t _nX;
		size_t _nY;

		double _sMinX;
		double _widthX;

		double _sMinY;
		double _widthY;

		realMatrix _data;
};

class lookupAmplitudeIntens_efficiency : public lookupAmplitudeIntens {
	public:
		lookupAmplitudeIntens_efficiency(std::shared_ptr<efficiencyFunction> efficiency, std::shared_ptr<kinematicSignature> kinSignature, const std::string& name, double sMinX, double widthX, double sMinY, double widthY, const realMatrix & intensities);

		virtual complex eval(const realVector& kin) const override;
	protected:
		std::shared_ptr<efficiencyFunction> _efficiency;
};
#endif//AMPLITUDE__
