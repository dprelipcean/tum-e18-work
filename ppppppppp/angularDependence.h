#ifndef ANGULARDEPENDENCE__
#define ANGULARDEPENDENCE__
#include<string>
#include<vector>
#include<complex>
#include<memory>
#include"kinematicSignature.h"

class angularDependence {
	public:
		angularDependence (std::shared_ptr<kinematicSignature> kinSignature, std::string name);

		virtual complex eval (const realVector& kin) const;

		std::string                         name()         const {return _name;}	
		size_t                              nKin()         const {return _kinSignature->nKin();}
		std::shared_ptr<kinematicSignature> kinSignature() const {return _kinSignature;}
	protected:
		std::shared_ptr<kinematicSignature> _kinSignature;
		std::string                         _name;
};

class sameMassZeroS : public angularDependence {
	public:
		sameMassZeroS (double fsMass);

		complex eval (const realVector& kin) const override;
	private:
		double _fsMass;
};

class sameMassOneP : public angularDependence {
	public:
		sameMassOneP (double fsMass);

		complex eval (const realVector& kin) const override;
	private:
		double _fsMass;
};

class sameMassTwoD : public angularDependence {
	public:
		sameMassTwoD (double fsMass);

		complex eval (const realVector& kin) const override;
	private:
		double _fsMass;
};

class sameMassZeroSnonRelativistic : public angularDependence {
	public:
		sameMassZeroSnonRelativistic (double fsMass);

		complex eval(const realVector& kin) const override;
	private:
		double _fsMass;
};

class sameMassOnePnonRelativistic : public angularDependence {
	public:
		sameMassOnePnonRelativistic (double fsMass);

		complex eval(const realVector& kin) const override;
	private:
		double _fsMass;
};

class sameMassTwoDnonRelativistic : public angularDependence {
	public:
		sameMassTwoDnonRelativistic (double fsMass);

		complex eval(const realVector& kin) const override;
	private:
		double _fsMass;
};

class arbitraryMass_S_nonRelativistic : public angularDependence {
	public:
		arbitraryMass_S_nonRelativistic(size_t isobarIndex, realVector fsMasses);

		complex eval(const realVector& kin) const override;

	private:
		size_t              _isobarIndex;
		realVector _fsMasses;
};

class arbitraryMass_P_nonRelativistic : public angularDependence {
	public:
		arbitraryMass_P_nonRelativistic(size_t isobarIndex, realVector fsMasses);

		complex eval(const realVector& kin) const override;

	private:
		size_t              _isobarIndex;
		realVector _fsMasses;
};

class ratioOfDependences : public angularDependence {
	public:
		ratioOfDependences (std::shared_ptr<angularDependence> numerator, std::shared_ptr<angularDependence> denominator);

		complex eval(const realVector& kin) const override;
	private:
		std::shared_ptr<angularDependence> _numerator;
		std::shared_ptr<angularDependence> _denominator;
};

class BELLE_S : public angularDependence {
	public: 
		BELLE_S(size_t isobarIndex, const realVector& fsMasses,  std::shared_ptr<kinematicSignature> kinSignature = std::make_shared<kinematicSignature>(2));
		complex eval(const realVector& kin) const override;
	private:
		size_t _isobarIndex;
};

class BELLE_P : public angularDependence {
	public: 
		BELLE_P(size_t isobarIndex, const realVector& fsMasses,  std::shared_ptr<kinematicSignature> kinSignature = std::make_shared<kinematicSignature>(2));
		complex eval(const realVector& kin) const override;

		bool setFSmasses(const realVector newMasses);
	private:
		size_t _isobarIndex;
		realVector _fsMasses;
};

class BELLE_D : public angularDependence {
	public: 
		BELLE_D(size_t isobarIndex, const realVector& fsMasses, std::shared_ptr<kinematicSignature> kinSignature = std::make_shared<kinematicSignature>(2));
		complex eval(const realVector& kin) const override;

		bool setFSmasses(const realVector newMasses);
	private:
		size_t _isobarIndex;
		realVector _fsMasses;
};

class binnedDalitzAmplitude : public angularDependence {
	public:
		binnedDalitzAmplitude(std::shared_ptr<angularDependence> angAmp, const realVector& binningXaxis, const realVector& binningYaxis);

		complex eval(const realVector& kin) const override;

		std::shared_ptr<angularDependence> getUnbinnedAmplitude() {return _angAmp;}

		realVector getBinnedKin(const realVector& kin) const;
	private:
		std::shared_ptr<angularDependence> _angAmp;
		realVector _binningX;
		realVector _binningY;
};
#endif// ANGULARDEPENDENCE__
