#ifndef GENERATOR__
#define GENERATOR__
#include<string>
#include<vector>
#include<memory>
#include"kinematicSignature.h"
class generator {
	public:
		generator ();

		virtual realVector generate () const;

		size_t                              nKin        ()         const {return _kinSignature->nKin();};
		std::shared_ptr<kinematicSignature> kinSignature()         const {return _kinSignature;}
		bool                                setMaxFail  (size_t n)       {_maxFail = n; return true;};
	protected:
		std::shared_ptr<kinematicSignature> _kinSignature;
		size_t                              _maxFail;
		mutable size_t                      _failCount;
};

class threeParticleMassGenerator : public generator {
	public:
		threeParticleMassGenerator (double initialMass,const realVector& fsMasses, std::shared_ptr<kinematicSignature> kinSig = std::make_shared<kinematicSignature>(1));

		realVector generate    ()                               const override;
		bool                isValidPoint(const realVector& kin) const;
	protected:
		double              _initialMass;
		realVector _fsMasses;
		
};

class dalitzScanGenerator : public threeParticleMassGenerator {
	public: 
		dalitzScanGenerator(double            initialMass, 
		                    const realVector& fsMasses, 
		                    const realVector& massesX,
		                    const realVector& massesY,
		                    std::shared_ptr<kinematicSignature> kinSig = std::make_shared<kinematicSignature>(1));
		bool       ended() const;
		realVector getNext () const;
		realVector generate() const override;

		void reset() {indexX = 0; indexY = 0;}
	protected:
		mutable size_t indexX;
		mutable size_t indexY;

		realVector massesX;
		realVector massesY;

};

#endif//GENERATOR__
