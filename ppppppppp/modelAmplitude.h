#ifndef MODELAMPLITUDE__
#define MODELAMPLITUDE__
#include<vector>
#include<complex>
#include<string>
#include<memory>
#include"amplitude.h"
#include"kinematicSignature.h"
class modelAmplitude : public amplitude {
	public:
		modelAmplitude (complexVector transitionAmplitudes, std::vector<std::shared_ptr<amplitude> > amplitudes, realVector normalizations, std::string name);
		complex eval   (const realVector& kin) const override;

		bool                                setTransitionAmplitude (size_t n, complex amp);
	protected:
		size_t                                   _nAmpl;
		complexVector       _transitionAmplitudes;
		std::vector<std::shared_ptr<amplitude> > _amplitudes;
		realVector                      _normalizations;
};
#endif// MODELAMPLITUDE__
