#ifndef ROOT_FIT_WRAPPER__
#define ROOT_FIT_WRAPPER__
#include<memory>
#include<complex>
#include<string>
#include"logLikelihood.h"

class ROOT_fit_wrapper {
	public:
		ROOT_fit_wrapper(std::shared_ptr<logLikelihoodBase> ll);

		static double DoEval(const double* xx);
//		static realVector DoGradient(const double* xx);

		std::pair<double, complexVector > fit(const complexVector& startVals);
		std::pair<double, complexVector > fit(const realVector& startVal);

		static std::shared_ptr<logLikelihoodBase> __ll;
	protected:
		size_t                             _minimizerStrategy;
		double                             _initStepSize;
		std::shared_ptr<logLikelihoodBase> _ll;
};
#endif//ROOT_FIT_WRAPPER__
