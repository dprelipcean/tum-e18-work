#include "ROOT_fit_wrapper.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include<string>

std::shared_ptr<logLikelihoodBase> ROOT_fit_wrapper::__ll = std::shared_ptr<logLikelihoodBase>(nullptr);

ROOT_fit_wrapper::ROOT_fit_wrapper(std::shared_ptr<logLikelihoodBase> ll) : _minimizerStrategy(0), _initStepSize(1.), _ll(ll) {} 

double ROOT_fit_wrapper::DoEval(const double* xx) {
	realVector grad(0);
	realVector vectorized(xx, xx + __ll->getNparFit());
	double retVal = __ll->nloptCall(vectorized, grad);
	return retVal;
}

//realVector ROOT_fit_wrapper::DoGradient(const double* xx) {
//	realVector grad(__ll->getNpar(),0.);
//	__ll->nloptCall(realVector(xx, xx + __ll->getNpar()), grad);
//	return grad;
//}

std::pair<double, complexVector > ROOT_fit_wrapper::fit(const complexVector& startVals) {
	realVector startPars = _ll->cutGradient(_ll->prodAmpsToFullParams(startVals)); // Use this, since no better option now... (will change the start parameters a bit)
	return fit(startPars);
}

std::pair<double, complexVector > ROOT_fit_wrapper::fit(const realVector& startPars) {
	__ll = _ll;
	ROOT::Math::Functor f(&ROOT_fit_wrapper::DoEval, __ll->getNparFit());
	std::shared_ptr<ROOT::Math::Minimizer> min(ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad"));
	min->SetStrategy(_minimizerStrategy);
	min->SetFunction(f);
	min->SetMaxIterations(   1000000);
	min->SetMaxFunctionCalls(1000000);
	for (size_t p = 0; p < _ll->getNparFit(); ++p) {
		std::string parName = "par_" + std::to_string(p);
		double stepSize = _initStepSize;
		min->SetVariable(p, parName, startPars[p], stepSize);
	}
	min->Minimize();
	const realVector result(min->X(),min->X() + _ll->getNparFit());
	complexVector finalProdAmps = _ll->makeProdAmpsFromFitPars(result);
	_ll->setFitParameters(result);
	return std::pair<double, complexVector > (_ll->eval(finalProdAmps), finalProdAmps);
}
