#ifndef LOGLIKELIHOOD__
#define LOGLIKELIHOOD__
#include<vector>
#include<memory>
#include<complex>
#include<iostream>

#include"amplitude.h"
#include"integrator.h"
#include"kinematicSignature.h"

#include"utils_general.h" // Include this explicitely, since otherwise it would be a circiular include (Instead of "utils.h").

class logLikelihoodBase {
	public:
		logLikelihoodBase (size_t nAmpl, std::shared_ptr<kinematicSignature> kinSignature, std::shared_ptr<integrator> integral);
	
		std::pair<double, complexVector > fitNlopt   (const realVector& parameters);
		std::pair<double, complexVector > fitNlopt   (const complexVector& parameters);
		virtual double                    nloptCall  (const realVector &x, realVector &grad) const;

		virtual double     eval       (const complexVector& prodAmps) const 
                                                                       {utils::makeError("logLikelihoodBase::eval(...)",
		                                                        "Base class eval called (" + std::to_string(prodAmps.size()) +" parameters), returning 0.");
		                                                        return 0.;}

		virtual realVector Deval      (const complexVector& prodAmps) const
		                                                       {utils::makeError("logLikelihoodBase::Deval(...)",
		                                                        "Base class Deval called (" + std::to_string(prodAmps.size()) +" parameters), returning empty vector.");
		                                                        return realVector();}

		virtual realMatrix DDeval     (const complexVector& prodAmps) const
		                                                       {utils::makeError("logLikelihoodBase::DDeval(...)",
		                                                        "Base class DDeval called (" + std::to_string(prodAmps.size()) +" parameters), returning empty vector.");
		                                                        return realMatrix();}

		virtual bool       loadDataPoints (const realMatrix& dataPoints, size_t maxNperEvent = 9) 
                                                                       {utils::makeError("logLikelihoodBase::loadDataPoints(...)","No loadDataPoints(...) method specified. Do not load " +
		                                                        std::to_string(dataPoints.size()) + " points (maxNperEvent = " +std::to_string(maxNperEvent) + ").");
		                                                        return false;}

		realVector                makeFinalGradient(const realVector& params, const realVector& fullGradient) const;
		realMatrix                makeFinalHessian(const realVector& params, const realMatrix& fullHessian) const;

		realVector                getFullParameters(const realVector& params) const;
		realVector                cutGradient( const realVector& params) const;

		complexVector fullParamsToProdAmps(const realVector& params                 ) const;
		realVector                prodAmpsToFullParams(const complexVector& prodAmps) const;


		realMatrix  DDconstrainedProdAmps(const realVector& params) const;
		realMatrix  DprodAmpsNumerical(const realVector& params, double delta = 1.e-5) const;

		bool                    setNcallsPrint    (size_t nCallPrint)        {_nCallsPrint = nCallPrint; return true;}
		bool                    setExtended       (bool flag)                {_extended = flag; return true;}
		bool                    fixParameter      (size_t nPar, double value);
		bool                    addCopyParameter  (size_t nPar, size_t copyFrom, int nScale = -1, double scaler = 1.);

		size_t                  getNparTot        ()                   const;
		size_t                  getNpar           ()                   const;
		virtual size_t          getNparFit        ()                   const {return getNpar();}
		size_t                  nAmpl             ()                   const {return _nAmpl;}
		size_t                  nCalls            ()                   const {return _nCalls;}
		bool                    resetNcalls       ()                         {_nCalls = 0; return true;}
		bool                    trackFit          (bool flag)                {_trackFit = flag; return true;}

		const std::vector<std::pair<size_t,double> >              getFixedParameters() const {return _fixedParameters;}
		const std::vector<std::tuple<size_t,size_t,int, double> > getCopyParameters () const {return _copyParameters ;}

		virtual complexVector                makeProdAmpsFromFitPars(const realVector& fitPars) const 
		                                                                                 {return fullParamsToProdAmps(getFullParameters(fitPars));}

		realMatrix getStoreParams()   const {return _storeParams;}
		realVector getFitParameters() const {return _directFitParameters;}
		bool       setFitParameters(const realVector& param) {_directFitParameters = param; return param.size() == getNparFit();}
	protected:
		bool                                                _extended;
		std::shared_ptr<kinematicSignature>                 _kinSignature;
		mutable size_t                                      _nCalls;
		size_t                                              _nCallsPrint;
		size_t                                              _nAmpl;
		size_t                                              _nPoints;
		size_t                                              _nScale; // Number of scale parameters
		double                                              _numDelta;
		std::shared_ptr<integrator>                         _integral;
		std::vector<std::pair<size_t,double> >              _fixedParameters;
		std::vector<std::tuple<size_t,size_t,int, double> > _copyParameters;
		realVector                                          _directFitParameters; // Stores the actual real parameters, that are adjusted by the fitter.
// //
		bool                                                _trackFit;
		mutable realMatrix                                  _storeParams;
};

class logLikelihood : public logLikelihoodBase {
	public:
		logLikelihood (std::vector<std::shared_ptr<amplitude> > amplitudes, std::shared_ptr<integrator> integral);

		virtual double     eval       (const complexVector& prodAmps)            const override;
		virtual realVector Deval      (const complexVector& prodAmps)            const override;
		virtual realMatrix DDeval     (const complexVector& prodAmps)            const override;

		virtual bool                       loadDataPoints (const realMatrix& dataPoints, size_t maxNperEvent = 9) override;

		size_t                             getSector(size_t a) const;
		bool                               setCoherenceBorders(sizeVector borders);

		bool                               setMaxII(double maxII) {_maximumIncoherentIntensity = maxII; return true;}
	protected:
		size_t                                   _nSect;
		double                                   _maximumIncoherentIntensity;
		std::vector<std::shared_ptr<amplitude> > _amplitudes;
		complexMatrix _points;
		sizeVector                               _amplitudeCoherenceBorders;
		std::vector<sizeVector >                 _contributingWaves;

};

class linearLogLikelihood : public logLikelihood {
	public:
		linearLogLikelihood(std::vector<std::shared_ptr<amplitude> > amplitudes, std::shared_ptr<integrator> integral, realMatrix linearTransformation);

		virtual double nloptCall (const realVector &x, realVector &grad) const override;

		size_t dimFull()      const {return _dimFull;}
		size_t dimTransform() const {return _dimTransform;}

		virtual size_t        getNparFit() const override {return _dimTransform;}
		virtual complexVector makeProdAmpsFromFitPars(const realVector& fitPars) const override;

		realVector getFullRealParametersLinear(const realVector& fitPars) const;

		double     linearEval  (const realVector& params) const;
		realVector linearDeval (const realVector& params) const;
		realMatrix linearDDeval(const realVector& params) const;

	protected:
		size_t     _dimFull;
		size_t     _dimTransform;
		realMatrix _transformation;

};

class logLikelihood_withPrior : public logLikelihood {
	public: 
		logLikelihood_withPrior (std::vector<std::shared_ptr<amplitude> > amplitudes, std::shared_ptr<integrator> integral);

		virtual double     eval       (const complexVector& prodAmps) const override;
		virtual realVector Deval      (const complexVector& prodAmps) const override;
		virtual realMatrix DDeval     (const complexVector& prodAmps) const override;

		bool addPriorDirection(double strength, const complexVector direction);

		bool setInterferencePriorStrength(double strength);

		double                    interferencePriorFunc(double coherent, double incoherent) const;
		std::pair<double, double> DinterferencePriorFunc(double coherent, double incoherent) const;
		realVector                DDinterferencePriorFunc(double coherent, double incoherent) const;

	protected:
		bool          _noNormWarn;
		double        _interferencePriorStrength;

		size_t        _nPrior;
		realVector    _priorStrengths;
		complexMatrix _priorDirections;

};

class logLikelihoodAllFree : public logLikelihoodBase {
	public:
		logLikelihoodAllFree (realVector binning, std::vector<std::shared_ptr<angularDependence> > freedAmplitudes, std::shared_ptr<integrator> integral);

		double     eval       (const complexVector& prodAmps)            const override;
		realVector Deval      (const complexVector& prodAmps)            const override;
		realMatrix DDeval     (const complexVector& prodAmps)            const override;

		bool                                        loadDataPoints (const realMatrix& dataPoints, size_t maxNperEvent = 9) override;
		std::pair<bool, std::pair<size_t, size_t> > findBin(const realVector& point) const;

	protected:
		size_t                     _nBins;
		realVector                 _binning;
		std::vector<sizeVector>    _eventsPerBin;
		std::vector<complexMatrix> _amplitudesInBin;
};
#endif//LOGLIKELIHOOD__
