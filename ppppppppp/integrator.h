#ifndef INTEGRATOR__
#define INTEGRATOR__
#include<string>
#include<vector>
#include<memory>
#include<complex>

#include"types.h"
#include"amplitude.h"
#include"generator.h"
#include"kinematicSignature.h"
#include"efficiencyFunction.h"
class integrator {
	public:
		integrator(size_t integralPoints, std::shared_ptr<generator> pointGenerator, const std::vector<std::shared_ptr<amplitude> >& amplitudes, std::shared_ptr<efficiencyFunction>& efficiency);

		bool                        integrate         ();
		bool                        loadIntegrals     (const std::string& psFileName, const std::string& accFileName);
		bool                        setIntegrals      (const complexMatrix& ps_integral, const complexMatrix& ac_integral);
		complexMatrix               getIntegralMatrix (bool accCorr = false)                     const;
		std::pair<bool, complex>    element           (size_t i, size_t j, bool accCorr = false) const;
		std::pair<bool, realVector> getNormalizations (bool accCorr = false) const;

		double     totalIntensity   (const complexVector& prodAmpl, bool accCorr = false) const;
		realVector DtotalIntensity  (const complexVector& prodAmpl, bool accCorr = false) const;
		realMatrix DDtotalIntensity (const complexVector& prodAmpl, bool accCorr = false) const;

		bool                                addIncoherentSector(std::shared_ptr<integrator> sector);

		bool                                isIntegrated ()                     const {return _isIntegrated;}
		bool                                setNpoints   (size_t n);
		bool                                writeToFile  (const std::string& fileName, bool accCorr = false) const;
		bool                                makeIntegrals(const std::string& ps_fileName, const std::string& ac_fileName, const std::string& ID_fileName);
		size_t                              nAmpl        ()                     const {return _nAmpl;}
		std::pair<bool, size_t>             getNpoints   ()                     const; 
		std::pair<bool, std::string>        getWaveName  (size_t i)             const;
		std::shared_ptr<kinematicSignature> kinSignature ()                     const {return _kinSignature;}
		sizeVector                 getCoherenceBorders()               const {return _amplitudeCoherenceBorders;}

		bool                                writeIDfile(const std::string& outFileName) const;
		bool                                checkIDfile(const std::string& inFileName)  const;

		bool                                setCoherenceBorders(sizeVector& borders);
		bool                                setNumLim(double val) {_numLim = val; return true;}
	protected:
		bool                                makeRealMatrices();

		bool                                     _isIntegrated;
		double                                   _numLim; //  Used in check of hermitianity and ID files
		std::shared_ptr<kinematicSignature>      _kinSignature;
		size_t                                   _nAmpl;
		size_t                                   _nPoints;
		sizeVector                               _amplitudeCoherenceBorders;
		std::vector<std::shared_ptr<amplitude> > _amplitudes;
		std::shared_ptr<generator>               _generator;
		std::shared_ptr<efficiencyFunction>      _efficiency;
		complexMatrix                            _integralMatrix;
		complexMatrix                            _accCorrIntegralMatrix;
		realMatrix                               _realIntegralMatrix;
		realMatrix                               _realAccCorrMatrix;
};
#endif//INTEGRATOR__

