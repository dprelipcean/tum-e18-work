#pragma once
#include"utils.h"
#include"massShape.h"

class splineMatcher {
	public:
		splineMatcher(const std::vector<std::shared_ptr<spline_deg> >& functionsFirstBin);

		static realMatrix getM(const std::vector<std::shared_ptr<spline_deg> >& functions);

		realMatrix getLastN() const;

		realMatrix getLastValueMatrix() const;

		realMatrix getTransformationMatrix() const {return _transformationMatrix;}

		bool addBin(const std::vector<std::shared_ptr<spline_deg> >& newFunctions);

	protected:
		size_t _nBin;   // Number of bins
		size_t _degree; // Highest derivative, that is still continuous

		std::vector<std::shared_ptr<spline_deg> > _functions; // Spline functions

		realMatrix _transformationMatrix; // The transformation matrix from fit parameter to coefficients for the single functions
		// first index: _degree^th derivative of the first bin low edge, _degree-1^th..., ampl value at fisrt bin low edge,
		//              ampl value at second bin low edge = ampl valueat first bin upper edge, ..., ampl value at last bin upper edge
		// second index: coefficient of first function of the first bin, second function first bin, ..., last function last bin
};
