#ifndef UTILS_FIT___TIF_SLITU
#define UTILS_FIT___TIF_SLITU
namespace utils {
	realVector makeSplineLocalFunction(size_t degree, const std::vector<std::shared_ptr<spline_deg> >& functions) {
// Creates the coefficients for the functions to give a d-times continuous differentiable
// Function, that is f^(i) = 0 at the beginning and the end
// Degree is the last continuous derivative (degree = 1 => f' is continuous)
		const size_t dim = (degree+2)*(degree+2);
		if (functions.size() != dim) {
			makeError("utils_fit::makeSplineLocalFunction(...)","Number iof functions does not match.");
			return realVector();
		}
		realMatrix theMatrix(dim, realVector(dim, 0)); // First index is the number of the constraint, second is the function index

		for (size_t c = 0; c < degree + 1; ++c ){
			// Initial (0,0,0,0,...) constraint
			// c is the constraint number and the derivative degree at the same time
			for (size_t f = 0; f < degree + 2; ++f) {
				// Start by evaluating the lower border (_parameters[0]) of the first function)
				theMatrix[f][c] = functions[f]->Deval(functions[f]->getParameter(0).second, c);
			}
		}
		for (size_t c = 0; c < degree + 1; ++c) {
			// (*) Constraint to force the total linear combination away from zero
			// Doint this here might not be the optimal choice, since the resultiung
			// functions tend to peak in the center of the bin range, bit it is the 
			// most simple and general way to go here
			for (size_t f =0; f < degree+2; ++f) {
				theMatrix[f][degree+1] = functions[f]->Deval(functions[f]->getParameter(1).second, 0);
			}
		}
		for (size_t m = 0; m < degree + 1; ++m) {
			// The matching constraints of bin m with bin m+1
			size_t d = 0; // Derivative index (has to start at zero, while the constraint index is set off.
			for (size_t c = (degree+1)*(m+1) + 1; c < (degree+1)*(m+2) + 1; ++c) {
				// One constrint per matching and derivative to be continuous
				for (size_t f = m*(degree+2); f < (m+1)*(degree+2); ++f) { // Function index for the functions in the lower bin
					size_t F = f + degree + 2; // Function index in the upper bin
					theMatrix[f][c] = functions[f]->Deval(functions[f]->getParameter(1).second, d); // Eval lower bin at upper border
					theMatrix[F][c] =-functions[F]->Deval(functions[F]->getParameter(0).second, d); // Eval upper bin at lower border (relative minus sign (sum has to be zero))
				}
				++d; // Increase the dericative for the next bin
			}
		}
		size_t d = 0; // Again us a different derivative index
		for (size_t c = (degree+1)*(degree+2)+1; c < dim; ++c) {
			for (size_t f = (degree+1)*(degree+2);f < dim; ++f) {
				theMatrix[f][c] = functions[f]->Deval(functions[f]->getParameter(1).second, d); // Evaluate at upper bin border
			}
			++d;
		}

		realMatrix inverted = inv(theMatrix);

		realVector theVector(dim,0.);
		theVector[degree+1] = 1.; // See (*)

		realVector theParams = dot(theVector, inverted);
		return theParams;
	}

	complexVector fft_filter(const complexVector& inValues, 
	                         const realVector&    norms, 
	                         const sizeVector&    nWaves,
	                         double               filterCoefficient,
	                         size_t               skipBatchesBelow) {
// This method Fourier transforms the start values of every freed wave (with a number of waves >= skipBatchesBelow).
// Then it removes all but the first 100*filterCoefficient% of the Fourier coefficients (Sets them to zero)
// In the end it transforms the remainung coefficients back. Ths gives an equally distributed, but smooth 
// set of starting values (they qre equally ditributed, but correllated). The normalizations are properly taken into account,
// To smoothe in amplitude level (as opposed to number of event level)
		const size_t dim = inValues.size();
		if (norms.size() != dim) {
			makeError("utils_fit::fft_filter(...)", "Dimension of normalization does not match.");
			throw;
		}
		complexVector normedValues(dim);
		for (size_t a = 0; a < dim; ++a) {
			normedValues[a] = inValues[a] * norms[a];
		}
		size_t count = 0;
		complexVector filtered(dim);
		for (size_t b : nWaves) {
			if (b < skipBatchesBelow) {
				for (size_t a = 0; a < b; ++a) {
					filtered[count + a] = normedValues[count + a];
				}
			} else {
				complexVector fft = FFT::FFT<complex>(complexVector(&normedValues[count], &normedValues[count + b]));
				size_t startZero = (size_t)(filterCoefficient * b);
				for (size_t z = startZero; z < b; ++z) {
					fft[z] = complex(0.,0.);
				}
				complexVector inverse = FFT::FFT_inv(fft);
				for (size_t a = 0; a < b; ++a) {
					filtered[count + a] = inverse[a];
				}
			}
			count += b;
		}
		if (count != dim) {
			makeError("utils_fit::fft_filter(...)", "Batch sizes don't add up to dimension of amplitudes (" + std::to_string(count) + " vs. " +std::to_string(dim) + ").");
			throw;
		}
		complexVector retVal(dim);
		for (size_t a = 0; a < dim; ++a) {
			if (norms[a] == 0.) {
				retVal[a] = complex(0.,0.);
			} else {
				retVal[a] = filtered[a] / norms[a];
			}
		}
		return retVal;
	}

	complexMatrix get_fft_prior_modes(const realVector& norms, 
	                                  const sizeVector& nWaves,
	                                  double            startPriorCoefficient,
	                                  size_t            skipBatchesBelow) {
// Give the modes, that do the fourier projections, to be used in priors.
// Take into account the normalizations properly.
		const size_t dim = norms.size();
		complexMatrix retVal;
		size_t countWaves = 0;
		for (size_t b : nWaves) {
			if (b < skipBatchesBelow) {
				countWaves += b;
				continue;
			}
			size_t startPrior = (size_t)(startPriorCoefficient * b);
			const double norm = std::pow((double) b, -.5);
			for (size_t m = startPrior; m < b; ++m) {
				retVal.push_back(complexVector(dim, complex(0.,0.)));
				for (size_t a = 0; a < b; ++a) {
					retVal.back()[countWaves+a] = norms[countWaves+a] * norm * std::exp(complex(0.,-2.*M_PI*a*m/b));
				}
			}		
			countWaves += b;
		}
		return retVal;
	}

	complexMatrix evaluateBinnedPolynomials(size_t nBins, const complexMatrix& coefficients, double xMin, double xMax); // Declare here.
	complexMatrix getChebyshevCoefficients(size_t n); 

	complexMatrix splineMatch(const realVector& norms, const sizeVector& nWaves) {
// Splines together the constant parts in the first half of norms to the linear parts in the second half
		const size_t dim = norms.size()/2;
		size_t count = 0;
		for (size_t nWave : nWaves) {
			count += nWave;
		}
		if (dim != count) {
			makeError("utils_fit::splineMatch(...)","Dimension ("+std::to_string(dim)+") does not match nWaves ("+std::to_string(count) +"=.");
		}
		complexMatrix retVal(2*dim, complexVector(dim + nWaves.size(), complex(0.,0.)));
		size_t countIn  = 0;
		size_t countOut = 0;
		for (size_t w = 0; w < nWaves.size(); ++w) {
			for ( size_t n = 0; n < nWaves[w]; ++n) {
				retVal[countOut][countIn]   = .5/norms[countOut];
				retVal[countOut][countIn+1] = .5/norms[countOut];

				retVal[dim+countOut][countIn]   = -.1/norms[countOut]; //x - .5
				retVal[dim+countOut][countIn+1] = .1/norms[countOut];
			}
		}
		return retVal;

	}

	complexMatrix get_smooth_transform(const realVector&  norms,
	                                   const sizeVector&  nWaves,
	                                   double             startPriorCoefficient,
	                                   size_t             skipBatchesBelow,
	                                   const std::string& smoothingMode) {
		makeInfo("utils_fit::get_smooth_transform(...)","Using smoothing mode '" + smoothingMode + "'.");
		size_t dim_in  = 0;
		size_t dim_out = 0;
		for (const size_t b : nWaves) {
			dim_in += b;
			if (b < skipBatchesBelow) {
				dim_out += b;
			} else {
				size_t stopModesAbove = (size_t)(startPriorCoefficient * b);
				dim_out += stopModesAbove;
			}
		}
		complexMatrix retVal(dim_in, complexVector(dim_out, complex(0.,0.)));
		size_t countWavesOut = 0;
		size_t countWavesIn  = 0;
		for (const size_t b : nWaves) {
			if (b < skipBatchesBelow) {
				for (size_t i = 0; i < b; ++i) {
					retVal[countWavesIn+i][countWavesOut+i] = complex(1.,0.);
				}
				countWavesIn  += b;
				countWavesOut += b;
				continue;
			}
			size_t stopModesAbove = (size_t)(startPriorCoefficient * b);
			complexMatrix basis(stopModesAbove, complexVector(b, complex(0.,0.)));
			if (smoothingMode == "fft") { // Fast fourier transform
				const double norm = std::pow((double) b, -.5);
				for (size_t a = 0; a < b; ++a) {
					for (size_t m = 0; m <  stopModesAbove; ++m) {
						basis[m][a] = norm * std::exp(complex(0.,2.*M_PI*a*m/b));
					}
				}
			} else if (smoothingMode == "hffft") { // Fast fourier transform of the double spectrum to avoid rapid change at 2pi = 0. (Half frequency fft)
				const double norm = std::pow((double) b, -.5); // Normalization for this method is not clear (normalize to double, or single spectrum)
				for (size_t a = 0; a < b; ++a) {
					for (size_t m = 0; m < stopModesAbove; ++m) {
						basis[m][a] = norm * std::exp(complex(0.,M_PI*a*m/b)); // The factor 2 is missing w.r.t. simple fft
					}
				}
			} else if (smoothingMode == "jumpfft") { // Simple fft with an additional mode to allow for a jump at 2pi = 0
				const double step = 2./b;
				const double norm = std::pow((double) b, -.5);
				for (size_t a = 0; a < b; ++a) {
					basis[0][a] = step * (a + .5) - 1.; // Goes from -1. to 1. as a goes from 0 to b
					for (size_t m = 0; m < stopModesAbove - 1; ++m) { // m-1, since one mode is the linear one
						basis[m+1][a] = norm * std::exp(complex(0.,2.*M_PI*a*m/b));
					}
				}
			} else if (smoothingMode == "cheby") { // Chebyshev polynomials
				basis = utils::evaluateBinnedPolynomials(b, utils::getChebyshevCoefficients(stopModesAbove-1), -1.,1.); 
				// Use stopModesAbove-1, since a polynomial of degree N, has N+1 degrees of freedom (0th order).
			} else {
				makeError("utils_fit::get_smooth_transform(...)","Unknown mode '" + smoothingMode + "' (Implemented: 'fft', 'hffft', 'jumpfft', 'cheby').");
				throw;
			}
			for (size_t a = 0; a < b; ++a) {
				for (size_t m = 0; m < stopModesAbove; ++m) {
					retVal[countWavesIn + a][countWavesOut + m] = basis[m][a]/norms[countWavesIn+a];
				}
			}
			countWavesIn  += b;
			countWavesOut += stopModesAbove;
		}
		return retVal;
	}

	realMatrix complex_to_real_matrix(const complexMatrix& c_matrix) {
//
// Transforms the complex matrix to a real valued one, so that M_c * v_c => M_r * v_r ( = sum_j M_r[i][j] v_r[j]), where v_r[2i] = Re(v_c[i]) and v_r[2i+1] = Im(v_c[i])
//
		const size_t dim_x = c_matrix.size();
		if (dim_x == 0) {
			makeWarning("utils_fit::complex_to_real_matrix(...)","Matrix of size zero received.");
			return realMatrix();
		}
		const size_t dim_y = c_matrix[0].size();
		for (const complexVector& c_line : c_matrix) {
			if (c_line.size() != dim_y) {
				makeWarning("utils_fit::complex_to_real_matrix(...)","Matrix with different sized lines encountered. Abort.");
				throw;
			}
		}
		realMatrix r_matrix = realMatrix(2*dim_x, realVector(2*dim_y, 0.));
		for (size_t i_x = 0; i_x < dim_x; ++i_x) {
			for (size_t i_y = 0; i_y < dim_y; ++i_y) {
				r_matrix[2*i_x  ][2*i_y  ] = c_matrix[i_x][i_y].real();
				r_matrix[2*i_x  ][2*i_y+1] =-c_matrix[i_x][i_y].imag();
				r_matrix[2*i_x+1][2*i_y  ] = c_matrix[i_x][i_y].imag();
				r_matrix[2*i_x+1][2*i_y+1] = c_matrix[i_x][i_y].real();
			}
		}
		return r_matrix;
	}

	bool updateBestLLfile(const std::string& llFileName, double newValue, const std::string& additionalInfo = "") { // Returns "true", <value> is better than the value in <llFileName> or the file doe not yet exist
		double val;
		bool updateFile = false;
		{
			std::ifstream fin(llFileName.c_str(), std::ifstream::in);
			if (fin >> val) {
				if (newValue < val) {
					updateFile = true;
				}
			} else {
				updateFile = true;
			}
			fin.close();
		}
		if (updateFile) {
			std::ofstream outFile;
			outFile.open(llFileName.c_str());
			outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1);
			outFile << newValue << std::endl;
			outFile << additionalInfo;
			outFile.close();
		}
		return updateFile;
	}

	template<typename T>
	void checkDerivatives(std::shared_ptr<T> ll, const complexVector& ampl, double delta, bool doSecond) {
		complexVector pa = ampl;
		const double evl = ll->eval(pa);
		const realVector Devl = ll->Deval(pa);
		for (size_t a = 0; a < pa.size(); ++a) {
			double dd = delta * pow(std::norm(pa[a]),.5);
			pa[a] += complex(dd,0.);
			std::cout << Devl[2*a  ] << " || " << (ll->eval(pa) - evl)/dd << std::endl;
			pa[a] += complex(-dd,dd);
			std::cout << Devl[2*a+1] << " || " << (ll->eval(pa) - evl)/dd << std::endl;
			pa[a] += complex(0.,-dd);
		}
		if (doSecond) {
			const realMatrix DDevl = ll->DDeval(pa);
			makeInfo("utils_fit::checkDerivatives(...)","Second derivative gotten. Checking numerically now...");
			realVector D;
			for (size_t a = 0; a < pa.size(); ++a) {
				double dd = delta * pow(std::norm(pa[a]),.5);
				pa[a] += complex(dd,0.);
				D = ll->Deval(pa);
				for (size_t b = 0; b < pa.size(); ++b) {
					if (a == b) {
						std::cout << " ! - ! - ! - ! - ! - ! - !" << std::endl;
					}
					std::cout << DDevl[2*a  ][2*b  ] << " ++ " << (D[2*b  ] - Devl[2*b  ])/dd << std::endl;
					std::cout << DDevl[2*a  ][2*b+1] << " ++ " << (D[2*b+1] - Devl[2*b+1])/dd << std::endl;
					if (a == b) {
						std::cout << " ! - ! - ! - ! - ! - ! - !" << std::endl;
					}
				}
				pa[a] += complex(-dd,dd);
				D = ll->Deval(pa);
				for (size_t b = 0; b < pa.size(); ++b) {
					if (a == b) {
						std::cout << " ! - ! - ! - ! - ! - ! - !" << std::endl;
					}
					std::cout << DDevl[2*a+1][2*b  ] << " ++ " << (D[2*b  ] - Devl[2*b  ])/dd << std::endl;
					std::cout << DDevl[2*a+1][2*b+1] << " ++ " << (D[2*b+1] - Devl[2*b+1])/dd << std::endl;
					if (a == b) {
						std::cout << " ! - ! - ! - ! - ! - ! - !" << std::endl;
					}
				}
				pa[a] += complex(0.,-dd);
			}
		}
	}

	complexMatrix evaluateBinnedPolynomials(size_t nBins, const complexMatrix& coefficients, double xMin, double xMax) {
		const size_t nPoly = coefficients.size();
		const double stepSize = (xMax - xMin)/nBins;
		complexMatrix retVal(nPoly, complexVector(nBins, complex(0.,0.)));
		for (size_t i = 0; i < nBins; ++i) {
			const double x = (i + .5) * stepSize + xMin;
			for (size_t p = 0; p < nPoly; ++p) {
				double x_to_n = 1.;
				for (size_t d = 0; d < coefficients[p].size(); ++d) {
					retVal[p][i] += coefficients[p][d] * x_to_n;
					x_to_n *= x;
				}
			}
		}
		return retVal;
	}

	complexMatrix getChebyshevCoefficients(size_t n) { // Complex for later purposes to be able to also include complex bases
		complexMatrix retVal(n+1);
		complexVector coefficientsNMinus2(1, complex(1,0.));
		retVal[0] = coefficientsNMinus2;
		if (n == 0) {
			
			return retVal;
		}
		complexVector coefficientsNMinus1(1, complex(0.,0.));
		coefficientsNMinus1.push_back(complex(1.,0.));
		retVal[1] = coefficientsNMinus1;
		if (n == 1) {
			return retVal;
		}
		complexVector coefficientsN;
		for (size_t degree = 2; degree < n+1; ++degree) {
			coefficientsN = complexVector(degree+1,complex(0.,0.)); // A polynomial of degree n has n+1 coefficients
			coefficientsN[0] = -coefficientsNMinus2[0];
			for (size_t i = 1; i < degree-1; ++i) {
				coefficientsN[i] = 2.*coefficientsNMinus1[i-1] - coefficientsNMinus2[i];
			}
			coefficientsN[degree-1] = 2.*coefficientsNMinus1[degree-2];
			coefficientsN[degree  ] = 2.*coefficientsNMinus1[degree-1];
			coefficientsNMinus2 = coefficientsNMinus1;
			coefficientsNMinus1 = coefficientsN;
			retVal[degree] = coefficientsN;
		}
		return retVal;
	}


// Use checkGradient and checkHessian with:
//	std::function<double(const realVector&)>       e = std::bind(&T::eval,   t, std::placeholders::_1);
//	std::function<realVector(const realVector&)>   d = std::bind(&T::Deval,  t, std::placeholders::_1);
//	std::function<realMatriox(const realVector&)> dd = std::bind(&T::DDeval, t, std::placeholders::_1);
// Where eval, Deval and DDeval are the function, gradient and hessian, and std::shared_ptr<T> t

	void checkGradient(double delta, const realVector& par, std::function<double(const realVector&)>& e, std::function<realVector(const realVector&)>& d) {
		const size_t dim        = par.size();
		realVector params       = par;
		const double centerEval = e(par);
		const realVector grad   = d(par);
		double maxDev           = 0.;
		for (size_t p = 0; p < dim; ++p) {
			params[p] += delta;
			double num = (e(params) - centerEval)/delta;
			maxDev = std::max(maxDev, fabs(num - grad[p]));
			std::cout << grad[p] << " || " << num << std::endl;
			params[p] -= delta;
		}
		makeInfo("utils_fit::checkGradient(...)","Maximum deviation is " + std::to_string(maxDev) + ".");
	}

	void checkHessian(double delta, realVector& par, std::function<realVector(const realVector&)>& d, std::function<realMatrix(const realVector&)> dd) {
		const size_t dim            = par.size();
		realVector params           = par;
		const realVector centerGrad = d(par);
		const realMatrix hessian    = dd(par);
		double maxDev               = 0.;
		for (size_t p = 0; p < dim; ++p) {
			params[p] += delta;
			const realVector grad = d(params);
			for (size_t q = 0; q < dim; ++q) {
				double num = (grad[q] - centerGrad[q])/delta;
				maxDev = std::max(maxDev, fabs(num - hessian[p][q]));
				std::cout << hessian[p][q] << " || " << num << std::endl;
			}
			params[p] -= delta;
			std::cout << "-----------------------------------" << std::endl;
		}
		makeInfo("utils_fit::checkHessian(...)","Maximum deviation is " + std::to_string(maxDev) + ".");
	}


	template<typename T>
	void checkNloptDerivatives(std::shared_ptr<T> ll, const realVector& params, double delta) {
		const size_t dim = params.size();
		realVector par = params;
		realVector grad(dim);
		realVector  empty;
		const double centerEval = ll->nloptCall(par, grad);
		makeInfo("utils_fit::checkNloptDerivatives(...)"," Compare numerical derivatives:");
		for (size_t p = 0; p < dim; ++p) {
			par[p] += delta;
			double eval = ll->nloptCall(par, empty);
			std::cout << "  " << p << " " << (eval - centerEval)/delta << " " << grad[p] << std::endl;
			par[p] -= delta;
		}
	}

	complexVector get_coherent_mass_shape(const realVector& s, const std::vector<std::shared_ptr<massShape> >& ampls, const complexVector& prod_amps_times_norm) {
		if (ampls.size() == 0) {
			makeError("utils_fit::get_coherent_mass_shape(...)","No amplitudes given.");
			throw;
		}
		if (ampls.size() != prod_amps_times_norm.size()) {
			makeError("utils_fit::get_coherent_mass_shape(...)","Amplitude size mismatch (" + std::to_string(ampls.size()) + " vs. " + std::to_string(prod_amps_times_norm.size()) + ").");
			throw;
		}
		complexVector retVal(s.size(), complex(0.,0.));
		for (size_t i = 0; i < s.size(); ++i) {
			for (size_t a = 0; a < ampls.size(); ++a) {
				retVal[i] += prod_amps_times_norm[a] * ampls[a]->eval(s[i]);
			}
		}
		return retVal;
	}

}
#endif//UTILS_FIT___TIF_SLITU
