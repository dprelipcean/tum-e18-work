#ifndef FFT___FFT
#define  FFT___FFT
#include<vector>


namespace FFT {
	template<typename T>
	complexVector FFT(const std::vector<T>& vals) {
		const size_t dim = vals.size();
		const double norm = std::pow((double) dim, -.5); // (*) Use a symmetric definition. has to be consistent with (**)
		complexVector retVal(dim, complex(0.,0.));
		for (size_t k = 0; k < dim; ++k) {
			for (size_t n = 0; n < dim; ++n) {
				retVal[k] += vals[n] * norm * std::exp(complex(0.,-2.*M_PI*k*n/dim));
			}
		}
		return retVal;
	}

	complexVector FFT_inv(const complexVector& vals) {
		const size_t dim = vals.size();
		const double norm = std::pow((double) dim, -.5); // (**) Use a symmetric definition. has to be consistent with (*)
		complexVector retVal(dim, complex(0.,0.));
		for (size_t k = 0; k < dim; ++k) {
			for (size_t n = 0; n < dim; ++n) {
				retVal[n] += vals[k] * norm * std::exp(complex(0.,2.*M_PI*k*n/dim));
			}
		}
		return retVal;
	}
}
#endif//FFT___FFT
