#ifndef UTILS_PHYS___SYHP_SLITU
#define UTILS_PHYS___SYHP_SLITU
namespace utils {

	bool isSanePoint(const realVector& point, const realVector& fsMasses) {
		const double masssum = fsMasses[0]*fsMasses[0] + fsMasses[1]*fsMasses[1] + fsMasses[2]*fsMasses[2];

		const double s1 = fsMasses[0]*fsMasses[0];
		const double s2 = fsMasses[1]*fsMasses[1];
		const double s3 = fsMasses[2]*fsMasses[2];

		const double& s12 = point[1];
		const double& s13 = point[2];
		const double  s23 = point[0] + masssum - point[1] - point[2];
		
		const double& s = point[0];

		const double Q2_12 = breakupMomentumSquared(s,s12, s3);
		if (Q2_12 < 0.) {
			return false;
		}
		const double q2_12 = breakupMomentumSquared(s12, s1, s2);
		if (q2_12 < 0.) {
			return false;
		}

		const double Q2_13 = breakupMomentumSquared(s,s13, s2);
		if (Q2_13 < 0.) {
			return false;
		}
		const double q2_13 = breakupMomentumSquared(s13, s1, s3);
		if (q2_13 < 0.) {
			return false;
		}

		const double Q2_23 = breakupMomentumSquared(s, s23, s1);
		if (Q2_23 < 0.) {
			return false;
		}
		const double q2_23 = breakupMomentumSquared(s23, s2, s3);
		if (q2_23 < 0.) {
			return false;
		}

		const double p1p2 = (s12 - s1 - s2)/2;
		const double p1p3 = (s13 - s1 - s3)/2;
//		const double p2p3 = (s23 - s2 - s3)/2;

		// Do s13 system
		{
			double E1 = sqrt(q2_13 + s1);
//			p1 + p3 = (sqrt(s13), 0, 0, 0) p2 = (E2, p2x, p2y, p2z) => (p1+p3)p2 = sqrt(s13) E2 = 2 E1 E2 // In the specific rest frame
//			p1p2 = (s12 - m1^2 - m2^2)/2; p2p3 - (s23 - m2^2 - m3^2)/2 => (p1+p3)p2 = (s12 + s23 - m1^2 - 2m2^2 - m3^2)/2 = (s - s13 - m2^2)/2 // Lorentz-invariant
//			=> E2 = (s - s13 - m2^2)/(4 E1)
			double E2 = (s - s13 - s2)/sqrt(s13)/2;

			double p1 = sqrt(q2_13);
			double p2 = sqrt(E2*E2 - s2);

			double cosT = (E1*E2 - p1p2)/p1/p2; // In The 13-isobar rest frame.

			if (cosT < -1. or cosT > 1.) {
				return false;
			}

		}{
		// Do s12 rest system
			double E1 = sqrt(q2_12 + s1);
			double E3 = (s - s12 - s3)/sqrt(s12)/2;

			double p1 = sqrt(q2_12);
			double p3 = sqrt(E3*E3 - s3);

			double cosT = (E1*E3 - p1p3)/p1/p3;

			if (cosT < -1. or cosT > 1.) {
				return false;
			}
		}{ 
		// Do the 23 system 
			double E1 = (s - s23 - s1)/sqrt(s23)/2;
			double E2 = sqrt(q2_23 + s2);
	
			double p1 = sqrt(E1*E1 - s1);
			double p2 = sqrt(q2_23);
	
			double cosT = (E1*E2 - p1p2)/p1/p2;

			if (cosT < -1. or cosT > 1.) {
				return false;
			}
		}

		return true;
	}

	realMatrix sanitizeBELLEdataPoints(const realMatrix& inData, const realVector& fsMasses) {
		if (fsMasses.size() != 3) {
			makeError("utils_phys::sanitizeBELLEdataPoints(...)","fsMasses has to have size 3.");
			throw;
		}
		const size_t nIn = inData.size();
		if (nIn == 0) {
			return realMatrix();
		}
		const size_t dim = inData[0].size();
		size_t count = 0;
		realMatrix retVal(nIn, realVector(dim, 0.));
		for (const realVector & point : inData) {
			if (isSanePoint(point, fsMasses)) {
				retVal[count][0] = point[0];
				retVal[count][1] = point[1];
				retVal[count][2] = point[2];
				++count;
			}
		}
		retVal.resize(count);
		makeInfo("utils_phys::sanitizeDataPoints(...)",std::to_string(count) + " events of " + std::to_string(nIn) + " events are sane.");

		return retVal;
	}

}
#endif//UTILS_PHYS___SYHP_SLITU
