#ifndef UTILS_GENERAL___LARENEG_SLITU
#define UTILS_GENERAL___LARENEG_SLITU

namespace utils {
	inline void makeError(const std::string& functionName, const std::string& errorMessage) {
		std::cout << "\033[1;31m" << functionName << ":\033[0m " << errorMessage << std::endl;
	}

	inline void makeWarning(const std::string& functionName, const std::string& warnMessage) {
		std::cout << "\033[1;33m" << functionName << ":\033[0m " << warnMessage << std::endl;
	}

	inline void makeInfo(const std::string& functionName, const std::string& infoMessage) {
		std::cout << "\033[1;32m" << functionName << ":\033[0m " << infoMessage << std::endl;
	}

	void opening() {
		std::cout << "This is ppppppppp, a\n Program for\n Pseudoscalar to\n Pseudoscalar,\n Pseudoscalar,\n Pseudoscalar\n Processes for\n Pushing our\n Paper\n Publication" << std::endl;
	}

	double random() {
		double retVal = ((double) rand() / (RAND_MAX));
		return retVal;
	}

	double random2() { 
		return 2*random() -1.;	
	}

	double breakupMomentumSquared(const double s, const double s1, const double s2) {
		return (s*s + s1*s1 + s2*s2 - 2*(s*s1 + s*s2 + s1*s2))/4/s;
	}

	double degToRad(double deg) {
		return deg * M_PI / 180.;
	}

	std::string to_string(const complex& val) {
		return "(" + std::to_string(val.real()) + "," + std::to_string(val.imag()) + ")";
	}

	bool checkComplexDouble() {
		complexVector vals = {complex(1.,2.), complex(3.,4.)};
		double* val = (double*)&vals[0];
		bool works = true;
		if (val[0] != 1.) {
			works = false;
		} 
		if (val[1] != 2.) {
			works = false;
		}
		if (val[2] != 3.) {
			works = false;
		}
		if (val[3] != 4.) {
			works = false;
		}
		if (!works) {
			makeError("utils_general::checkComplexDouble()","Complex double may not be reinterpreted as re, im array.");
		}
		return works;
	}

	std::vector<std::string> splitString(const std::string& inputString, const char& splitChar) {
		std::vector<std::string> retVal;
		std::string part = "";
		for (const char& c : inputString) {
			if (c == splitChar) {
				if (part != "") {
					retVal.push_back(part);
					part = "";
				}
			} else {
				part += c;
			}
		}
		if (part != "") {
			retVal.push_back(part);
		}
		return retVal;
	}

	bool stringStartsWith(const std::string& inString, const std::string& startString) {
		if (inString.size() < startString.size()) {
			return false;
		}
		for (size_t i = 0; i < startString.size(); ++i) {
			if (inString[i] != startString[i]) {
				return false;
			}
		}
		return true;
	}

}
#endif//UTILS_GENERAL___LARENEG_SLITU

