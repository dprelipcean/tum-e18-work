#ifndef UTILS_RW___WR_SLITU
#define UTILS_RW___WR_SLITU
namespace utils {
	template<typename T> 
	std::vector<std::vector<T> > reshape(const std::vector<T>& inVector, size_t dimX, size_t dimY) {
		std::vector<std::vector<T> > retVal(dimY, std::vector<T>(dimX, 0.));
		size_t countX = 0;
		size_t countY = 0;
		for (const T& val : inVector) {
			retVal[countY][countX] = val;
			++countX;
			if (countX == dimX) {
				countX = 0;
				++countY;
			}
		}
		if (countX != 0) {
			makeError("utils_rw::reshape(...)"," not finish at line end.");
			throw;
		}
		if (countY != dimY) {
			makeError("utils_rw::reshape(...)","Did not find specified number of lines.");
			throw;
		}
		return retVal;	
	}

	template<typename T, typename U>
	std::pair<std::vector<T>, std::vector<U> > readTwoTypeValuesFromTextFile(const std::string& inFileName, bool allowZeroLength = false) {
		std::vector<T> retVal1;
		std::vector<U> retVal2;
		T v;
		U w;
		std::ifstream fin(inFileName.c_str(), std::ifstream::in);
		while (true) {
			if (fin>>v) {
				retVal1.push_back(v);
			} else {
				break;
			}
			if (fin>>w) {
				retVal2.push_back(w);
			} else {
				makeError("utils_rw::readTwoTypeValuesFromTextFile(...)","Could parse only one value.");
				throw;
			}
		}
		if (retVal1.size() == 0 and !allowZeroLength) {
			makeError("utils_rw::readTwoTypeValuesFromTextFile(...)","No values could be read from '" + inFileName + "'.");
			throw;
		}
		fin.close();
		return std::pair<std::vector<T>, std::vector<U> >(retVal1, retVal2);
	}

	template<typename T>
	std::vector<T> readValuesFromTextFile(const std::string& inFileName, bool allowZeroLength = false) {
		std::vector<T> retVal;
		T v;
		std::ifstream fin(inFileName.c_str(), std::ifstream::in);
		while (fin>>v) {
			retVal.push_back(v);
		}
		if (retVal.size() == 0 and !allowZeroLength) {
			makeError("utils_rw::readValuesFromTextFile(...)","No values could be read from '" + inFileName + "'.");
			throw;
		}
		fin.close();
		return retVal;
	}

	template<typename T>
	std::pair<bool, std::vector<std::vector<T> > > readMatrixFromTextFile(const std::string& inFileName, size_t dim) {
		std::vector<std::vector<T> > retVal(dim, std::vector<T>(dim));
		size_t col = 0;
		size_t lin = 0;
		T c;
		std::ifstream fin(inFileName.c_str(), std::ifstream::in);
		while (fin>>c) {
			if (lin >= dim || col >= dim) {
				makeError("utils_rw::readMatrixFromTextFile(...)","Number of input values exceeds matrix dimensions.");
				return std::pair<bool, std::vector<std::vector<T> > >(false, std::vector<std::vector<T> >());
			}

			retVal[lin][col] = c;
			col += 1;
			if (col == dim) {
				col  = 0;
				lin += 1;
			}
		}
		if (col != 0 || lin != dim) {
			makeError("utils_rw::readMatrixFromTextFile(...)","Mismatch in matrix dimensions: ("+std::to_string(col)+","+std::to_string(lin)+
			                                               ")= (col,lin) != (0,dim) = (0,"+std::to_string(dim)+").");
			return std::pair<bool, std::vector<std::vector<T> > >(false, std::vector<std::vector<T> >());
		}
		fin.close();
		return std::pair<bool, std::vector<std::vector<T> > >(true, retVal);
	}


	template<typename T>
	bool mergeMatrixFile(const std::string& inFileName, const std::string& outFileName,const std::vector<sizeVector >& mergeMap) {
// // // // // Set up & read
		size_t inDim = 0;
		size_t countChunks = 0;
		sizeVector linChunks;
		for (const sizeVector& mergeList : mergeMap) {
			for(const size_t chunk : mergeList) {
				inDim += chunk;
				++countChunks;
				linChunks.push_back(chunk);
			}
		}
		std::pair<bool, std::vector<std::vector<T> > > inMatrix = readMatrixFromTextFile<T>(inFileName, inDim);

// // // // // Merge
		std::vector<std::vector<T> > outMatrix(countChunks, std::vector<T>(countChunks, 0.) );
		size_t count1 = 0;
		for (size_t c1 = 0; c1 < countChunks; ++c1) {
			size_t chunk1 = linChunks[c1];
			for (size_t i1 = 0; i1 < chunk1; ++i1) {
				size_t count2 = 0;
				for (size_t c2 = 0; c2 < countChunks; ++c2) {
					size_t chunk2 = linChunks[c2];
					for (size_t i2 = 0; i2 < chunk2; ++i2) {
						outMatrix[c1][c2] += inMatrix.second[count1 + i1][count2 + i2];
					}
					count2 += chunk2;
				}
			}
			count1 += chunk1;
		}
// // // // // Write
		std::ofstream outFile;
		outFile.open (outFileName.c_str());
		outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1);
		for (size_t i = 0; i < countChunks; ++i) {
			for (size_t j = 0; j < countChunks; ++j) {
				outFile << outMatrix[i][j] << " ";
			}
			outFile << std::endl;
		}
		outFile.close();
		return true;
	}
}
#endif//UTILS_RW___WR_SLITU
