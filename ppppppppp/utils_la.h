#pragma once
#include<Eigen/Dense>
namespace utils {
	realMatrix inv(const realMatrix& matrix) {
		const size_t dim = matrix.size();
		for (const realVector& line : matrix) {
			if (line.size() != dim) {
				makeError("utils_la::inv(realMatrix)","Not a square matrix.");
				return realMatrix();
			}
		}
		Eigen::MatrixXd mat(dim,dim);
		for (size_t i = 0; i < dim; ++i) {
			for (size_t j = 0; j < dim; ++j) {
				mat(i,j) = matrix[i][j];
			}
		}
		mat = mat.inverse();
		realMatrix retVal(dim, realVector(dim, 0.));
		for (size_t i = 0; i < dim; ++i) {
			for (size_t j = 0; j < dim; ++j) {
				retVal[i][j] = mat(i,j);
			}
		}
		return retVal;
	}

	realVector dot(const realMatrix& M, const realVector& v) {
		const size_t dim_cnt = v.size(); // Contracted dimension
		const size_t dim_ret = M.size(); // Return vector dimension
		for (size_t i = 0; i < dim_ret; ++i) {
			if (M[i].size() != dim_cnt) {
				makeError("utils_la::dot(realMatrix, realVector)","Line " + std::to_string(i) + " of the matrix has the wrong size: " + std::to_string(M[i].size()) + " instead of " + std::to_string(dim_cnt) + " = v.size().");
				return realVector();
			}
		}
		realVector retVal(dim_ret, 0.);
		for (size_t i = 0; i < dim_ret; ++i) {
			for (size_t j = 0; j < dim_cnt; ++j) {
				retVal[i] += M[i][j] * v[j];
			}
		}
		return retVal;
	}

	realVector dot(const realVector& v, const realMatrix& M) {
		const size_t dim_cnt = v.size(); // Contracted dimension
		if (M.size() != dim_cnt) {
			makeError("utils_la::dot(realVector, realMatrix)","Matrix has the wrong size: " + std::to_string(M.size()) + " instead of " + std::to_string(dim_cnt) + " = v.size().");
			return realVector();
		}
		if (dim_cnt == 0) {
			return realVector(); // Just to make it work. Sizes have to match anyway
		}
		const size_t dim_ret = M[0].size(); // Return vector dimension
		for (size_t i = 0; i < dim_cnt; ++i) {
			if (M[i].size() != dim_ret) {
				makeError("utils_la::dot(realMatrix, realVector)","Line " + std::to_string(i) + " of the matrix has the wrong size: " + std::to_string(M[i].size()) + " instead of " + std::to_string(dim_ret) + " = M[0].size().");
				return realVector();
			}
		}
		realVector retVal (dim_ret, 0.);
		for (size_t i = 0; i < dim_cnt; ++i) {
			for (size_t j = 0; j < dim_ret; ++j) {
				retVal[j] += v[i]*M[i][j];
			}
		}
		return retVal;
	}

	realMatrix unit(size_t dim) {
		realMatrix retVal(dim, realVector(dim, 0.));
		for (size_t i = 0; i < dim ; ++i) {
			retVal[i][i] = 1.;
		}
		return retVal;
	}

	realMatrix mergeBlockMatrix(const std::vector<realMatrix>& matrices) {
		size_t dim_i = 0;
		size_t dim_j = 0;
		for (const realMatrix& matrix : matrices) {
			const size_t dim_i_mat = matrix.size();
			if (dim_i_mat == 0) {
					makeError("untils_la::mergeBlockMatrix(matrices)","Encountered zero dimension in i, cannot determine correspondoing j dimension.");
					return realMatrix();
			}
			const size_t dim_j_mat = matrix[0].size();
			for (const realVector& line : matrix) {
				if (line.size() != dim_j_mat) {
					makeError("untils_la::mergeBlockMatrix(matrices)","Received matrix lines of varying sizes.");
					return realMatrix();
				}
			}
			dim_i += dim_i_mat;
			dim_j += dim_j_mat;
		}
		realMatrix retVal(dim_i, realVector(dim_j, 0.));
		size_t count_i = 0;
		size_t count_j = 0;
		for (const realMatrix& matrix : matrices) {
			const size_t dim_i_mat = matrix.size();
			const size_t dim_j_mat = matrix[0].size();
			for (size_t i = 0; i < dim_i_mat;++i) {for (size_t j = 0; j < dim_j_mat; ++j) {
				retVal[count_i+i][count_j+j] = matrix[i][j];
			}}
			count_i += dim_i_mat;
			count_j += dim_j_mat;
		}
		return retVal;
	}

	template <typename T> 
	bool normLinearTransform(std::vector<std::vector<T> >& matrix, realVector norms) {
		for (std::vector<T>& line : matrix) {
			if (line.size() != norms.size()) {
				makeWarning("utils_la::normLinearTransform(...)","Line has incorrect size: " + std::to_string(line.size()) + " (Should be "+std::to_string(norms.size()) +").");
			}
			for (size_t i = 0; i < line.size(); ++i) {
				line[i] /= norms[i];
			}
		}
		return true;
	}

	realMatrix double_real_part_matrix(const realMatrix& matrix) {
		if (matrix.size() == 0) {
			makeWarning("utils_la::double_real_part_matrix","Matrix of zero size received.");
			return realMatrix();
		}
		const size_t dimX = matrix.size();
		const size_t dimY = matrix[0].size();
		realMatrix retVal(2*dimX, realVector(2*dimY, 0.));
		for (size_t i = 0; i < dimX; ++i) {
			if (matrix[i].size() != dimY) {
				makeWarning("utils_la::double_real_part_matrix","Non-square matrix received.");
				retVal[2*i  ].resize(2*matrix[i].size());
				retVal[2*i+1].resize(2*matrix[i].size());
			}
			for (size_t j = 0; j < matrix[i].size(); ++j) {
				retVal[2*i  ][2*j  ] = matrix[i][j];
				retVal[2*i+1][2*j+1] = matrix[i][j];
			}
		}
		return retVal;
	}
	
	template <typename T>
	std::vector<std::vector<T> > transpose(const std::vector<std::vector<T> >& matrix) {
		const size_t dimX = matrix.size();
		if (dimX == 0) {
			makeWarning("utils_la::transpose","Zero-sized matrix received.");
			return std::vector<std::vector<T> >();
		}
		const size_t dimY = matrix[0].size();
		if (dimY == 0) {
			makeWarning("utils_la::transpose","Zero-sized lines in matrix received.");
			return std::vector<std::vector<T> >();
		}
		std::vector<std::vector<T> > retVal(dimY, std::vector<T>(dimX));
		for (size_t i = 0; i < dimX; ++i) {
			if (matrix[i].size() != dimY) {
				makeError("utils_la::transpose","Matrix lines of varying size found. Cannot transpose.");
				return std::vector<std::vector<T> >();
			}
			for (size_t j = 0; j < dimY; ++j) {
				retVal[j][i] = matrix[i][j];
			}
		}
		return retVal;
	}

}
