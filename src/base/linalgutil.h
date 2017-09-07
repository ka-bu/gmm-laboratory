#ifndef LINALGUTIL_H
#define LINALGUTIL_H

#include "../base.h"

namespace linalg
{
	
	
	void inv_sqrt(Matrix& m);
	
	/**
	 * Computes the Cholesky ("Scholeskie") decomposition of a symmetric, positive definit matrix M.
	 * That is, it computes a lower triangle matrix G such that M = G*G^T. The matrix G is stored in M.
	 * @param m symmetric, positive definit matrix
	 * @return lower triangle G (where M = G*G^T) is stored in M
	 */
	void cholesky(Matrix& m);
	
	/**
	 * Tests whether the given matrix is symmetric and positive definit.
	 * @param m a matrix
	 * @return true iff the given matrix M is symmetric and positive definit.
	 */
	bool spd(Matrix m, bool verbose = false);
	
	bool allspd(std::vector<Matrix> const& matrices);
		
	/**
	 * Computes the Gram–Schmidt process.
	 * @param m a matrix containing the vectors to be orthonormalized as column vectors
	 * @return orthonormalized vectors are stored in M as column vectors
	 */
	void gramschmidt(Matrix& m);
	
	/**
	 * Computes the inverse of the given lower triangle matrix m by backwards substitution and stores the result
	 * in the lower triangle of m.
	 * @param m lower triangle matrix
	 * @return inverse of m is stored in m
	 */
	void ltrinv(Matrix& m);
	
	/**
	 * Computes the minimum and maximum eigenvalue of m.
	 * @param m square matrix
	 * @param minEw minimum eigenvalue
	 * @param maxEW maximum eigenvalue
	 */
	void minmaxEigenvalue(Matrix m, fp_type& minEW, fp_type& maxEW);
	
}

#endif
