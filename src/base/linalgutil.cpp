#include "linalgutil.h"

#include "../settings/settings.h"

#include <iostream>


void linalg::inv_sqrt(Matrix& m)
{
	linalg::cholesky(m);
	linalg::ltrinv(m);
}



	
void linalg::cholesky(Matrix& m)
{
	idx_type n = m.rows();
	if (m.cols()!=n)
		gmmlab_throw("no square matrix");

	for (idx_type i=0; i<n; ++i)
	{
		fp_type sum;
		for (idx_type j=0; j<i; ++j)
		{
			sum = m(i,j);
			for (idx_type k=0; k<j; ++k)
				sum -= m(i,k)*m(j,k);
			m(i,j) = sum/m(j,j);
		}

		sum = m(i,i);
		for (idx_type k=0; k<i; ++k)
		{
			sum -= m(i,k)*m(i,k);
			m(k,i) = 0.0;	// delete the entries above the diagonal
		}
		if (sum > 0.0)
			m(i,i) = sqrt(sum);
		else
		{
			if (commonSettings().verbose)
				std::cout << "linalg::cholesky(Matrix& m) - m(" << i << "," << i << ") := sqrt(" << sum << ") !!!" << std::endl;
			gmmlab_throw("Matrix is not positive definite.");
		}
	}
}


bool linalg::spd(Matrix m, bool verbose)
{
	idx_type n = m.rows();
	if (m.cols()!=n)
		return false;

	// compute the Cholesky decomposition
	for (idx_type i=0; i<n; ++i)
	{
		fp_type sum;
		for (idx_type j=0; j<i; ++j)
		{
			sum = m(i,j);
			for (idx_type k=0; k<j; ++k)
				sum -= m(i,k)*m(j,k);
			m(i,j) = sum/m(j,j);
		}

		sum = m(i,i);
		for (idx_type k=0; k<i; ++k)
		{
			sum -= m(i,k)*m(i,k);
			m(k,i) = 0.0;	// delete the entries above the diagonal
		}

		if (sum > 0.0)
			m(i,i) = sqrt(sum);
		else
		{
			if (verbose)
				std::cout << "linalg::spd(Matrix m) - m(" << i << "," << i << ") := sqrt(" << sum << ") !!!" << std::endl;
			return false;
		}
	}

	return true;
}

bool linalg::allspd(std::vector<Matrix> const& matrices)
{
	bool ret = true;
	std::size_t size = matrices.size();
	for (std::size_t i=0; i<size; ++i)
		ret = ret && linalg::spd(matrices[i]);
	return ret;
}

void linalg::gramschmidt(Matrix& m)
{
	idx_type r = m.rows();
	idx_type c = m.cols();

	for (idx_type i=0; i<c; ++i)
	{
		for (idx_type j=0; j<i; ++j)
		{
			fp_type s = 0;
			for (size_t k=0; k<r; ++k)
				s += m(k,j)*m(k,i);
			for (size_t k=0; k<r; ++k)
				m(k,i) -= s*m(k,j);
		}

		fp_type normalizer = 0;
		for (idx_type k=0; k<r; ++k)
			normalizer += m(k,i)*m(k,i);
		normalizer = sqrt(normalizer);

		for (idx_type k=0; k<r; ++k)
			m(k,i) /= normalizer;
	}
}


void linalg::ltrinv(Matrix& m)
{
	idx_type n = m.rows();
	// compute inverse of lower triangle matrix m by backwards substitution and
	// store the result in the upper triangle of m
	for (idx_type i=n; i-- > 0; )
	{
//		inv(i,i) = 1/lt(i,i);
		m(i,i) = 1/m(i,i);
		for (idx_type j=i+1; j<n; ++j)
		{
			fp_type sum = 0;
			for (idx_type k=i; k<j; ++k)
			{
//				sum += lt(j,k)*inv(k,i);
				sum += m(j,k)*m(i,k);
			}
//			inv(j,i) = -sum*inv(j,j);
			m(i,j) = -sum*m(j,j);
		}
	}

	// transpose the resulting matrix and wipe out the upper triangle
	for (idx_type i=0; i<n; ++i)
		for (idx_type j=0; j<i; ++j)
		{
			m(i,j) = m(j,i);
			m(j,i) = 0.0;
		}
}


void linalg::minmaxEigenvalue(Matrix m, fp_type& minEW, fp_type& maxEW)
{
   if (m.cols()!=m.rows())
		gmmlab_throw("no square matrix");
   Eigen::SelfAdjointEigenSolver<Matrix> eigensolver;
   // compute only the eigenvalues, not eigenvectors
   eigensolver.compute(m,Eigen::EigenvaluesOnly);
   if (eigensolver.info() != Eigen::Success)
        gmmlab_throw("eigensolver failed");
   minEW = eigensolver.eigenvalues().minCoeff();
   maxEW = eigensolver.eigenvalues().maxCoeff();
}

