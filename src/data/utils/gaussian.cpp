#include "gaussian.h"

#include <cmath>

#include "../../base/linalgutil.h"

SingleGaussian::SingleGaussian(Vector const& m, Matrix const& cov) : mean(m)
//, covariance(cov)
{
	this->cholesky = cov;
	linalg::cholesky(this->cholesky);
	this->inverseCholesky = this->cholesky;
	linalg::ltrinv(this->inverseCholesky);

	idx_type d = m.size();
	this->logSqrt = d*log(2*M_PI)/2;
	for (idx_type i=0; i<d; ++i)
		this->logSqrt += log(this->cholesky(i,i));
}

fp_type SingleGaussian::density(Vector const& x) const
{
	return exp(-0.5*this->squaredMahalanobis(x) - logSqrt);	
}

fp_type SingleGaussian::nll(Vector const& x) const
{
	return 0.5*this->squaredMahalanobis(x) + logSqrt;	
}

fp_type SingleGaussian::squaredMahalanobis(Vector const& x) const
{
	// transform the given point
	idx_type d = x.size();
	Vector y = this->inverseCholesky*(x-this->mean);
	fp_type qf = 0.0;
	for (idx_type i=0; i<d; ++i)
		qf += y[i]*y[i];

	return qf;	
}

std::ostream& operator<<(std::ostream& os, SingleGaussian const& sg)
{
	os	<< "single gaussian with" << std::endl
		<< ".\tmean=" << sg.mean << std::endl
		<< "\tcholesky=" << sg.cholesky << std::endl
		<< "\tinverseCholesky=" << sg.inverseCholesky << std::endl
		<< "\tlogSqrt=" << sg.logSqrt << std::endl;
	return os;
}

