#ifndef GAUSSIAN_H
#define GAUSSIAN_H

#include "../../base.h"
#include "../../base/parameters.h"

#include <random>
#include <boost/concept_check.hpp>

/**
 * @brief Multivariate normal distribution.
 */
class SingleGaussian
{
public:

	SingleGaussian(Vector const&, Matrix const&);
	
	/**
	 * evaluates the density of the multivariate normal distribution at the given data point x.
     */
	fp_type density(Vector const& x) const;
	
	/**
	 * computes the negative log-likelihood of the density at the given data point x.
     */
	fp_type nll(Vector const& x) const;

	fp_type squaredMahalanobis(Vector const& x) const;
	
	template<typename RndEngine> Vector draw(RndEngine&) const;
	
private:
	friend std::ostream& operator<<(std::ostream&, SingleGaussian const&);
	Vector mean;
//	Matrix const& covariance;
	Matrix cholesky, inverseCholesky;
	fp_type logSqrt; // log of sqrt of determinant of 2*pi*covariance
};


template<typename RndEngine>
Vector SingleGaussian::draw(RndEngine& re) const
{
	idx_type d = this->mean.size();
	std::normal_distribution<> nd(0,1);

	Vector sample(d);
	for (idx_type i=0; i<d; ++i)
		sample[i] = nd(re);
	Vector prod = this->cholesky*sample;

	for (idx_type i=0; i<d; ++i)
		sample[i] = prod[i]+this->mean[i];
		
	return sample;
}

std::ostream& operator<<(std::ostream&, SingleGaussian const&);


#endif // ifndef GAUSSIAN_H