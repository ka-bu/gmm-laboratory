#ifndef GAUSSIANMIXTURE_H
#define GAUSSIANMIXTURE_H

#include "gaussian.h"
#include "../../base.h"
#include "../../base/parameters.h"
#include "../../base/commonutil.h"

#include <random>
#include <boost/concept_check.hpp>



/**
 * @brief Gaussian mixture model distribution.
 */
class GaussianMixture
{
public:

	GaussianMixture(Parameters const&);

	/**
	 * evaluates the density of the GMM distribution at the given point x.
     */
	virtual fp_type density(Vector const& x) const;

	/**
	 * computes the negative log-likelihood of the density at the given point x.
     */
	virtual fp_type nll(Vector const& x) const;
	
	virtual fp_type minNLL(Vector const& x) const;
	
	virtual fp_type minSquaredMahalanobis(Vector const& x) const;
	
	template<typename RndEngine> Vector draw(RndEngine&) const;
	
	template<typename RndEngine> commonutil::CompleteData drawCompleteData(RndEngine& re) const;

private:
	Parameters gmmdesc;
	std::vector<SingleGaussian> gaussians;
};

template<typename RndEngine>
Vector GaussianMixture::draw(RndEngine& re) const
{
	std::size_t k = this->gmmdesc.weights.size();
	
	assert(this->gaussians.size()==k);	
	
	if (k==0)
		throw "There are no gaussians in this mixture!";

	std::uniform_real_distribution<> urd(0, 1);
	fp_type r = urd(re);
	
	std::size_t i = 0;
	fp_type w = this->gmmdesc.weights[0];
	while (r>w&&i<k)
	{
		r -= w;
		w = this->gmmdesc.weights[++i];
	}
				
	return this->gaussians[i].draw(re);
}

template<typename RndEngine>
commonutil::CompleteData GaussianMixture::drawCompleteData(RndEngine& re) const
{
	std::size_t k = this->gmmdesc.weights.size();
	
	assert(this->gaussians.size()==k);	
	
	if (k==0)
		throw "There are no gaussians in this mixture!";

	std::uniform_real_distribution<> urd(0, 1);
	fp_type r = urd(re);
	
	std::size_t i = 0;
	fp_type w = this->gmmdesc.weights[0];
	while (r>w&&i<k)
	{
		r -= w;
		w = this->gmmdesc.weights[++i];
	}
				
	commonutil::CompleteData cd;
	cd.point = this->gaussians[i].draw(re);
	cd.source = i;
	return cd;
}

#endif 
