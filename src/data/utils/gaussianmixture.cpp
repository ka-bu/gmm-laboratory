#include "gaussianmixture.h"

#include <cmath>

#include "../../base/linalgutil.h"

GaussianMixture::GaussianMixture(Parameters const& desc)
{
	this->gmmdesc = desc;
	std::size_t k = desc.weights.size();
	
	assert(desc.means.cols()==k && desc.covariances.size()==k);
	
	for (std::size_t i=0; i<k; ++i)
		this->gaussians.push_back(SingleGaussian(desc.means.col(i), desc.covariances[i]));
}

fp_type GaussianMixture::density(Vector const& x) const
{
	std::size_t k = this->gmmdesc.weights.size();
	
	assert(this->gaussians.size()==k);
	
	fp_type sum = 0;
	for (std::size_t i=0; i<k; ++i)
		sum += this->gmmdesc.weights[i]*this->gaussians[i].density(x);
		
	return sum;
}

fp_type GaussianMixture::nll(Vector const& x) const
{
	return -log(this->density(x));
}

fp_type GaussianMixture::minNLL(Vector const& x) const
{
	std::size_t k = this->gmmdesc.weights.size();
	
	assert(this->gaussians.size()==k);	
	
	fp_type ret;
	if (k>0)
	{
		ret = this->gaussians[0].nll(x);
		for (std::size_t i=1; i<k; ++i)
		{
			fp_type c = this->gaussians[i].nll(x);
			if (ret>c)
				ret = c;
		}
	}
	else
		ret = 0;
		
	return ret;
}

fp_type GaussianMixture::minSquaredMahalanobis(Vector const& x) const
{
	std::size_t k = this->gmmdesc.weights.size();
	
	assert(this->gaussians.size()==k);	
	
	fp_type ret;
	if (k>0)
	{
		ret = this->gaussians[0].squaredMahalanobis(x);
		for (std::size_t i=1; i<k; ++i)
		{
			fp_type c = this->gaussians[i].squaredMahalanobis(x);
			if (ret>c)
				ret = c;
		}
	}
	else
		ret = 0;
		
	return ret;
}

