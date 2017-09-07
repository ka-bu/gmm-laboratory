#include "parameters.h"

bool Parameters::empty() const
{
	return (this->weights.size()==0 && this->covariances.size()==0 && this->means.cols()==0);
}

bool Parameters::operator==(Parameters const& rhs) const
{
	return (
			// check size
			this->weights.size()==rhs.weights.size()
			&& this->means.rows()==rhs.means.rows()
			&& this->means.cols()==rhs.means.cols()
			&& this->covariances.size()==rhs.covariances.size()
			// we assume that the covariance matrices of this and rhs have the same size if the means have the same size
			// compare content componentwise
			&& this->weights==rhs.weights
			&& this->means==rhs.means
			&& this->covariances==rhs.covariances);

}

std::ostream& operator<<(std::ostream& os, Parameters const& desc)
{
	std::size_t k = desc.means.cols();

	assert((desc.weights.size()==k && desc.covariances.size()==k) || (desc.weights.size()==0 && desc.covariances.size()==0));

	if(desc.weights.size() == k)
	{
		os << "gaussian mixture with " << k << " components" << std::endl;
		for(std::size_t i=0; i<k; ++i)
		{
			os << i+1 << ". weight = " << desc.weights[i] << std::endl
				<< "mean:" << std::endl << desc.means.col(i) << std::endl;
			if(desc.covariances.size()>0)
				os << "covariance:" << std::endl << desc.covariances[i] << std::endl;
			if(desc.kappa.size()>0)
			        os << "kappa:" << std::endl << desc.kappa[i] << std::endl << std::endl;
		}
	}
	else
	{
		os <<  k << " means" << std::endl;
		for(std::size_t i=0; i<k; ++i)
			os << i+1 << ". mean:" << std::endl << desc.means.col(i) << std::endl << std::endl;
		
	}
	
	return os;
}

idx_type Parameters::components() const
{
	return this->means.cols(); // also works for "kmeans" algorithms!
}