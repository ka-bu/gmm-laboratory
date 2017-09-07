#include "similarities.h"

#include "../algorithm/gmm/utils/gmmutil.h"

fp_type MahalanobisDis::operator()(Vector const& m1, Matrix const& c1, Vector const& m2, Matrix const& c2)
{
	Matrix mixcov = c1+c2;
	Vector x = m1-m2;
	return x.transpose()*mixcov*x;
}


Vector GMMDis::operator()(Parameters const& gmm1, Parameters const& gmm2){
	
	std::size_t k = gmm1.weights.size();

	if(k != gmm2.weights.size())
		gmmlab_throw("gmms have different number of components"); 
	
	// match components by shortest distance between means
	Vector match(k);
	for(size_t i=0; i<k; ++i)
	{
		Vector distances(k);
		for(size_t j=0; j<k; ++j)
		{
			if(j==i)
				continue;
			Vector tmp = gmm1.means.col(i)-gmm2.means.col(j);
			distances(j) = abs(tmp.norm());
		}
		size_t m;
		distances.minCoeff(&m);
		match(i) = m;
	}

	// output = tuple of weight-, mean-, covariance-distance
	Vector output = Vector::Zero(3);
	for(size_t l = 1; l<k; ++l)
	{
		output(0) += abs(gmm2.weights(l)-gmm1.weights(match(l)));
		Vector tmpVector = gmm2.means.col(l)-gmm1.means.col(match(l));
		output(1) += tmpVector.squaredNorm();
		Matrix tmpMatrix = gmm2.covariances.at(l)-gmm1.covariances.at(match(l));
		output(2) += tmpMatrix.squaredNorm();
	}

	return output;
}

fp_type AverageLinkage::operator()(const Vector& sum1, idx_type count1, const Vector& sum2, idx_type count2)
{
	if(count1 + count2 <= 0)
		gmmlab_throw("count1+count2 <= 0")
	return 1./(count1 + count2) * sum1.dot(sum2);
}

