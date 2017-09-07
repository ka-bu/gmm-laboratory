#include "datautil.h"

#include <cmath>

#include "../../base/linalgutil.h"


void datautil::stats(commonutil::DataSet const& ds, DataStats& dataStats )
{
// spreads, maxSpread, totalWeight, minWeight
	dataStats.spreads = (ds.points.rowwise().maxCoeff()-ds.points.rowwise().minCoeff());
	dataStats.maxSpread = dataStats.spreads.maxCoeff();
	dataStats.totalWeight = ds.weights.sum();
	if (dataStats.totalWeight<=0)
		gmmlab_throw("Total input weight less or equal to zero.");
	dataStats.minWeight = dataStats.totalWeight;
	idx_type n = ds.n();
	for (idx_type i=0; i<n; ++i)
		if (ds.weights[i]>0 && ds.weights[i]<dataStats.minWeight)
			dataStats.minWeight = ds.weights[i];
}

Matrix datautil::drawFromSG(Vector const& mean, Matrix const& cov, std::mt19937& gen, idx_type n, bool isCholesky)
{
	Matrix covchol = cov;
	if (!isCholesky)
		linalg::cholesky(covchol);

	idx_type d = mean.size();
	std::normal_distribution<> nd(0,1);

	Matrix samples(d,n);
	for (idx_type j = 0; j < n; ++j)
	{
		for (idx_type i = 0; i < d; ++i)
			samples(i,j) = nd(gen);
		samples.col(j) = covchol*samples.col(j) + mean;
	}
	
	return samples;
}



