#ifndef KMEANSUTIL_H
#define KMEANSUTIL_H

#include "../../../base.h"
#include <iostream>

#include "../../gmm/utils/gmmutil.h"
#include "../../../base/parameters.h"

namespace kmeansutil
{


	/**
	 * wrap some means into a GMMDesc. That is, the means are used as centers of the components,
	 * while the weights and covariances are estimated by considering the kmeans clusters given by the means.
	 */
	Parameters wrapMeans(commonutil::DataSet const& input, Matrix const& means, bool computeWeightAndCovar, bool verbose = false);
  
	/**
	 * computes the sum of the squared distances between a point and its closest mean (weighted by the weight of the point).
	 */
	fp_type kmeanscost(commonutil::DataSet const& data, Matrix const& means);
	

	Vector squaredDistances(Matrix const& points, Vector const& point,  std::vector<idx_type> const& indices_of_input_points_to_be_used = std::vector<idx_type>());
	
	fp_type minSquaredDistance(Matrix const& points, Vector const& point);
	
	Vector minSquaredDistances(Matrix const& points, Matrix const& means, std::vector<idx_type> const& indices_of_input_points_to_be_used = std::vector<idx_type>());
	
	idx_type indexByAdaptiveDensities(commonutil::DataSet const& data, Matrix const& means, fp_type ratio, std::mt19937& gen, std::vector<idx_type> const& indices_of_input_points_to_be_used = std::vector<idx_type>());
	
	
	/**
	 * returns the assignment of each point to its nearest mean.
	 */
	std::vector<idx_type> kmeansPartition(Matrix const& points, Matrix const& means, std::vector<idx_type> const& indices_of_input_points_to_be_used = std::vector<idx_type>());
	
	/**
	 * returns the index of the mean that has the smalles squared distance to the point.
	 */
	idx_type nearest(Vector const& point, Matrix const& means);
	
	/**
	 * the density of a point is its squared distance to its nearest center (normalized by 
	 * the sum of all squared distances * 2) plus (1-ratio)/2 times uniform distribution.
	 */
	Vector adaptiveDensities(commonutil::DataSet const& data, Matrix const& means, fp_type const ratio);

	/**
	 * the density of a point is its squared distance  to its nearest center (normalized by 
	 * the sum of all squared distances * 2) plus (1-ratio)/2 times uniform distribution.
	 */
	Vector adaptiveDensities(commonutil::DataSet const& data, Matrix const& means, fp_type const ratio, std::vector<idx_type> const& indices_of_input_points_to_be_used);
	
}


#endif
