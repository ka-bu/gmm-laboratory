#include "kmeansutil.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <ctime>

#include "../../gmm/utils/gmmutil.h"
#include "../../../base/linalgutil.h"
#include "../../../settings/settings.h"



Parameters kmeansutil::wrapMeans(commonutil::DataSet const& input, Matrix const& means, bool computeWeightAndCovar, bool verbose)
{
	Parameters desc;

	std::size_t k = means.cols();

	if (k==0)
		gmmlab_throw("No means given.");

	idx_type n = input.points.cols();
	idx_type d = input.points.rows();

	if (n==0 || d==0)
		gmmlab_throw("Input is empty.");

	fp_type totalWeight = input.weights.sum();
	if (totalWeight<=0)
		gmmlab_throw("Total input weight less or equal to zero.");

	desc.means = means;
	
	if(computeWeightAndCovar){
	
	  // assign each point to its nearest center
	  Matrix pmatrix(k,n);
	  for (idx_type i=0; i<k; ++i)
		  pmatrix.row(i) = (input.points.colwise()-means.col(i)).colwise().squaredNorm();
	  std::vector<idx_type> indices(n);
	  for (idx_type i=0; i<n; ++i)
		  pmatrix.col(i).minCoeff(&indices[i]);
	  
	  // estimate weights and covariances
	  desc.weights = Vector::Zero(k);
	  desc.covariances = std::vector<Matrix>(k, Matrix::Zero(d,d));
	  std::vector<fp_type> sphericals(k, 0);
	  for (std::size_t i=0; i<n; ++i)
	  {
		  desc.weights[indices[i]] += input.weights[i];
		  Vector y = input.points.col(i) - desc.means.col(indices[i]);
		  desc.covariances[indices[i]].noalias() += input.weights[i]*(y*y.transpose());
		  sphericals[indices[i]] += input.weights[i]*y.squaredNorm();
	  }
	  for (std::size_t i=0; i<k; ++i)
	  {
		  if (desc.weights[i] > 0)
		  {
			  desc.covariances[i] = desc.covariances[i]/desc.weights[i];
			  if (!linalg::spd(desc.covariances[i]))
			  {
				  desc.covariances[i] = Matrix::Identity(d,d)*(sphericals[i]/(desc.weights[i]*d));
				  if (sphericals[i]/(desc.weights[i]*d) <= 0)
				  {
					  desc.covariances[i] = Matrix::Identity(d,d);
					  if (verbose)
						  std::cout << "replaced non-spd covariance " << i+1 << " with unit sphere covariance." << std::endl;
				  }
				  else
					  if (verbose)
						  std::cout << "replaced non-spd covariance " << i+1 << " with spherical covariance." << std::endl;
			  }
		  }
		  else
		  {
			  desc.weights[i] = totalWeight/k;
			  desc.covariances[i] = Matrix::Identity(d,d);
			  if (verbose)
				  std::cout << "kmeansutil::kmeansToGMM() - replaced empty covariance " << i+1 << " with unit sphere covariance." << std::endl;
		  }

	  }

	  // finally normalize the weights
	  desc.weights /= desc.weights.sum();
	
	}
	
	return desc;
}

fp_type kmeansutil::kmeanscost(commonutil::DataSet const& data, Matrix const& means)
{
	Vector v = kmeansutil::minSquaredDistances(data.points, means);
	idx_type n = data.points.cols();

	assert(v.size() == n && data.weights.size() == n);

	for (idx_type i = 0; i < n; ++i)
		v[i] *= data.weights[i];
	
	return v.sum();
}



			
				
// const fp_type kmeansutil::minSquaredDistance(Matrix const& points)
// {
// 	fp_type min = 0;
// 	idx_type n = points.cols();
// 
// 	for (idx_type i = 1; i < n; ++i)
// 	{
// 		Vector sn = (points.rightCols(n - i).colwise() - points.col(i - 1)).colwise().squaredNorm();
// 		idx_type m = sn.size();
// 
// 		for (idx_type j = 0; j < m; ++j)
// 			if (sn[j] > 0 && (min == 0 || sn[j] < min))
// 				min = sn[j];
// 	}
// 
// 	return min;
// }

Vector kmeansutil::squaredDistances(Matrix const& points, Vector const& point, std::vector<idx_type> const& indices_of_input_points_to_be_used)
{
	if(indices_of_input_points_to_be_used.empty())
		return (points.colwise() - point).colwise().squaredNorm();
	
	// else

	idx_type n = indices_of_input_points_to_be_used.size();
	Vector distances(n);
	for(idx_type i=0; i<n; ++i)
		distances(i) = (points.col(indices_of_input_points_to_be_used.at(i)) - point).squaredNorm();
	return distances;
}

				
fp_type kmeansutil::minSquaredDistance(Matrix const& points, Vector const& point)
{
	Vector dists = kmeansutil::squaredDistances(points, point);
	idx_type index;
	dists.minCoeff(&index);
	return dists[index];
}


idx_type kmeansutil::indexByAdaptiveDensities(commonutil::DataSet const& data, Matrix const& means, fp_type ratio, std::mt19937& gen, std::vector<idx_type> const& indices_of_input_points_to_be_used)
{
	const bool indices_given = !indices_of_input_points_to_be_used.empty();
	
	idx_type index;
	if(means.cols() == 0)
	{
		unsigned int n = indices_given ? indices_of_input_points_to_be_used.size() : data.points.cols();
		std::uniform_int_distribution<idx_type> uid(0,n-1);
		index = uid(gen);
		if(indices_given)
			index = indices_of_input_points_to_be_used.at(index);
		return index;
	}
	
	Vector densities = adaptiveDensities(data, means, ratio, indices_of_input_points_to_be_used);
	index = commonutil::randomIndex<std::mt19937>(densities,gen);
	if(indices_given)
		index = indices_of_input_points_to_be_used.at(index);
	return index;
}


Vector kmeansutil::minSquaredDistances(Matrix const& points, Matrix const& means, std::vector<idx_type> const& indices_of_input_points_to_be_used)
{
	const bool indices_given = !indices_of_input_points_to_be_used.empty();
	const idx_type n = indices_given ? indices_of_input_points_to_be_used.size() : points.cols();
	Vector dists = Vector::Zero(n);
	for (idx_type i = 0; i < n; ++i)
		dists[i] = minSquaredDistance(means, points.col( indices_given ? indices_of_input_points_to_be_used.at(i) : i));
	return dists;
}


std::vector<idx_type> kmeansutil::kmeansPartition(Matrix const& points, Matrix const& means, std::vector<idx_type> const& indices_of_input_points_to_be_used)
{
	if(indices_of_input_points_to_be_used.empty())
	{
		idx_type n = points.cols();
		std::vector<idx_type> partition(n);
		for(idx_type i=0; i<n; ++i)
			partition.at(i) = (kmeansutil::nearest(points.col(i), means));
		return partition;
	}
	else
	{
		idx_type n = indices_of_input_points_to_be_used.size();
		std::vector<idx_type> partition(n);
		for(idx_type i=0; i<n; ++i)
			partition.at(i) = (kmeansutil::nearest(points.col(indices_of_input_points_to_be_used.at(i)), means));
		return partition;
	}
}


idx_type kmeansutil::nearest(Vector const& point, Matrix const& means)
{
	idx_type index;
	kmeansutil::squaredDistances(means, point).minCoeff(&index);
	return index;
}


Vector kmeansutil::adaptiveDensities(commonutil::DataSet const& data, Matrix const& means, fp_type const ratio)
{
	assert (ratio>=0 && ratio<=1);

	Vector densities = kmeansutil::minSquaredDistances(data.points, means);
	idx_type n = data.points.cols();

	for (idx_type i=0; i<n; ++i)
		if (!(densities[i]>0)) // test if <0 or NaN
			densities[i] = 0;

	fp_type sum = densities.sum();
	if (sum>0)
		densities /= sum;

	fp_type iidDens = fp_type(1)/n;
	for (idx_type j=0; j<n; ++j)
	{
		densities[j] = ratio*densities[j]+(1-ratio)*iidDens;
		densities[j] *= data.weights[j];
	}

	sum = densities.sum();
	if (sum>0)
		densities /= sum;

	return densities;
}

Vector kmeansutil::adaptiveDensities(commonutil::DataSet const& data, Matrix const& means, fp_type const ratio, std::vector<idx_type> const& indices_of_input_points_to_be_used)
{
	assert (ratio>=0 && ratio<=1);

	Vector densities = kmeansutil::minSquaredDistances(data.points, means, indices_of_input_points_to_be_used);
	idx_type n =indices_of_input_points_to_be_used.size();

	for (idx_type i=0; i<n; ++i)
		if (!(densities[i]>0)) // test if <0 or NaN
			densities[i] = 0;

	fp_type sum = densities.sum();
	if (sum>0)
		densities /= sum;

	fp_type iidDens = fp_type(1)/n;
	for (idx_type j=0; j<n; ++j)
	{
		densities[j] = ratio*densities[j]+(1-ratio)*iidDens;
		densities[j] *= data.weights[j];
	}

	sum = densities.sum();
	if (sum>0)
		densities /= sum;

	return densities;
}


