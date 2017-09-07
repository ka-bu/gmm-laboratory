#include "gmmutil.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <ctime>
#include <cmath>

#include "../../kmeans/utils/kmeansutil.h"
// #include "../../../data/utils/gaussian.h"
#include "../../../base/linalgutil.h"
#include "../../../base/commonutil.h"
#include "../../../settings/settings.h"


Parameters gmmutil::meansToGMM(commonutil::DataSet const& input, Matrix const& means, bool spherical, bool verbose, std::vector<idx_type> const& indices_of_input_points_to_be_used)
{
	
	Parameters desc;

	idx_type k = means.cols();

	if (k==0)
		gmmlab_throw("No means given.");

	bool indices_given = !indices_of_input_points_to_be_used.empty();
	idx_type n = indices_given ? indices_of_input_points_to_be_used.size() : input.points.cols();
	idx_type d = input.points.rows();

	if (n==0 || d==0)
		gmmlab_throw("Input is empty.");

	fp_type totalWeight = input.weights.sum();
	if (totalWeight<=0)
		gmmlab_throw("Total input weight less or equal to zero.");

	
// 	Matrix pmatrix(k,n);
// 	for (idx_type i=0; i<k; ++i)
// 		pmatrix.row(i) = (input.points.colwise()-means.col(i)).colwise().squaredNorm();
// 
// 	std::vector<idx_type> indices(n);
// 	for (idx_type i=0; i<n; ++i)
// 		pmatrix.col(i).minCoeff(&indices[i]);
	
	std::vector<idx_type> indices = kmeansutil::kmeansPartition(input.points, means, indices_of_input_points_to_be_used);
	
	desc.weights = Vector::Zero(k);
	desc.means = Matrix::Zero(d,k);
	desc.covariances = std::vector<Matrix>(k, Matrix::Zero(d,d));
	std::vector<fp_type> sphericals(k, 0);

	for (std::size_t i=0; i<n; ++i)
	{
		idx_type index = indices_given ? indices_of_input_points_to_be_used.at(i) : i;
		desc.weights[indices[i]] += input.weights[index];
		desc.means.col(indices[i]).noalias() += input.weights[index]*input.points.col(index);
	}

	for (std::size_t i=0; i<k; ++i)
		if (desc.weights[i] > 0)
			desc.means.col(i) = desc.means.col(i)/desc.weights[i];
		else
			desc.means.col(i) = means.col(i);

	if(spherical)
	{
		for (std::size_t i=0; i<n; ++i)
		{
			idx_type index = indices_given ? indices_of_input_points_to_be_used.at(i) : i;
			Vector y = input.points.col(index) - desc.means.col(indices[i]);
			sphericals[indices[i]] += input.weights[index]*y.squaredNorm();
		}
		for (std::size_t i=0; i<k; ++i)
		{
			if (desc.weights[i] > 0)
			{
				desc.covariances[i] = Matrix::Identity(d,d)*(sphericals[i]/(desc.weights[i]*d));
				if (!linalg::spd(desc.covariances[i], verbose))
				{
					desc.covariances[i] = Matrix::Identity(d,d);
					if (verbose)
						std::cout << "replaced non-spd covariance " << i+1 << " with unit sphere covariance." << std::endl;
				}
			}
			else
			{
				desc.weights[i] = totalWeight/k;
				desc.covariances[i] = Matrix::Identity(d,d);
				if (verbose)
					std::cout << "replaced empty covariance " << i+1 << " with unit sphere covariance." << std::endl;
			}

		}
	}
	else
	{
		for (std::size_t i=0; i<n; ++i)
		{
			idx_type index = indices_given ? indices_of_input_points_to_be_used.at(i) : i;
			Vector y = input.points.col(index) - desc.means.col(indices[i]);
			desc.covariances[indices[i]].noalias() += (input.weights[i]*y)*y.transpose();
			sphericals[indices[i]] += input.weights[index]*y.squaredNorm();
		}

		for (std::size_t i=0; i<k; ++i)
		{
			if (desc.weights[i] > 0)
			{
				desc.covariances[i] = desc.covariances[i]/desc.weights[i];
				if (!linalg::spd(desc.covariances[i], verbose))
				{
					desc.covariances[i] = Matrix::Identity(d,d)*(sphericals[i]/(desc.weights[i]*d));
					if (!linalg::spd(desc.covariances[i], verbose))
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
					std::cout << "replaced empty covariance " << i+1 << " with unit sphere covariance." << std::endl;
			}

		}
	}
	
	// finally normalize the weights
	desc.weights /= desc.weights.sum();
	
	return desc;
}


fp_type gmmutil::nll(commonutil::DataSet const& data, Parameters const& desc, std::vector<idx_type> const& indices_of_input_points_to_be_used)
{	
	const bool indices_given = !indices_of_input_points_to_be_used.empty();
	
	Vector v = gmmNLL(data.points, desc, indices_of_input_points_to_be_used);
	idx_type n = v.size();
	
	for (idx_type i=0; i<n; ++i)
		v[i] *= data.weights[ indices_given ? indices_of_input_points_to_be_used.at(i) : i];
	
	return v.sum();
}

fp_type gmmutil::nlcdl(commonutil::DataSet const& data, Parameters const& desc)
{		
	idx_type n = data.points.cols();
	idx_type k = desc.components();
	
	Matrix nlls(k,n);
	for(idx_type i=0; i<k; ++i)
		nlls.row(i) = (gmmutil::gaussianDensity(data.points, desc.means.col(i), desc.covariances.at(i) ).array() - log(desc.weights(i))).matrix();
		
	Vector min = nlls.colwise().minCoeff();
	return min.sum();
}

const fp_type gmmutil::minSquaredDistance(Matrix const& points)
{
	fp_type min = 0;
	idx_type n = points.cols();
	for (idx_type i=1; i<n; ++i)
	{
		Vector sn = (points.rightCols(n-i).colwise()-points.col(i-1)).colwise().squaredNorm();
		idx_type m = sn.size();
		for (idx_type j=0; j<m; ++j)
			if (sn[j]>0&&(min==0||sn[j]<min))
				min = sn[j];
	}
	return min;
}

Vector gmmutil::gaussianDensity(Matrix const& points, Vector const& mean, Matrix covariance, std::vector<idx_type> const& indices_of_input_points_to_be_used)
{
	linalg::cholesky(covariance);
	Vector v = covariance.diagonal();
	idx_type d = v.size();
	fp_type logSqrt = d*log(2*M_PI)/2;
	for (idx_type i=0; i<d; ++i)
		logSqrt += log(v(i));
	linalg::ltrinv(covariance);
	
	
	if(indices_of_input_points_to_be_used.empty())
		v.noalias() = (covariance*(points.colwise()-mean)).colwise().squaredNorm();
	else
	{
		idx_type n = indices_of_input_points_to_be_used.size();
		v.resize(n);
		for(idx_type i=0; i<n; ++i)
			v(i) = (covariance*(points.col(indices_of_input_points_to_be_used.at(i))-mean)).squaredNorm();
	}
	
	
	idx_type n = v.size();
	for (idx_type i=0; i<n; ++i)
		v(i) = exp(-0.5*v(i)-logSqrt);
	
	return v;
}

Vector gmmutil::gmmDensity(Matrix const& points, Parameters const& desc, std::vector<idx_type> const& indices_of_input_points_to_be_used)
{
	
	idx_type n = indices_of_input_points_to_be_used.empty() ? points.cols() : indices_of_input_points_to_be_used.size();
	idx_type k = desc.weights.size();
	
	assert(desc.means.cols()==k && desc.covariances.size()==k);
	
	Vector sum = Vector::Zero(n);
	for (idx_type i=0; i<k; ++i)
		sum.noalias() += desc.weights[i]*gaussianDensity(points, desc.means.col(i), desc.covariances.at(i), indices_of_input_points_to_be_used);
	
	return sum;
}


Vector gmmutil::gaussianNLL(Matrix const& points, Vector const& mean, Matrix covariance, const std::vector<idx_type>& indices_of_input_points_to_be_used)
{
	linalg::cholesky(covariance);
	Vector v = covariance.diagonal();
	idx_type d = v.size();
	fp_type logSqrt = d*log(2*M_PI)/2;
	for (idx_type i=0; i<d; ++i)
		logSqrt += log(v(i));
	linalg::ltrinv(covariance);
	
	if(indices_of_input_points_to_be_used.empty())
		v.noalias() = (covariance*(points.colwise()-mean)).colwise().squaredNorm();
	else
	{
		unsigned int n = indices_of_input_points_to_be_used.size();
		v.resize(n);
		for(idx_type i=0; i<n; ++i)
			v(i) = (covariance*(points.col(indices_of_input_points_to_be_used.at(i))-mean)).squaredNorm();
	}
	
	d = v.size();
	for (idx_type i=0; i<d; ++i)
		v(i) = 0.5*v(i)+logSqrt;
	
	return v;
	
// 	Vector nlls = gaussianDensity(points, mean, covariance);
// 	idx_type n = nlls.size();
// 	for(idx_type i=0; i<n; ++i)
// 		nlls(i) = -log(nlls(i));
// 	return nlls;
// 	
}


Vector gmmutil::gmmNLL(Matrix const& points, Parameters const& desc, const std::vector<idx_type> & indices_of_input_points_to_be_used)
{
	bool indices_given = !indices_of_input_points_to_be_used.empty();
	
	idx_type k = desc.weights.size();
	idx_type n = indices_given ? indices_of_input_points_to_be_used.size() : points.cols();
	
	if (k==0)
		return Vector::Zero(n);
	
	Vector densities = gmmDensity(points, desc, indices_of_input_points_to_be_used);
	for (size_t i=0; i<n; ++i)
	{
		fp_type density = densities[i];
		if (density>0 && density<std::numeric_limits<fp_type>::infinity())
			densities[i] = -log(density);
		else
		{
			Vector nlls(k);
			Vector point;
			if(indices_given)
				point = points.col(indices_of_input_points_to_be_used.at(i));
			else 
				point = points.col(i);
			for (idx_type j=0; j<k; ++j)
				nlls[j] = gaussianNLL(point, desc.means.col(j), desc.covariances.at(j))[0] - log(desc.weights[j]);
			idx_type index;
			densities[i] = nlls.minCoeff(&index);
			if(commonSettings().verbose)
				if (density<=0)
					std::cout << "gmmutil::gmmNLL() - Density of point " << i << " was zero!!! Approximated by cost in component " << index+1 << std::endl;
				else
					std::cout << "gmmutil::gmmNLL() - Density of point " << i << " was infinity!!! Approximated by cost in component " << index+1 << std::endl;
		}
	}
	return densities;
}


Vector gmmutil::gaussianMahalanobisOfInvSqrtCovar(Matrix const& points, Vector const& mean, Matrix const& sqrt_of_inverted_covariance)
{
	return (sqrt_of_inverted_covariance*(points.colwise()-mean)).colwise().squaredNorm();
}


Vector gmmutil::gaussianMahalanobis(Matrix const& points, Vector const& mean, Matrix covariance, const std::vector<idx_type> & indices_of_input_points_to_be_used)
{
	
	linalg::cholesky(covariance);
	linalg::ltrinv(covariance);
	
	if(!indices_of_input_points_to_be_used.empty())
	{
		idx_type n = indices_of_input_points_to_be_used.size();	
		Vector output = Vector(n);
		for(idx_type i=0; i<n; ++i)
			output(i) = (covariance*(points.col(indices_of_input_points_to_be_used.at(i))-mean)).squaredNorm();
		return output;
	}
	
	// else
	
	return (covariance*(points.colwise()-mean)).colwise().squaredNorm();
	
}

Vector gmmutil::gmmMixMahalanobis(Matrix const& points, Parameters const& desc)
{
	idx_type n = points.cols();
	idx_type k = desc.weights.size();
	
	Vector sum = Vector::Zero(n);
	for (std::size_t i=0; i<k; ++i)
		sum.noalias() += desc.weights[i]*(-gmmutil::gaussianMahalanobis(points, desc.means.col(i),desc.covariances.at(i))).array().exp().matrix();
	
	for (size_t i=0; i<n; ++i)
		if (sum[i]>0 && sum[i]<std::numeric_limits<fp_type>::infinity())
			sum[i] = -log(sum[i]);
		else
		{
			Vector nlls(k);
			for (idx_type j=0; j<k; ++j)
				nlls[j] = gaussianMahalanobis(points.col(i), desc.means.col(j), desc.covariances.at(j))[0];
			if (sum[i]<=0)
				sum[i] = nlls.minCoeff();
			else
				sum[i] = nlls.maxCoeff();
		}
		
	return sum;
}

Vector gmmutil::gmmMinMahalanobis(Matrix const& points, Parameters const& desc, const std::vector<idx_type> & indices_of_input_points_to_be_used)
{
	idx_type n = indices_of_input_points_to_be_used.empty() ? points.cols() : indices_of_input_points_to_be_used.size();	
	idx_type k = desc.weights.size();
	
	assert(k>0 && n>0);
	
	Matrix mahalanobis(k,n);
	
	for (idx_type i=0; i<k; ++i)
		mahalanobis.row(i) = gmmutil::gaussianMahalanobis(points, desc.means.col(i),desc.covariances.at(i), indices_of_input_points_to_be_used);
	
	return mahalanobis.colwise().minCoeff();
}



idx_type gmmutil::indexByAdaptiveDensities(commonutil::DataSet const& data, Parameters const& desc, fp_type const ratio, std::mt19937& gen, std::vector<idx_type> const& indices_of_input_points_to_be_used)
{
	bool indices_given = !indices_of_input_points_to_be_used.empty();
		
	idx_type index;
	
	if(desc.weights.size() == 0)
	{
		idx_type n = indices_given ? indices_of_input_points_to_be_used.size() : data.points.cols();
		std::uniform_int_distribution<idx_type> uid(0,n-1);
		index = uid(gen);
		if(indices_given)
			index = indices_of_input_points_to_be_used.at(index);
	}
	else
	{
		Vector densities = adaptiveDensities(data, desc, ratio, indices_of_input_points_to_be_used);
		index = commonutil::randomIndex<std::mt19937>(densities,gen);
		if(indices_given)
			index = indices_of_input_points_to_be_used.at(index);
	}
	
	return index;
}

Vector gmmutil::adaptiveDensities(commonutil::DataSet const& data, Parameters const& desc, fp_type const ratio, std::vector<idx_type> const& indices_of_input_points_to_be_used)
{
	bool indices_given = !indices_of_input_points_to_be_used.empty();
	
	assert (ratio>=0 && ratio<=1);
	
	Vector densities = gmmutil::gmmMinMahalanobis(data.points, desc, indices_of_input_points_to_be_used);
	idx_type n = densities.size();
	
	for (idx_type i=0; i<n; ++i)
		if (!(densities[i]>0)) // test if <0 or NaN
			densities[i] = 0;
		
	fp_type sum = densities.sum();
	if (sum>0)
		densities /= sum;
	
	fp_type iidDens = fp_type(1)/n;
	for (idx_type j=0; j<n; ++j)
	{
		densities[j] = ratio*densities[j]+(1.-ratio)*iidDens;
		if(indices_given)
			densities[j] *= data.weights[indices_of_input_points_to_be_used.at(j)];
		else
			densities[j] *= data.weights[j];
	}
	
	sum = densities.sum();
	if (sum>0)
		densities /= sum;
	
	return densities;
}


std::vector<idx_type> gmmutil::gaussPartition(Matrix const& points, Parameters const& gmmdesc)
{
	const idx_type n = points.cols();
	const idx_type k = gmmdesc.means.cols();
	std::vector<idx_type> partition;
	// initial solution: component 0
	Vector highestDensities = gmmutil::gaussianDensity(points, gmmdesc.means.col(0), gmmdesc.covariances[0]);
	for(idx_type i=0; i<n; ++i)
		partition.push_back((idx_type)0);
	Vector tmpDensities = Vector(n);
	// iterate over components 1 to k-1
	for(idx_type j=1; j<k; ++j)
	{
		tmpDensities = gmmutil::gaussianDensity(points, gmmdesc.means.col(j), gmmdesc.covariances[j]);
		for(idx_type i=0; i<n; ++i)
		{
			if(highestDensities(i) < tmpDensities(i))
			{
				partition.at(i) = j;
				highestDensities(i) = tmpDensities(i);
			}
		}
	}
	return partition;
}


std::vector<idx_type> gmmutil::gmmPartition(Matrix const& pmatrix)
{
	const idx_type n = pmatrix.cols();
	const idx_type k = pmatrix.rows();
	
	if(n==0||k==0)
		gmmlab_throw("gmmutil::gmmPartition() - pmatrix must not be empty.");
	
	std::vector<idx_type> partition;
	idx_type index;
	for(idx_type i=0; i<n; ++i)
	{
		pmatrix.col(i).maxCoeff(&index);
		partition.push_back(index);
	}
	
	return partition;
}


Matrix gmmutil::pmatrix(commonutil::DataSet const& input, Parameters gmmdesc)
{
	const idx_type k = gmmdesc.means.cols();
	const idx_type n = input.points.cols();
	Matrix pmatrix(k,n);
	for (idx_type i=0; i<k; ++i)
		pmatrix.row(i) = gmmdesc.weights[i]*gmmutil::gaussianDensity(input.points, gmmdesc.means.col(i), gmmdesc.covariances[i]);
	for (idx_type j=0; j<n; ++j)
		pmatrix.col(j) *= input.weights[j];
	return pmatrix;
}


Parameters gmmutil::optimalGaussian(commonutil::DataSet const& input, bool only_spherical, std::vector<idx_type> const& indices_of_input_points_to_be_used)
{
	bool indices_given = !indices_of_input_points_to_be_used.empty();
	
	const idx_type n = indices_given ? indices_of_input_points_to_be_used.size() : input.points.cols();
	const idx_type d = input.points.rows();
	
	fp_type totalWeight =  0;
	for(idx_type i=0; i<n; ++i)
		totalWeight += input.weights( indices_given ? indices_of_input_points_to_be_used.at(i) : i );
	
	Parameters output;
	
	output.means = Matrix::Zero(d,1);	
	for(idx_type i=0; i<n; ++i)
	{
		idx_type index =  indices_given ? indices_of_input_points_to_be_used.at(i) : i;
		output.means.col(0) += input.weights(index)*input.points.col(index);
	}
	output.means.col(0) = output.means.col(0) / totalWeight;
	
	Matrix covar = Matrix::Zero(d,d);
	fp_type spherical = 0;
	for(idx_type i=0; i<n; ++i)
	{
		idx_type index =  indices_given ? indices_of_input_points_to_be_used.at(i) : i;
		Vector diff = input.points.col(index)-output.means.col(0);
		if(!only_spherical)
			covar.noalias() += (input.weights(index)*diff)*diff.transpose();
		spherical += input.weights(index)*diff.squaredNorm();
	}
	covar = covar/totalWeight;
	
	if (only_spherical || !linalg::spd(covar))
	{
		fp_type factor = spherical/(totalWeight*d);
		if(std::isfinite(factor) && factor>0)
			covar = Matrix::Identity(d,d)*(spherical/(totalWeight*d));
		else
			covar = Matrix::Identity(d,d);
		
	}
	
	output.covariances.push_back(covar);
	output.weights = Vector(1);
	output.weights(0) = 1.;
	return output;
}



Matrix gmmutil::partialEM(commonutil::DataSet const& input, Parameters& gmmdesc, idx_type numSteps, bool verbose)
{
	// 	std::cout << "partialEM call " <<std::endl;
	
	const idx_type k = gmmdesc.weights.size();
	assert(gmmdesc.means.cols()==k && gmmdesc.covariances.size()==k);
	
	if (k==0)
		gmmlab_throw("No or empty initial solution given.");
	
	const idx_type n = input.points.cols();
	const idx_type d = input.points.rows();
	
	if (n==0 || d==0)
		gmmlab_throw("Input is empty.");
	
	fp_type totalWeight = input.weights.sum();
	if (totalWeight<=0)
		gmmlab_throw("Total input weight less or equal to zero.");
	
	fp_type minWeight = totalWeight;
	for (idx_type i=0; i<n; ++i)
		if (input.weights[i]>0 && input.weights[i]<minWeight)
			minWeight = input.weights[i];
		
		// probabilities/responsibilities regarding the first k-1 components
		Matrix pmatrix(k,n);
	for (idx_type i=0; i<k-1; ++i)
		pmatrix.row(i) = gmmdesc.weights[i]*gmmutil::gaussianDensity(input.points,
									     gmmdesc.means.col(i), gmmdesc.covariances[i]);
		for (idx_type j=0; j<n; ++j)
			pmatrix.col(j) *= input.weights[j];
		pmatrix.row(k-1) = Vector::Zero(n);
	
	// sum of the probabilities regarding the first k-1 components
	Vector pmatrixSums(n);
	for(idx_type j=0; j<n; ++j)
		pmatrixSums(j) = pmatrix.col(j).sum();
	fp_type pmatrixSum = pmatrixSums.sum();
	
	// weight of the k-th component
	fp_type alpha = gmmdesc.weights(k-1);
	
	for (unsigned int step=0; step<numSteps; ++step)
	{
		// only update the probabilities regarding the k-th component
		pmatrix.row(k-1) = alpha*gmmutil::gaussianDensity(input.points, gmmdesc.means.col(k-1), gmmdesc.covariances[k-1]);
		for (idx_type j=0; j<n; ++j)
			pmatrix(k-1, j) *= input.weights[j];
		
		// update probabilities
		Vector normalizedKthRow = Vector(n);
		for(idx_type j=0; j<n; ++j)
		{
			fp_type sum = alpha*pmatrix(k-1,j)+(1.-alpha)*pmatrixSums(j);
			if(sum<=0)
			{
				if (verbose)
					std::cout << "gmmutil::partialEM() - In round "<< step << " no gaussian is responsible for point " << j << std::endl;
				normalizedKthRow(j) = alpha;
			}
			else
				normalizedKthRow(j) = pmatrix(k-1,j)/sum;
		}
		
		// partition
		std::vector<idx_type> partition = gmmutil::gmmPartition(pmatrix);
		fp_type resp = 0;
		for(idx_type j=0; j<n; ++j)
			if(partition.at(j) == (k-1))
				resp += normalizedKthRow(j);
			
		if(resp > 0.)
		{
			// mean
			Vector tmpMean = Vector::Zero(d);	
			for(idx_type j=0; j<n; ++j)	
				if(partition.at(j) == k-1)
					tmpMean += normalizedKthRow(j) * input.points.col(j);
				tmpMean /= resp;
			// covar
			Matrix tmpCovar = Matrix::Zero(d,d);
			for(idx_type j=0; j<n; ++j)
				if(partition.at(j) == k-1) 
				{
					Vector y = input.points.col(j)-tmpMean;
					tmpCovar += (normalizedKthRow(j) * y) * y.transpose();
				}
				tmpCovar = tmpCovar / resp;
			
			// if the resulting covar is not positive definite, return the last valid gmm
			if(!linalg::spd(tmpCovar))
			{
				if (verbose)
					std::cout << "gmmutil::partialEM() - In round "<< step << " the update of "<< (k-1) <<"-th covariance resulted in semi-definite matrix. "
					<< "Stop and return last valid gmm."<< std::endl;
				return pmatrix;
			}
			
			// ... otherwise store the update
			alpha = resp/(pmatrixSum+pmatrix.row(k-1).sum()); // instead of dividing by "N" as in the paper, which does not necessarily result in alpha in [0,1]
			gmmdesc.means.col(k-1) = tmpMean;
			gmmdesc.covariances[k-1] = tmpCovar;
			
			
		}
		else
		{
			if (verbose)
				std::cout << "gmmutil::partialEM() - In round "<< step << " the "<< (k-1) <<"-th component is (in expectation) not responsible for any point.\n"
				<< "gmmutil::partialEM() - Stop and return last valid gmm. " << std::endl;
			
			// if the resulting component is not responsible for any point, return the last gmm 
			return pmatrix;
		}	
			
	}
	
	// adjust the weights!
	gmmdesc.weights *= (1.-alpha);
	gmmdesc.weights(k-1) = alpha;
	
	// update pmatrix
	pmatrix *= (1.-alpha);
	pmatrix.row(k-1) = alpha*gmmutil::gaussianDensity(input.points, gmmdesc.means.col(k-1), gmmdesc.covariances[k-1]);
	for (idx_type j=0; j<n; ++j)
		pmatrix(k-1, j) *= input.weights[j];
	
	// 	std::cout << "returned pmatrix = " << pmatrix.rows() << ", " << pmatrix.cols() << std::endl;
	
	return pmatrix;
}


idx_type gmmutil::minNLL(commonutil::DataSet const& data, std::vector<Parameters> const& models)
{
	const idx_type num = models.size();
	Vector nlls(num);
	for (idx_type i=0; i<num; ++i)
		nlls[i] = nll(data, models[i]);
	idx_type index;
	nlls.minCoeff(&index);
	return index;
}

Parameters gmmutil::appendGMMDesc(const Parameters& first, Parameters second, fp_type const weight)
{

//std::cout << "appendGMMDescs: Merge " << first.components() << " and " << second.components() << std::endl;

	Parameters merged;
	idx_type first_k = first.components();
	idx_type second_k = second.components();
		
	if(first_k == 0)
		
		merged = second;
	
	else if(second_k > 0)
	{
		merged = first;
		idx_type d = second.means.rows();
		second.weights *= weight;
		merged.weights.conservativeResize(first_k+second_k);
		merged.means.conservativeResize(d,first_k+second_k);
		for(idx_type i=0; i<second_k; ++i)
		{
			merged.weights(first_k+i) = second.weights(i);
			merged.means.col(first_k+i) = second.means.col(i);
			merged.covariances.push_back(second.covariances.at(i));
		}
		merged.weights /= merged.weights.sum();
	}
	
//std::cout << "appendGMMDescs: => " << merged.components() << std::endl;
	
	return merged;
}

Parameters gmmutil::removeComponent(const Parameters& gmmdesc, idx_type k, bool redistributeWeight)
{

	idx_type with_k = gmmdesc.components();
	
	if(k<0 || k>=with_k)
		gmmlab_throw("k out of range.");
	if(gmmdesc.weights.size()==0)
		gmmlab_throw("Empty gmmdesc given.");
	
	idx_type d = gmmdesc.means.rows();
	
	// weight which has to be destributed among the remaining components
	fp_type weight = gmmdesc.weights(k);
	
	// remove the k-th component
	Parameters without;
	without.means.resize(d,with_k-1);
	without.weights.resize(with_k-1);
	for(idx_type i=0; i<with_k-1; ++i)
	{
		idx_type index = (i<k) ? i : i+1;
		without.weights(i) = gmmdesc.weights(index);
		without.means.col(i) = gmmdesc.means.col(index);
		without.covariances.push_back(gmmdesc.covariances.at(index));
	}
	
	if(redistributeWeight)
	{
		// normalize
		without.weights /= without.weights.sum();
		
		// distribute released weight ...
		for(idx_type i=0; i<with_k-1; ++i)
			without.weights(i) += without.weights(i) * weight;
		
		// ... and normalize again
		without.weights /= without.weights.sum();
	}
	
	return without;
}

Parameters gmmutil::getComponent(const Parameters& gmmdesc, idx_type k)
{
	if(k<0 || k >= gmmdesc.components())
		gmmlab_throw("k out of range");
	idx_type d = gmmdesc.means.rows();
	Parameters component;
	component.means.resize(d,1);
	component.means.col(0) = gmmdesc.means.col(k);
	component.weights.resize(1);
	component.weights(0) = gmmdesc.weights(k);
	component.covariances.push_back(gmmdesc.covariances.at(k));
	return component;
}

Parameters gmmutil::mergeComponents(Parameters const& gmmdesc, idx_type k, idx_type l)
{
	if(gmmdesc.components() <= k || gmmdesc.components() <= l || k<0 || l<0)
		gmmlab_throw("Given gmmdesc has no such components.");
	
	if(k > l)
	{
		idx_type tmp = k;
		k = l;
		l = tmp;
	}
	
	Parameters component = gmmutil::getComponent(gmmdesc, l);
	
	Parameters merged = gmmutil::removeComponent(gmmdesc, l, false);
	
//std::cout << "MERGECOMPONENTS: sum weights = " << gmmdesc.weights.sum() << std::endl;
	
	fp_type sum_of_weights = merged.weights(k) + component.weights(0);
	Vector merged_mu = merged.means.col(k);
	merged.means.col(k) = (merged.weights(k) * merged.means.col(k) + component.weights(0) * component.means.col(0))/sum_of_weights;
	Vector diff_merged = (merged_mu - merged.means.col(k));
	Vector diff_component = (component.means.col(0) - merged.means.col(k));
	merged.covariances.at(k) = merged.weights(k) * ( merged.covariances.at(k) + diff_merged * diff_merged.transpose() );
	merged.covariances.at(k).noalias() += component.weights(0) * ( component.covariances.at(0) + diff_component * diff_component.transpose() );	
	merged.covariances.at(k) /= sum_of_weights; 
	merged.weights(k) = sum_of_weights;
	
	return merged;
}


fp_type gmmutil::getSeparation(const Vector& mean1, const fp_type trace1, const Vector& mean2, const fp_type trace2)
{
	return(mean1-mean2).norm() / std::sqrt(std::max(trace1, trace2));
}

fp_type gmmutil::getSeparation(Parameters const& gmmdesc)
{
	idx_type k = gmmdesc.means.cols();
	idx_type sqrtd = sqrt(gmmdesc.means.rows());
	
	//Vector maxSqrtEigenvalues = Vector::Zero(k);
	Vector traces = Vector::Zero(k);
	fp_type separation;
	
	if(k==1)
		return 1.;
	
	for(idx_type i=0; i<k; ++i)
	{
		//Eigen::SelfAdjointEigenSolver<Matrix> eigensolver(gmmdesc.covariances.at(i));
		//if (eigensolver.info() != Eigen::Success)
		//	gmmlab_throw("Eigensolver failed.");
	    //maxSqrtEigenvalues(i) = std::sqrt(eigensolver.eigenvalues().maxCoeff());
		traces(i) = gmmdesc.covariances.at(i).trace();
		
		for(idx_type j=0; j<i; ++j)
		{
			//fp_type tmp = (gmmdesc.means.col(i)-gmmdesc.means.col(j)).norm() / (sqrtd * std::max(maxSqrtEigenvalues(i),maxSqrtEigenvalues(j)));
			fp_type tmp = gmmutil::getSeparation(gmmdesc.means.col(i), traces(i), gmmdesc.means.col(j), traces(j));
			if( (i==1&&j==0) || tmp < separation)
				separation = tmp;
		}
	}
	return separation;	
}


