#ifndef GMMUTIL_H
#define GMMUTIL_H

#include "../../../base.h"
#include "../../../base/commonutil.h"
#include "../../../base/parameters.h"

namespace gmmutil
{

	/**
	 * creates a GMM based on a set of k means and the clustering induced by these means 
	 * (i.e. each point is assigned to its nearest mean) as follows:
	 * - The mean of a component is set to the centroid of the corresponding  cluster.
	 * - The covariance of a component is set to the optimal covariance of the (k=1)-problem
	 *   that is given by the corresponding cluster and its mean.
	 * - The weights are estimated by the fraction of points assigned to the corresponding
	 *   mean.
	 */
	Parameters meansToGMM(commonutil::DataSet const& input, Matrix const& means, bool spherical, bool verbose, std::vector<idx_type> const& indices_of_input_points_to_be_used = std::vector<idx_type>());
	
	/**
	 * returns the minimum squared (pairwise) distance between the given points.
	 */
	const fp_type minSquaredDistance(Matrix const&);

	/**
	 * returns the density of each of the points as defined by the gaussian distribution (which is defined by the given mean and covariance).
	 */
	Vector gaussianDensity(Matrix const& points, Vector const& mean, Matrix covariance, const std::vector<idx_type>& indices_of_input_points_to_be_used = std::vector<idx_type>());
	
	/**
	 * returns the densitiy of each of the points as defined by the given mixture model.
	 */
	Vector gmmDensity(Matrix const& points, Parameters const& desc, const std::vector<idx_type>& indices_of_input_points_to_be_used = std::vector<idx_type>());
	
	/**
	 * returns the negative-log-likelihood of each of the points w.r.t. a gaussian distribution.
	 */
	Vector gaussianNLL(Matrix const& points, Vector const& mean, Matrix variance, const std::vector<idx_type>& indices_of_input_points_to_be_used = std::vector<idx_type>());
	
	/**
	 * returns the negative-log-likelihood of each of the points w.r.t. the gaussian mixture model.
	 */
	Vector gmmNLL(Matrix const& points, Parameters const& gmmdesc, const std::vector<idx_type>& indices_of_input_points_to_be_used = std::vector<idx_type>());
	
	
	fp_type nll(commonutil::DataSet const& input, Parameters const& gmmdesc, const std::vector<idx_type>& indices_of_input_points_to_be_used = std::vector<idx_type>());

	fp_type nlcdl(commonutil::DataSet const& data, Parameters const& desc);
	
	/**
	 * returns the mahalanobis distance between the given gaussian and each of the points.
	 */
	Vector gaussianMahalanobisOfInvSqrtCovar(Matrix const& points, Vector const& mean, Matrix const& sqrt_of_inverted_covariance);
		
	/**
	 * returns the mahalanobis distance between the given gaussian and each of the points with index i, where i is contained in cur_dataset.
	 */
	Vector gaussianMahalanobis(Matrix const& points, Vector const& mean, Matrix covariance, const std::vector<idx_type>& indices_of_input_points_to_be_used = std::vector<idx_type>());
	
	/**
	 * returns the density of each of the points, which is given by a mixture of Mahalanbos distances.
	 * The weights of this mixture are the weights of the GMM, while the means and covariances define the
	 * single Mahalanobis distances.
	 */
	Vector gmmMixMahalanobis(Matrix const& points, Parameters const& gmmdesc);
	

	/**
	 * returns the minimal Mahalanobis distance to a mean w.r.t. the respective covariance matrix
	 * for each of the given points.
	 */
	Vector gmmMinMahalanobis(Matrix const& points, Parameters const& gmmdesc, const std::vector<idx_type>& indices_of_input_points_to_be_used = std::vector<idx_type>());
	

	
	idx_type indexByAdaptiveDensities(commonutil::DataSet const& data, Parameters const& desc, fp_type const ratio, std::mt19937& gen, const std::vector<idx_type>& indices_of_input_points_to_be_used = std::vector<idx_type>());
	

	/**
	 * returns densities defined by factor * gmmMinMahalanobis + (1-factor) * uniform.
	 */
	Vector adaptiveDensities(commonutil::DataSet const& input, Parameters const& gmmdesc, fp_type const ratio, const std::vector<idx_type>& indices_of_input_points_to_be_used = std::vector<idx_type>());

	/**
	 * returns the optimal (k=1)-solution. 
	 * Guarantees that the returned covariance matrix is positive definite: In case the computed optimal covariance is not positive definite it
	 * returns an approximated spherical covariance or the identity matrix.
	 */	
	Parameters optimalGaussian(commonutil::DataSet const& input, bool only_spherical = false, const std::vector<idx_type>& indices_of_input_points_to_be_used = std::vector<idx_type>());
	
	/**
	 * returnes a matrix whose (i,j)-th entry contains w(x_j) * w_i * N(x_j | mu_i, Sigma_i).
	 * 
	 */
	Matrix pmatrix(commonutil::DataSet const& input, Parameters gmmdesc);
	
	/**
	 * executes a partial EM that only updates the k-th component of the given gmm (cf. "Efficient Greedy Learning of Gaussian Mixture Models", Verbeek et al.)
	 * and returns the p-matrix w.r.t. the resulting gmmdesc.
	 */
	Matrix partialEM(commonutil::DataSet const& input, Parameters& gmmdesc, idx_type numSteps, bool verbose=0);
	
	
	/**
	 * returns the assignment of each point to its likeliest component, i.e. argmax_j=1..k{ w_j N(x | \mu_j,\sigma_j) }.
	 */
	std::vector<idx_type> gmmPartition(Matrix const& pmatrix);
	
	/**
	 * returns the assignment of each point to the component with highest density i.e. argmax_j=1..k{ N(x | \mu_j,\sigma_j) }.
	 */
	std::vector<idx_type> gaussPartition(Matrix const& points, Parameters const& gmmdesc);
	
	
// 	/**
// 	 * Given some means this method computes an empirical estimate of the variances and the weights.
// 	 * This is done by computing the clusters induced by the means (i.e., each point is assigned to its nearest cluster) and
// 	 */
// 	void empiricalCovarianceAndWeights(DataSet const& data, Vector const& means, std::vector<Matrix>& covariances, Vector& weights);

	idx_type minNLL(commonutil::DataSet const&, std::vector<Parameters> const&);
	
	/**
	 * appends the second gmm to the first gmm and modifies the weights as follows: First, all weights of the second gmm will
	 * be mulitplied with the given weight. Then, the weights of the complete gmm will be normalized.
	 */
	Parameters appendGMMDesc(Parameters const& first, Parameters second, fp_type const weight);
	
	/**
	 * removes the k-th component from the given gmmdesc.
	 */
	Parameters removeComponent(Parameters const& gmmdesc, idx_type k, bool redistributeWeight = true);
	
	/**
	 * merges the k-th component of the given gmmdesc with the given component (a single Gaussian).
	 */
	Parameters mergeComponents(Parameters const& merged, idx_type k, idx_type l);
	
	/**
	 * returns a GMMDesc containing only the k-th component of the given gmmdesc.
	 */
	Parameters getComponent(Parameters const& gmmdesc, idx_type k);
	
	/**
	 * computes the separation as defined by Dasgupta, i.e.
	 *   ||mu_i - mu_j|| / sqrt( max(trace(Sigma_i), trace(Sigma_j) ) ) .
	 */
	fp_type getSeparation(Vector const& mean1, fp_type const trace1, Vector const& mean2, fp_type const trace2);
	
	/**
	 * computes the separation as defined by Dasgupta, i.e.
	 * min_ij { ||mu_i - mu_j|| / sqrt( max(trace(Sigma_i), trace(Sigma_j) ) ) }.
	 */
	fp_type getSeparation(Parameters const& gmmdesc);
}


#endif
