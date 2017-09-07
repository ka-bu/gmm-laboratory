#ifndef DATAGENUTIL_H
#define DATAGENUTIL_H

#include <boost/concept_check.hpp>
#include "../algorithm/gmm/utils/gmmutil.h"

namespace datageneration
{
	
	static const std::string ME_NO_ME = "ME_0";
	static const std::string ME_BLUR_CORRELATED = "ME_BC";
	static const std::string ME_BLUR_DIFFICULT = "ME_BD";
	static const std::string ME_UNIFORM_NOISE = "ME_UN";
	static const std::string ME_GAUSSIAN_NOISE = "ME_G";
	
	struct BoundingBox
	{
		Vector min;
		Vector max;
	};
	
	/**
	 * generates a random Gaussian mixture model.
	 */
	Parameters generateRandomGMM(idx_type d, idx_type k, std::mt19937& gen, unsigned int weightExp, fp_type spread = 1);
	
	/**
	 * generates a random mixture model.
	 * For each generated covariance matrix this method ensures that the proportion between its maximum and minimum sqrt(eigenvalue) is in [minSqrtEWProportion,maxSqrtEWProportion]
	 * and that the minimum eigenvalue is in [minSqrtEW, maxSqrtEW]. 
	 * The entries of the means are uniformly distributed in a cube and then scaled such that the separation is equal to the given separation value.
	 */
	Parameters generateRandomGMMWithUniformMeans(idx_type d, idx_type k, std::mt19937& gen,
						   fp_type separation,
						   fp_type minSqrtEW, fp_type maxSqrtEW, 
						   fp_type minSqrtEWProportion, fp_type maxSqrtEWProportion,
						   fp_type weightExp,
						   fp_type minSqrtEWExp = -1);
		
	/**
	 * generates a random mixture model.
	 * For each generated covariance matrix this method ensures that the proportion between its maximum and minimum sqrt(eigenvalue) is in [minSqrtEWProportion,maxSqrtEWProportion]
	 * and that the minimum eigenvalue is in [minSqrtEW, maxSqrtEW]. 
	 * The entries of the means are uniformly distributed in [minMean, maxMean].
	 */
	Parameters generateRandomGMMWithUniformMeans(idx_type d, idx_type k, std::mt19937& gen,
						     fp_type minMean, fp_type maxMean,
						     fp_type minSqrtEW, fp_type maxSqrtEW, 
						     fp_type minSqrtEWProportion, fp_type maxSqrtEWProportion,
						     fp_type weightExp,
						     fp_type minSqrtEWExp);
	

	
	/**
	 * generates a dataset from according to a given mixture model.
	 */
	void generateInputFromGMM(commonutil::DataSet&, Parameters const&, idx_type, std::mt19937& gen);

	void generateInputFromMPEMM(commonutil::DataSet&, Parameters const&, idx_type, std::mt19937& gen);
			
	
	/**
	 * noise of a data point correlates with the covariance matrix of the component that the data point has been drawn from.
	 */
	void generateCorrelatedNoisyInputFromGMM(commonutil::DataSet& target, Parameters const& desc, idx_type n, fp_type sqrtEWFraction, std::mt19937& gen);
	
	
	/**
	 * noise of a data point is determined by a gaussian with spherical covariance matrix. 
	 * The eigenvalue of this matrix is drawn at random.
	 * Let min/max the minimum/maximum over all maximum eigenvalues of the covariance matrices of the truth.
	 * Then the eigenvalue is chosen to be min + i/k*(max-min) with probability proportional to 2^i, where k is the number of components of the truth.
	 */
	void generateDifficultNoisyInputFromGMM(commonutil::DataSet& target, Parameters const& desc, idx_type n, std::mt19937& gen);
	
	/**
	 * adds noise points uniformly at random in an enlarged bounding box of the given data.
	 * The enlarged bounding box has sidelengths which are 1.2 times as large as the side lengths of the actual bounding box
	 * and has the same center.
	 */
	void addUniformNoise(commonutil::DataSet& dataset, idx_type n, std::mt19937& gen);
	
	/**
	 * draws points uniformly at random from the bounding box of the given points and projects them to a random subspace with dimension >=1 and <=(d-1), the dimension of the points in the given data set.
	 */
	void addUniformProjectedToSubspaceNoise(commonutil::DataSet& dataset, idx_type n, std::mt19937& gen, bool verbose);
	
	/**
	 * adds points that are randomly drawn according to k "shape"-distribution by calling one of datagenutil::addBox, datagenutil::addArc,... .
	 * Each of the implemented shapes is drawn with equal probability.
	 */
	void addUniformShapes(commonutil::DataSet& dataset, idx_type d, idx_type n, std::mt19937& gen, idx_type k, fp_type genWeightFloatingExp,
		           fp_type minMean, fp_type maxMean, fp_type minDiameter, fp_type maxDiameter);
	
	/**
	 * draws point uniformly at random from 6 boxes that are arranged such that they form a bridge.
	 * Each of these subboxes contain 1/6 of the total number of points. The whole bridge fits in the given interval.
	 */
	void addUniformBridge(commonutil::DataSet& dataset, idx_type d, idx_type n, std::mt19937& gen, Vector anchor, Vector const& interval);
	
	/**
	 * draws points randomly as follows:
	 * 1. draws 2 dimensions uniformly at random
	 * 2. for these 2 dimension: it proceeds just at datagenutil::addSphere()
	 * 3. for the remaining dimensions: draws entries uniformly at random from the interval [minRadius, maxRadius]
	 */
	void addUniformArc(commonutil::DataSet& dataset, idx_type d, idx_type n, std::mt19937& gen, Vector center, fp_type minRadius, fp_type maxRadius);
	
	/**
	 * draws points uniformly at random from the box that is defined by anchor+{0,1}^d*interval, where d is the dimension of the interval and center.
	 */
	void addUniformBox(commonutil::DataSet& dataset, idx_type d, idx_type n, std::mt19937& gen,  Vector const& anchor, Vector const& interval);
	
	/**
	 * draws points randomly as follows:
	 * 1. draws a random radius from [minRadius, maxRadius]
	 * 2. draws a point uniformly at random from the sphere with the drawn radius and zero-vector as center
	 * 2. shifts the drawn points by the center-vector
	 */
	void addUniformSphere(commonutil::DataSet& dataset, idx_type d, idx_type n, std::mt19937& gen, Vector center, fp_type minRadius, fp_type maxRadius);
	
	/**
	 * returns the minimum and maximum entries that points in the given dataset have.
	 */
	BoundingBox boundingBox(commonutil::DataSet const& dataset);
	
	/**
	 * multiplies the vector returned by datagenutil::getExpIncrWeights() times the given number and rounds the resulting entries.
	 * It is ensured that all vector entries sum up to the number by setting the last entry of the vector to  number - sum of the remaining entries.
	 */
	Vector getExpIncrNumbers(idx_type k, fp_type exponent, idx_type number);
	
	/**
	 * returnes a vector (2^0, 2^(exponent), 2^(2*exponent),...,2^(k*exponent)) normalized by the sum of its entries.
	 * I.e., the entries of the returned vector are exponentially increasing, are in [0,1], and sum up to 1.
	 */
	Vector getExpIncrWeights(idx_type k, fp_type exponent);
};

#endif // GENUTIL_H
