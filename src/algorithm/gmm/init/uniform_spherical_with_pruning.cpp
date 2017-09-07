#include "uniform_spherical_with_pruning.h"

#include "uniform_spherical.h"
#include "../../../base.h"

const std::string UniformSphericalWithPruningID::CLASSTAG = "UniformSphericalWithPruning";

UniformSphericalWithPruningID::UniformSphericalWithPruningID(unsigned int oversamplingFactor, uint32_t s) : seed(s), oversamplingFactor(oversamplingFactor)
{
	std::stringstream sstream;
	sstream << "UniformSphericalWithPruning " << "_o" << oversamplingFactor << "_i" << s ;
	this->nametag = sstream.str();
}

void UniformSphericalWithPruningID::parseDescriptionString(std::vector< InitDescription* >& descriptions, std::string descriptionstring)
{
	std::istringstream iss(descriptionstring);
	std::vector<unsigned int> oList = readIntList(iss);
	std::vector<unsigned int> sList = Description::readIntList(iss);
	for(idx_type s=0; s<sList.size(); ++s)
		for(idx_type ok = 0; ok<oList.size(); ++ok)
			descriptions.push_back(new UniformSphericalWithPruningID(oList[ok], sList[s]));
}


Parameters UniformSphericalWithPruningID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	std::mt19937 gen(seed);
	return initializer::uniformSphericalGMMwithPruning(input, k, gen, oversamplingFactor);
}


Parameters initializer::uniformSphericalGMMwithPruning(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, unsigned int oversamplingFactor)
{
	initializer::check(k, input.points.rows(), input.points.cols());
	
	idx_type oversamplingk = (idx_type)(round(oversamplingFactor*1.*k*log(oversamplingFactor*1.*k)));
		
	if(oversamplingk < k)
	  gmmlab_throw("initutil::uniformSphericalGMMwithPruning() - oversamplingk should be larger than k.");
	
	Parameters oversampled = initializer::uniformSphericalGMM(input, oversamplingk, gen);
	return initializer::prune(oversampled, k, gen);
}

Parameters initializer::prune(Parameters& desc, idx_type k, std::mt19937& gen)
{
	if(desc.weights.size() == k)
		return desc;
	if(desc.weights.size() < k)
		gmmlab_throw("gmm has less than k components");
	
	idx_type d = desc.means.rows();

	// 1. remove center estimates whose mixing weights are below w_T = 1/(4k) 
	idx_type desc_k = desc.weights.size();
	idx_type index;
	while(desc.weights.size()>k)
	{
		if(desc.weights.minCoeff(&index) < 1/(4.*desc_k))
		{
			commonutil::erase(desc.weights,index);
			desc.covariances.erase(desc.covariances.begin()+index);
			commonutil::eraseColumn(desc.means, index);
		}
		else
			break;
	}

	desc_k = desc.weights.size();
	if (desc_k>k)
	{
		Parameters prunedGMM;
		prunedGMM.weights = Vector::Constant(k, 1.0/k);
		prunedGMM.means = Matrix::Zero(d, k);
		prunedGMM.covariances.clear();

		// 2. prune the center estimates (adaption of Gonzales (1985))
		//    and set the mixing weights to w_i = 1/k

		// 2.a Choose one of these centers arbitrarily
		std::uniform_int_distribution<> uid(0, desc.weights.size()-1);
		index = uid(gen);

		prunedGMM.means.col(0) = desc.means.col(index);
		prunedGMM.covariances.push_back(desc.covariances.at(index));

		// 2.b Pick the remaining k-1 centers iteratively: pick the center farthest from the ones picked so far,
		//     i.e., pick he point x with $min_{y\in S} dist(x,y)$
		fp_type diffcovar, dist;
		for(idx_type i=1; i<k; ++i)
		{
			Vector distances = Vector(desc_k);
			for (idx_type j=0; j<desc_k; ++j)
			{
				diffcovar = desc.covariances.at(j)(0,0) - prunedGMM.covariances.at(0)(0,0);
				if(diffcovar <= 0)
				  diffcovar = FP_INFINITE;
				  //gmmlab_throw("difference of the covariances is <= 0.");
				distances[j] = (desc.means.col(j)-prunedGMM.means.col(0)).norm()/diffcovar;
				
				for (idx_type l=1; l<i; ++l)
				{
				  	diffcovar = desc.covariances.at(j)(0,0) - prunedGMM.covariances.at(l)(0,0);
					if(diffcovar <= 0)
					  diffcovar = FP_INFINITE;
					  //gmmlab_throw("difference of the covariances is <= 0.");
					dist = (desc.means.col(j)-prunedGMM.means.col(l)).norm()/diffcovar;

					if (dist<distances[j])
						distances[j]=dist;
				}
			}

			// pick the farthest
			distances.maxCoeff(&index);

			// insert into pruned GMM
			prunedGMM.means.col(i)=desc.means.col(index);
			prunedGMM.covariances.push_back(desc.covariances.at(index));
		}

		
		//std::cout << "prunedGMM = \n" << prunedGMM << std::endl;
		
		return prunedGMM;
	}
	return desc;
}