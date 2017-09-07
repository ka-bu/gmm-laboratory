#include "lloyds.h"

#include <iostream>
#include <ctime>
#include <iostream>

#include "../../base/linalgutil.h"
#include "../gmm/utils/gmmutil.h"

#include "../../settings/settings.h"
#include "utils/kmeansutil.h"

const std::string LloydsAD::CLASSTAG = "Lloyds";

Algorithm* LloydsAD::create(const commonutil::DataSet& input) const
{
	return new Lloyds(input, commonSettings().verbose, this->seed, this->ratio);
}

LloydsAD::LloydsAD(unsigned int st, fp_type ratio, uint32_t seed):seed(seed),ratio(ratio)
{
	this->steps = st;
	std::stringstream sstream;
	sstream << "Lloyds";
	if(ratio>0)
		sstream << "_r" << ratio << "_sd" << seed;
	sstream << "_s" << st;
	this->nametag = sstream.str();
}

void LloydsAD::parseDescriptionString(std::vector< AlgoDescription* >& descriptions, const std::string descriptionstring, const bool first)
{
	std::istringstream iss(descriptionstring);
	std::vector<unsigned int> stList = readIntList(iss);
	std::vector<fp_type> ratioList = readFpTypeList(iss);
	std::vector<unsigned int> sdList = readIntList(iss);
	
	
	if (first)	
		for (std::size_t st=0; st<stList.size(); ++st)
			for (std::size_t sd=0; sd<sdList.size(); ++sd)
				for (std::size_t r=0; r<ratioList.size(); ++r)
					descriptions.push_back(new LloydsAD(stList[st], ratioList[r], sdList[sd]));
	else
	{
		if(stList.size()!=1)
			gmmlab_throw("multiple step values forbidden for next:...!");
		if(sdList.size()!=1)
			gmmlab_throw("multiple seed values forbidden for next:...!");
		if(ratioList.size()!=1)
				gmmlab_throw("multiple ratio values forbidden for next:...!");
		descriptions.push_back(new LloydsAD(stList[0], ratioList[0], sdList[0]));
	}
}

// =================================================================================================================================


Lloyds::Lloyds(commonutil::DataSet const& input, bool verbose, uint32_t seed,  fp_type ratio)
	: RandomAlgorithm(input, verbose, seed), ratio(ratio)
{
}

Lloyds::Lloyds(commonutil::DataSet const& input, bool verbose, std::mt19937& gen,  fp_type ratio)
	: RandomAlgorithm(input, verbose, 0), ratio(ratio)
{
	this->gen = gen;
}


void Lloyds::init(const Matrix& means)
{
	Parameters gmmdesc;
	gmmdesc.means = means;
	Algorithm::init(gmmdesc);
}


void Lloyds::run(unsigned int numSteps)
{

	clock_t start, end;
	start = clock();

	std::size_t k = this->desc.means.cols();

	if (k==0)
		gmmlab_throw("No or empty initial solution given.");

	idx_type n = this->input->points.cols();
	idx_type d = this->input->points.rows();

	if (n==0 || d==0)
		gmmlab_throw("Input is empty.");

	fp_type totalWeight = this->input->weights.sum();
	if (totalWeight<=0)
		gmmlab_throw("Total input weight less or equal to zero.");

	std::uniform_real_distribution<> uid(0, 1);
	
	std::vector<idx_type> indices(n);
	
	for (unsigned int step=0; step<numSteps; ++step)
	{
	  		
		// assign the points to the centers (either probabilistically or deterministically)
		bool probabilistic = uid(this->gen) < this->ratio;
		
		if(probabilistic)
		{
		  
		      Matrix pmatrix(k,n);
		  
		      std::cout << "ratio = " << this->ratio << std::endl;
		      gmmlab_throw("Not implemented properly yet.");
		      
		      // compute the squared distance between each point and each center
		      // i.e., pmatrix(i,j)=distance of the j-th point to the i-th center
		      for (idx_type i=0; i<k; ++i)
			      pmatrix.row(i) = kmeansutil::squaredDistances(this->input->points,this->desc.means.col(i));

		      //for(idx_type i=0; i<pmatrix.cols(); ++i)
		      // pmatrix.col(i) = pmatrix.col(i)/pmatrix.col(i).sum();
		      
		      for (idx_type i=0; i<n; ++i)
		      {
			      indices[i] = commonutil::randomIndex(pmatrix.col(i), this->gen);			
			      assert(indices[i]<k);
		      }
		 
		}
		else
		{
			for (idx_type i=0; i<n; ++i)
				// note (...).minCoeff(&indices[i])
				kmeansutil::squaredDistances(this->desc.means, this->input->points.col(i)).minCoeff(&(indices[i]));
		}

		this->desc.weights = Vector::Zero(k);
		this->desc.means = Matrix::Zero(d,k);

		// for each cluster: sum up the weights and the points
		for (std::size_t i=0; i<n; ++i)
		{
			this->desc.weights[indices[i]] += this->input->weights[i];
			this->desc.means.col(indices[i]).noalias() += this->input->weights[i]*this->input->points.col(i);
		}		
		
		// for each cluster: divide the sum of points by the weight of the cluster.
		// if the weight of a cluster is zero (i.e., there are no points in this cluster), then choose a new center
		// uniformly at random
		for (std::size_t i=0; i<k; ++i)
			if (this->desc.weights[i] > 0)
				this->desc.means.col(i) = this->desc.means.col(i)/this->desc.weights[i];
			else
			{
				this->desc.means.col(i) = this->input->points.col(commonutil::randomIndex(this->input->weights, this->gen));
 				if (commonSettings().verbose)
 					std::cout << "Lloyds::run() - replaced empty mean " << i+1 << " by random input point" << std::endl;
			}

		// if the result shall be displayed, then compute the covariance and weights of the resulting clusters
		if(!commonSettings().ndisplay)
			this->desc = kmeansutil::wrapMeans(*this->input, this->desc.means, !commonSettings().ndisplay);
	}

	this->change = true;

	end = clock();
	this->runtime += double(end-start)/CLOCKS_PER_SEC;
}

Lloyds* Lloyds::toLloyds(Algorithm* a)
{
	return dynamic_cast<Lloyds*>(a);
}

