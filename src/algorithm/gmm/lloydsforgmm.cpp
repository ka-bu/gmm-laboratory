#include "lloydsforgmm.h"

#include <iostream>
#include <ctime>
#include <cmath>

#include "../../base/linalgutil.h"
#include "../../settings/settings.h"
#include "utils/gmmutil.h"


const std::string LloydsForGMMAD::CLASSTAG = "LloydsForGMM";

LloydsForGMMAD::LloydsForGMMAD(unsigned int st, uint32_t sd, fp_type ratio, bool spherical ) : seed(sd), ratio(ratio), spherical(spherical)
{
	this->steps = st;
	std::stringstream sstream;
	sstream << "LloydsForGMM" << "_s" << st << "_r" << ratio << "_sph" << spherical << "_a" << sd;
	this->nametag = sstream.str();
}

void LloydsForGMMAD::parseDescriptionString(std::vector<AlgoDescription*>& descriptions, const std::string descriptionstring, const bool first)
{
	std::istringstream iss(descriptionstring);
	std::vector<unsigned int> stList = Description::readIntList(iss);
	fp_type llr = Description::readFpType(iss);
	bool spherical = (Description::readInt(iss)==1);
	std::vector<unsigned int> sdList = Description::readIntList(iss);
	
	if (first)
	{
		for (std::size_t st=0; st<stList.size(); ++st)
			for (std::size_t sd=0; sd<sdList.size(); ++sd)
				descriptions.push_back(new LloydsForGMMAD(stList[st], sdList[sd], llr, spherical));
	}
	else
	{
		if (stList.size()!=1)
			gmmlab_throw("multiple step values forbidden for next:...!");
		if (sdList.size()!=1)
			gmmlab_throw("multiple seed values forbidden for next:...!");
		descriptions.push_back(new LloydsForGMMAD(stList[0], sdList[0], llr, spherical));
		
	}
}

Algorithm* LloydsForGMMAD::create(commonutil::DataSet const& input) const
{
	return new LloydsForGMM(input, commonSettings().verbose, this->seed, this->ratio, this->spherical);
}

// =================================================================================================================================


LloydsForGMM::LloydsForGMM(commonutil::DataSet const& ds, bool verbose, uint32_t seed, fp_type ratio, bool spherical)
	: RandomAlgorithm(ds, verbose, seed), ratio(ratio), spherical(spherical)
{
}

LloydsForGMM::LloydsForGMM(const commonutil::DataSet& ds, bool verbose, std::mt19937& gen, fp_type ratio, bool spherical): RandomAlgorithm(ds, verbose, 0),  spherical(spherical), ratio(ratio)
{
	this->gen = gen;
	
	//std::cout << "ratio = " << this->ratio << std::endl;
}



void LloydsForGMM::init(Parameters const& desc)
{
	this->desc = desc;
	this->change = true;
}

void LloydsForGMM::run(unsigned int numSteps)
{
	clock_t start, end;
	start = clock();

	std::size_t k = this->desc.weights.size();

	assert(this->desc.means.cols()==k && this->desc.covariances.size()==k);

	if (k==0)
		gmmlab_throw("No or empty initial solution given.");

	idx_type n = this->input->points.cols();
	idx_type d = this->input->points.rows();

	if (n==0 || d==0)
		gmmlab_throw("Input is empty.");

	fp_type totalWeight = this->input->weights.sum();
	if (totalWeight<=0)
		gmmlab_throw("Total input weight less or equal to zero.");

	Matrix pmatrix(k,n);
	std::uniform_real_distribution<> urd(0, 1);
	for (unsigned int step=0; step<numSteps; ++step)
	{
		for (idx_type i=0; i<k; ++i)
			pmatrix.row(i) = this->desc.weights[i]*gmmutil::gaussianDensity(this->input->points,
				this->desc.means.col(i), this->desc.covariances[i]);

		for (idx_type j=0; j<n; ++j)
			if (pmatrix.col(j).sum()<=0)
				pmatrix.col(j) = this->desc.weights;

		bool probabilistic = false;
		if(this->ratio > 0)
			probabilistic = urd(this->gen)<this->ratio;

		std::vector<idx_type> indices(n);
		for (idx_type i=0; i<n; ++i)
		{
			if (probabilistic)
				indices[i] = commonutil::randomIndex(pmatrix.col(i), this->gen);
			else
				pmatrix.col(i).maxCoeff(&indices[i]);

			assert(indices[i]<k);
		}

		this->desc.weights = Vector::Zero(k);
		this->desc.means = Matrix::Zero(d,k);
		std::vector<Matrix> covs = this->desc.covariances;
		this->desc.covariances = std::vector<Matrix>(k, Matrix::Zero(d,d));

		for (std::size_t i=0; i<n; ++i)
		{
			this->desc.weights[indices[i]] += this->input->weights[i];
			this->desc.means.col(indices[i]).noalias() += this->input->weights[i]*this->input->points.col(i);
		}

		for (std::size_t i=0; i<k; ++i)
			if (this->desc.weights[i] > 0)
				this->desc.means.col(i) = this->desc.means.col(i)/this->desc.weights[i];
			else
			{
				this->desc.means.col(i) = this->input->points.col(commonutil::randomIndex(this->input->weights, this->gen));
				if (this->verbose)
					std::cout << "LloydsForGMM::run() - replaced empty mean " << i+1 << " by random input point" << std::endl;
			}

		if(this->spherical)
		{
			
			std::vector<fp_type> sphericals(k,0);
			
			for (std::size_t i=0; i<n; ++i)
			{
				Vector y = this->input->points.col(i) - desc.means.col(indices[i]);
				sphericals[indices[i]] += input->weights[i]*y.squaredNorm();
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
							std::cout << "LloydsForGMM::run() - replaced non-spd covariance " << i+1 << " with unit sphere covariance." << std::endl;
					}
					else
						if (verbose)
							std::cout << "LloydsForGMM::run() - replaced non-spd covariance " << i+1 << " with spherical covariance." << std::endl;
				}
				else
				{
					desc.weights[i] = totalWeight/k;
					desc.covariances[i] = Matrix::Identity(d,d);
					if (verbose)
						std::cout << "LloydsForGMM::run() - replaced empty covariance " << i+1 << " with unit sphere covariance." << std::endl;
				}

			}
		
		}
		else
		{
			for (std::size_t i=0; i<n; ++i)
			{
				Vector y = this->input->points.col(i) - this->desc.means.col(indices[i]);
				this->desc.covariances[indices[i]].noalias() += this->input->weights[i]*(y*y.transpose());
			}

			for (std::size_t i=0; i<k; ++i)
			{
				if (this->desc.weights[i] > 0)
				{
					this->desc.covariances[i] = this->desc.covariances[i]/this->desc.weights[i];
					if (!linalg::spd(this->desc.covariances[i]))
					{
						this->desc.covariances[i] = (this->desc.covariances[i]+covs[i])/2.0;
						if (!linalg::spd(this->desc.covariances[i]))
						{
							this->desc.covariances[i] = covs[i];
							if (this->verbose)
								std::cout << "LloydsForGMM::run() - replaced non-spd covariance " << i+1 << " with previous covariance." << std::endl;
						}
						else
							if (this->verbose)
								std::cout << "LloydsForGMM::run() - mixed non-spd covariance " << i+1 << " with previous covariance." << std::endl;
					}
				}
				else
				{
					this->desc.weights[i] = totalWeight/k;

					Vector sqDist = (this->desc.means.colwise()-this->desc.means.col(i)).colwise().squaredNorm();
					fp_type maxd = sqDist.maxCoeff();
					for (std::size_t j=0; j<k; ++j)
						if (this->desc.means.col(i)==this->desc.means.col(j))
							sqDist[j] = maxd;

					fp_type factor = (1.0/(2*d))*sqDist.minCoeff();
					if (!std::isfinite(factor) || factor <= 0) // i.e., factor*I_d would not be positive definite
					{
						this->desc.covariances[i].noalias() = (1.0/(2*d))*Matrix::Identity(d,d);
						if (this->verbose)
							std::cout << "LloydsForGMM::run() - replaced zero covariance " << i+1 << " by scaled identity matrix, factor = 1.0/(2*d)" << std::endl;
					}
					else
					{
						this->desc.covariances[i].noalias() = factor*Matrix::Identity(d,d);
						if (this->verbose)
							std::cout << "LloydsForGMM::run() - replaced zero covariance " << i+1 << " by scaled identity matrix, factor = " << sqDist.minCoeff() << "/(2*d)" << std::endl;
					}
				}

			}
		}

		// finally normalize the weights
		this->desc.weights /= this->desc.weights.sum();
	}

	this->change = true;

//	std::cout << "solution after " << numSteps << " step(s): " << this->desc << std::endl;

	end = clock();
	this->runtime += double(end-start)/CLOCKS_PER_SEC;
}

LloydsForGMM* LloydsForGMM::toLloydsForGMM(Algorithm* a)
{
	return dynamic_cast<LloydsForGMM*>(a);
}

