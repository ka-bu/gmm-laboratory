#include "emforgmm.h"

#include <algorithm>
#include <iostream>
#include <cmath>

#include "utils/gmmutil.h"
#include "../../base/linalgutil.h"
#include "../../settings/settings.h"



const std::string EMforGMMAD::CLASSTAG = "EM";

std::ostream& operator<<(std::ostream& os, const EMforGMMAD& emforgmmag)
{
	os << EMforGMMAD::CLASSTAG << "(" << emforgmmag.getNameTag() << ", " << emforgmmag.getSteps() << ", " << ((emforgmmag.getNext()!=NULL)?"hasNext":"noNext") << ")" << std::endl;
}


EMforGMMAD::EMforGMMAD(unsigned int st, uint32_t sd, bool spherical) : seed(sd), spherical(spherical)
{
	this->steps = st;
	std::stringstream sstream;
	sstream << "EM" << "_s" << st << (spherical?"_sph":"") << "_a" << sd;
	this->nametag = sstream.str();
}

Algorithm* EMforGMMAD::create(commonutil::DataSet const& input) const
{
	return new EMforGMM(input, commonSettings().verbose, this->seed, spherical);
}

void EMforGMMAD::parseDescriptionString(std::vector<AlgoDescription*>& descriptions, const std::string descriptionstring, const bool first)
{
	std::istringstream iss(descriptionstring);
	std::vector<unsigned int> stList = Description::readIntList(iss);
	bool spherical = (Description::readInt(iss) == 1);
	std::vector<unsigned int> sdList = Description::readIntList(iss);
	
	if (first)
	{
		for (std::size_t st=0; st<stList.size(); ++st)
			for (std::size_t sd=0; sd<sdList.size(); ++sd)
				descriptions.push_back(new EMforGMMAD(stList[st], sdList[sd], spherical));
	}
	else
	{
		if (stList.size()!=1)
			gmmlab_throw("multiple step values forbidden for next:...!");
		if (sdList.size()!=1)
			gmmlab_throw("multiple seed values forbidden for next:...!");
		descriptions.push_back(new EMforGMMAD(stList[0], sdList[0], spherical));
	}
	
	std::cout << "   created " << descriptions.size() << " instances of the " << EMforGMMAD::CLASSTAG << " algorithm" << std::endl;
}

// =================================================================================================================================


EMforGMM::EMforGMM(commonutil::DataSet const& ds, bool v, uint32_t s, bool spherical)
	: RandomAlgorithm(ds,v,s), spherical(spherical)
{
	this->initializeTotalAndMinWeight();
}

EMforGMM::EMforGMM(commonutil::DataSet const& ds, bool v, std::mt19937& gen, bool spherical)
	: RandomAlgorithm(ds,v,0), spherical(spherical)
{
	this->gen = gen;
	this->initializeTotalAndMinWeight();
}


void EMforGMM::initializeTotalAndMinWeight()
{
	bool indices_given = !this->indices_of_input_points_to_be_used.empty();
	
	// total weight of all points
	idx_type n = this->input->n();
	this->totalWeight = 0;
	if(indices_given)
		for(idx_type i=0; i<n; ++i)
			totalWeight += this->input->weights(this->indices_of_input_points_to_be_used.at(i));
	else
		totalWeight = this->input->weights.sum();
	if (totalWeight<=0)
		gmmlab_throw("Total input weight less or equal to zero.");

	// minimum weight of a point
	this->minWeight = totalWeight;
	if(indices_given)
	{
		idx_type index;
		for (idx_type i=0; i<n; ++i)
		{
			index = this->indices_of_input_points_to_be_used.at(i);
			if (this->input->weights[index]>0 && this->input->weights[index]<minWeight)
				minWeight = this->input->weights[index];
		}
	}
	else
	{
		for (idx_type i=0; i<n; ++i)
			if (this->input->weights[i]>0 && this->input->weights[i]<minWeight)
				minWeight = this->input->weights[i];
	}
}


void EMforGMM::setIndicesOfInputPointsToBeUsed(std::vector< idx_type > indices)
{
	this->indices_of_input_points_to_be_used = indices;
	this->initializeTotalAndMinWeight();
}


void EMforGMM::run(unsigned int numSteps)
{
	clock_t start, end;
	start = clock();

	idx_type k = this->desc.weights.size();
	bool indices_given = !this->indices_of_input_points_to_be_used.empty();
	idx_type n = indices_given ? this->indices_of_input_points_to_be_used.size() : this->input->points.cols();
	idx_type d = this->input->points.rows();
	
	if (n==0 || d==0)
		gmmlab_throw("Input is empty.");
	assert(this->desc.means.cols()==k && this->desc.covariances.size()==k);
	if (k==0)
		gmmlab_throw("No or empty initial solution given.");
	
	// run the em algorithm
	Matrix pmatrix(k,n);
	for (unsigned int step=0; step<numSteps; ++step)
	{
		// update the pmatrix
		
		// ... first compute the densities N(x|mu_i,Sigma_i)
		for (idx_type i=0; i<k; ++i)
			pmatrix.row(i) = this->desc.weights[i]*gmmutil::gaussianDensity(this->input->points,
				this->desc.means.col(i), this->desc.covariances[i], indices_of_input_points_to_be_used);
			
		// .. then normalize them to obtain the posterior probability p(i|x_j,theta)
		for (idx_type j=0; j<n; ++j)
		{
			fp_type sum = pmatrix.col(j).sum();

			if (sum<=0)
			{
				errorHandling("EMforGMM::run() - Lonely point:");
				if (this->verbose)
					std::cout << "EMforGMM::run() - no gaussian is responsible for point " << j << std::endl;
				pmatrix.col(j) = this->desc.weights;
			}
			else
				for (idx_type i=0; i<k; ++i)
					pmatrix(i,j) /= sum;
				// NOTE: pmatrix.col(j) /= sum;	// ist wie *= 1/sum und führt bei kleinen Einträgen zu fehlern
		}

		// ... multiply this probabilities with the weights, obtaining w_j*p(i|x_j,theta)
		if(indices_given)
			for(idx_type j=0; j<n; ++j)
				pmatrix.col(j) *= this->input->weights[indices_of_input_points_to_be_used.at(j)];
		else
			for (idx_type j=0; j<n; ++j)
				pmatrix.col(j) *= this->input->weights[j];

		// compute updates
		for (idx_type i=0; i<k; ++i)
		{
			// resp = normalization term for the update of the i-th component
			fp_type resp = pmatrix.row(i).sum();
			
			// only update mean and covariance, if the gaussian has non-negligible responsibility
			if (resp>minWeight)
			{
				// update weight
				this->desc.weights[i] = resp/totalWeight; // may become zero, no problem

				// update mean
				Vector p = Vector::Zero(d);
				if(indices_given)
					for (idx_type j=0; j<n; ++j)
						p.noalias() += pmatrix(i,j)*this->input->points.col(this->indices_of_input_points_to_be_used.at(j));
				else
					for (idx_type j=0; j<n; ++j)
						p.noalias() += pmatrix(i,j)*this->input->points.col(j);
				this->desc.means.col(i) = p/resp;
//std::cout << "EMforGMM::run() - EM mean " << i << " is " << this->desc.means.col(i) << std::endl;
				
				// update covariance
				Matrix m = Matrix::Zero(d,d);
				if(!this->spherical)
				{
					for (idx_type j=0; j<n; ++j)
					{
						Vector y;
						if(indices_given)
							y = this->input->points.col(this->indices_of_input_points_to_be_used.at(j)) - this->desc.means.col(i);
						else
							y = this->input->points.col(j) - this->desc.means.col(i); //ALS COLWISE NACH OBEN ZIEHEN?
						
						m.noalias() += pmatrix(i,j)*(y*y.transpose());						
					}
					m = m/resp;
				}
				else
				{
					fp_type sph = 0;
					
					for (std::size_t j=0; j<n; ++j)
					{
						Vector y =  this->input->points.col(j) - desc.means.col(i);
						sph += pmatrix(i,j)*input->weights[j]*y.squaredNorm();
					}
					m.noalias() = Matrix::Identity(d,d)*(sph/(resp*d));
				}
				
				
				
				// error handling: if covariance is not spd, then...
				if (linalg::spd(m))
					this->desc.covariances[i] = m;
				else
				{
					errorHandling("EMforGMM::run() - Covariance matrix is not SPD:");
					// ... try to mix it with the old covariance or just don't update the covariance
					m = (this->desc.covariances[i]+m)/2.0;
					if (linalg::spd(m))
					{
						this->desc.covariances[i] = m;
						if (this->verbose)
							std::cout << "EMforGMM::run() - matrix not spd, mixing with old covariance of gaussian " << i+1 << std::endl;
					}
					else
						if (this->verbose)
							std::cout << "EMforGMM::run() - matrix not spd, keeping old covariance " << i+1 << std::endl;
				}
			}
			else
			{
				// if the gaussian has negligible responsibility ...
				
				errorHandling("EMforGMM::run() - Some Cluster has negligible weight:");
				
				// ... set the weight to 1/k
				this->desc.weights[i] = 1.0/k;

				// ... choose a new mean randomly
				Vector random_mean;
				if(indices_given)
				{
					std::uniform_int_distribution<unsigned int> uid(0,n-1);
					random_mean = this->input->points.col(this->indices_of_input_points_to_be_used.at(uid(this->gen)));
				}
				else
					random_mean = this->input->points.col(commonutil::randomIndex(this->input->weights, this->gen));
				this->desc.means.col(i) = random_mean; 
				
				// notify
 				if (this->verbose)
					std::cout << "EMforGMtestlab::run() - Cluster " << i+1 << " hast negligible weight. "
						  << "                  Replace mean of the empty cluster by a random input point and re-estimate covariance by heuristic." << std::endl;

				// ... initialize the covariance matrix
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
						std::cout << "EMforGMM::run() - replaced covariance " << i+1 << " by scaled identity matrix, factor = 1.0/(2*d)" << std::endl;
				}
				else
				{
					this->desc.covariances[i].noalias() = factor*Matrix::Identity(d,d);
					if (this->verbose)
						std::cout << "EMforGMM::run() - replaced covariance " << i+1 << " by scaled identity matrix, factor = " << sqDist.minCoeff() << "/(2*d)" << std::endl;
				}
				
			}
		}

		// finally normalize the weights (maybe weights were replaced and the sum is not 1)
		this->desc.weights /= this->desc.weights.sum();
	}

	this->change = true;

	end = clock();
	this->runtime += double(end-start)/CLOCKS_PER_SEC;
}

EMforGMM* EMforGMM::toEMforGMM(Algorithm* a)
{
	return dynamic_cast<EMforGMM*>(a);
}

