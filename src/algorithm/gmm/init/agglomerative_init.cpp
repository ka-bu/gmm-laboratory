#include "agglomerative_init.h"

#include "../../../base.h"
#include "../../../base/similarities.h"
#include <boost/algorithm/string.hpp>

const std::string AgglomerativeInitID::CLASSTAG = "AgglomerativeInit";

AgglomerativeInitID::AgglomerativeInitID(fp_type sampleSizeFactor, uint32_t s) : seed(s), sampleSizeFactor(sampleSizeFactor)
{
	std::stringstream sstream;
	sstream << "AgglomerativeInit_" << "_sample" << sampleSizeFactor << "_i" << s;
	this->nametag = sstream.str();
}

Parameters AgglomerativeInitID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	std::mt19937 gen(seed);
	fp_type sampleSizeFactor = this->sampleSizeFactor;
	if(this->sampleSizeFactor == 0)
		sampleSizeFactor = 1./k;
	return initializer::sampleAgglomerativeMeansToGMM(input, k, sampleSizeFactor, 0, gen);
}

void AgglomerativeInitID::parseDescriptionString(std::vector< InitDescription* >& descriptions, std::string descriptionstring, bool first)
{
	std::istringstream iss(descriptionstring);
	std::vector<unsigned int> sList = readIntList(iss);
	std::vector<fp_type> initSampleSizeFactorList = readFpTypeList(iss);
	for(idx_type s=0; s<sList.size(); ++s)
		for(idx_type sf = 0; sf<initSampleSizeFactorList.size(); ++sf)
			descriptions.push_back(new AgglomerativeInitID(initSampleSizeFactorList[sf],sList[s]));
}


Parameters initializer::sampleAgglomerativeMeansToGMM(const commonutil::DataSet& input, const idx_type k, const fp_type sampleSizeFactor, const bool precompute, std::mt19937& gen)
{
	idx_type d = input.points.rows();
	idx_type n = input.points.cols();
	
	Matrix means;
	
	idx_type sampleSize = (int) round(n*sampleSizeFactor);
	if(sampleSizeFactor < n)
	{
		commonutil::DataSet sample;
		commonutil::get_uniform_subset(input, sample, sampleSize,gen);
		means = initializer::agglomerativeMeans(sample, k, precompute);
	}
	else
		means = initializer::agglomerativeMeans(input, k, precompute);
	
	return gmmutil::meansToGMM(input, means, false, false);
}

// TODO: only works with unweighted points
Matrix initializer::agglomerativeMeans(commonutil::DataSet const& input, const idx_type k, const bool precompute)
{
	Parameters desc;

	if (k==0)
		gmmlab_throw("Empty mixture requested.");

	idx_type n = input.points.cols();
	idx_type d = input.points.rows();

	if (n==0 || d==0)
		gmmlab_throw("Input is empty.");


	// initialize partition with trivial n-clustering
	std::vector<std::vector<idx_type>> partition;
	std::vector<Vector> partition_sums;
	for(idx_type i=0; i<n; ++i)
	{
	  std::vector<idx_type> cluster;
	  cluster.push_back(i);
	  partition.push_back(cluster);
	  partition_sums.push_back(input.points.col(i));
	}

	// distance measure
	AverageLinkage averageLinkageDis;
	

	// precompute dissimilarities
	Matrix dis;
	if (precompute)
	{
		dis = Matrix::Zero(n,n);
		for (idx_type i=0; i<n; ++i)	
			for (idx_type j=0; j<=i; ++j)
				dis(i,j) = averageLinkageDis(partition_sums.at(i), partition.at(i).size(), partition_sums.at(j), partition.at(j).size());
	}

	for (idx_type r=n; r>k; --r)
	{				
		idx_type first = -1;
		idx_type second = -1;
		fp_type min = FP_INFINITE;
		for (idx_type i=0; i<r; ++i)
			for (idx_type j=0; j<i; ++j) // j < i
			{
				//std::cout << "i=" << i << ", j=" << j << std::endl;
				if (i!=j)
				{
					fp_type d;
					if (precompute)
						d = dis(i,j);
					else
						d = averageLinkageDis(partition_sums.at(i), partition.at(i).size(), partition_sums.at(j), partition.at(j).size());

					if ((i==1 && j==0) || d < min)
					{
						min = d;
						first = j;
						second = i;
					}
				}
			}
		
		// merge clusters (note: first < second)
		// 1. update sufficient statistics 
		partition_sums.at(first) = partition_sums.at(first)+partition_sums.at(second);
		// 2. move points from second cluster to the first
		while(!partition.at(second).empty())
		{
		  idx_type point = partition.at(second).back();
		  partition.at(second).pop_back();
		  partition.at(first).push_back(point);
		} 
		// overwrite second cluster with the last stored cluster and remove the last cluster  (note: first < second)
		partition.at(second) = partition.back();
		partition_sums.at(second) = partition_sums.back();
		partition.pop_back();
		partition_sums.pop_back();
		
		if (precompute)
		{
			// cluster with index second has been overwritten by the last cluster (with index r-1)
			for(idx_type i=0; i<r-1; i++)
			{
				if(i!=second)
				{
					fp_type d = dis(r-1,i);
					if(i < second)
						dis(second, i) = d;
					else if(second < i)
						dis(i,second) =  d;
				}
			}
			
			// compute distances wrt the newly formed cluster which is stored at index first
			for (idx_type i=0; i<r-1; i++)
			{
				if(i!=first){
					fp_type d = averageLinkageDis(partition_sums.at(first), partition.at(first).size(), partition_sums.at(i), partition.at(i).size());
					if(i < first)
						dis(first, i) = d;
					else if(first < i)
						dis(i,first) =  d;
				}
			}
		}
	}
	
	Matrix means(d,k);
	for(idx_type i=0; i<k; ++i)
		means.col(i) = partition_sums.at(i) / partition.at(i).size();
	
	
	return means;
}
