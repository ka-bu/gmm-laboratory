#include "gonzalez_kwedlo.h"

#include "../../../base.h"

const std::string GonzalezKwedloID::CLASSTAG = "GonzalezKwedlo";

GonzalezKwedloID::GonzalezKwedloID(uint32_t s, fp_type sampleSizeFactor) : seed(s), sampleSizeFactor(sampleSizeFactor)
{
	std::stringstream sstream;
	sstream << "GonzalezKwedlo" << "_i" << s <<"_sample" << sampleSizeFactor ;
	this->nametag = sstream.str();
}

Parameters GonzalezKwedloID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	std::mt19937 gen(seed);
	//(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, fp_type sampleFactor)
	return initializer::gonzalezKwedlo(input, k, gen, sampleSizeFactor);
}

void GonzalezKwedloID::parseDescriptionString(std::vector< InitDescription* >& descriptions, std::string descriptionstring)
{
	std::istringstream iss(descriptionstring);
	std::vector<fp_type> initSampleSizeFactorList = readFpTypeList(iss);
	std::vector<unsigned int> sList = Description::readIntList(iss);
	for(idx_type s=0; s<sList.size(); ++s)
		for(idx_type sf = 0; sf<initSampleSizeFactorList.size(); ++sf)
			descriptions.push_back(new GonzalezKwedloID(sList[s],initSampleSizeFactorList[sf]));
}


Parameters initializer::gonzalezKwedlo(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, fp_type sampleFactor)
{  		
	idx_type d = input.points.rows();
	idx_type n = input.points.cols();
	
	initializer::check(k, d, n);
	
	idx_type sampleSize = (int) round(n*sampleFactor);
	commonutil::DataSet samples;
	if(sampleSize < n)
		commonutil::get_uniform_subset(input, samples, sampleSize, gen);

	Vector sample;
	Parameters gmmdesc = gmmutil::optimalGaussian(input);
	Matrix covar = gmmdesc.covariances.at(0);
	fp_type trace = covar.trace()/(10.*d*k);	
	
	Parameters desc;
	desc.means.resize(0,0);
	desc.weights.resize(0);

		
	std::uniform_real_distribution<> urd(0,1);
	std::uniform_int_distribution<> uid(0,input.points.cols()-1);
	
	Matrix randMatrix(d,d);
	Matrix randOrthonormalMatrix(d,d);
	Vector eigenvalues(d);
	
	for (idx_type i=0; i<k; ++i)
	{
		
		if (i==0)
			// draw first point uniformly
			sample = input.points.col(commonutil::randomIndex(input.weights, gen)); 
		else
		{
			// draw next point w.r.t. current mixture
			idx_type index;
			if(sampleSize < n)
			{
				Vector densities = gmmutil::adaptiveDensities(samples, desc, 1.);	
				densities.maxCoeff(&index);
				sample = samples.points.col(index); 
			}
			else
			{
				Vector densities = gmmutil::adaptiveDensities(input, desc, 1.);	
				densities.maxCoeff(&index);
				sample = input.points.col(index); 
			}
			
		}
		
		desc.means.conservativeResize(d,i+1);
		desc.means.col(i) = sample;
		
		desc.weights.conservativeResize(i+1);
		desc.weights(i) = urd(gen);
			
		// random covariance
		for (idx_type i=0; i<d; ++i)
			for (idx_type j=0; j<d; ++j)
				randMatrix(i,j)=urd(gen);
		Eigen::HouseholderQR<Matrix> qr(randMatrix);
		randOrthonormalMatrix = qr.householderQ();
		commonutil::fill(eigenvalues,gen, 1.,10.);
		fp_type tmp = trace/eigenvalues.sum();
		eigenvalues = tmp * eigenvalues;
		randMatrix = randOrthonormalMatrix.transpose()*eigenvalues.asDiagonal()*randOrthonormalMatrix;
			
		desc.covariances.push_back(randMatrix);
		
	}
	

	desc.weights /= desc.weights.sum();

	return desc;
}
