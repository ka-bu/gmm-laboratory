#include "gonzalez_forgmm.h"

#include "../../../base.h"

const std::string GonzalezForGMMID::CLASSTAG = "GonzalezForGMM";

GonzalezForGMMID::GonzalezForGMMID(uint32_t s, bool use2GMM, bool spherical, fp_type sampleSizeFactor) : seed(s), use2GMM(use2GMM), sampleSizeFactor(sampleSizeFactor), spherical(spherical)
{
	std::stringstream sstream;
	sstream << "GonzalezForGMM";
	if(spherical)
		sstream << "_spherical";
	sstream << "_i" << s << "_2gmm" << use2GMM <<"_sample" << sampleSizeFactor ;
	this->nametag = sstream.str();
}

Parameters GonzalezForGMMID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	std::mt19937 gen(seed);
	//(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, bool use2GMM, fp_type sampleFactor)
	return initializer::gonzalezForGMM(input, k, gen, use2GMM, spherical, sampleSizeFactor);
}

void GonzalezForGMMID::parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring)
{
	std::istringstream iss(descriptionstring);
	std::vector<unsigned int> initUse2GMMList = readIntList(iss);
	bool initSphericalFlag = (readInt(iss)==1);
	std::vector<fp_type> initSampleSizeFactorList = readFpTypeList(iss);
	std::vector<unsigned int> sList = Description::readIntList(iss);
	
	for (idx_type s=0; s<sList.size(); ++s)
		for(idx_type u = 0; u<initUse2GMMList.size(); ++u)
			for(idx_type sf = 0; sf<initSampleSizeFactorList.size(); ++sf)
				descriptions.push_back(new GonzalezForGMMID(sList[s],(initUse2GMMList[u]>0), initSphericalFlag, initSampleSizeFactorList[sf]));
}

Parameters initializer::gonzalezForGMM(const commonutil::DataSet& input, idx_type k, std::mt19937& gen, bool use2GMM, bool spherical, fp_type sampleFactor)
{
	
// 	std::cout << "sampleFactor = " << sampleFactor << std::endl;	
// 	std::cout << "use2GMM = " << use2GMM << std::endl;	
	
	idx_type d = input.points.rows();
	idx_type n = input.points.cols();
	
	initializer::check(k, d, n);
		
	idx_type sampleSize = (int) round(n*sampleFactor);
	commonutil::DataSet samples;
	if(sampleSize < n)
		commonutil::get_uniform_subset(input, samples, sampleSize, gen);

	Parameters desc;
	Vector sample;
	
	fp_type trace;
	if(!use2GMM)
	{
		Parameters gmmdesc = gmmutil::optimalGaussian(input);
		Matrix covar = gmmdesc.covariances.at(0);
		trace = covar.trace()/(10.*d*k);
		
		desc.means.resize(0,0);
		desc.weights.resize(0);
	}
		
	Matrix randMatrix(d,d);
	Matrix randOrthonormalMatrix(d,d);
	Vector eigenvalues(d);
	std::uniform_real_distribution<> urd(0,1);
	
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
		
		if(use2GMM)
		{
			Matrix tmpMeans = Matrix::Zero(d,i+1);
			if(i>0)
				tmpMeans.block(0,0,d,i) = desc.means;
			tmpMeans.col(i) = sample;
			desc = gmmutil::meansToGMM(input, tmpMeans, spherical, false);
		}
		else
		{
			desc.means.conservativeResize(d,i+1);
			desc.means.col(i) = sample;
			
			desc.weights.conservativeResize(i+1);
			
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
	}
	
	if(!use2GMM)
	{
		// estimate weights
		desc.weights = Vector::Zero(k);
		std::vector<idx_type> partition = gmmutil::gaussPartition(input.points, desc);
		for(idx_type n=0; n<input.points.cols(); ++n)
			++desc.weights(partition.at(n));		
		desc.weights /= input.points.cols();
	}

	return desc;
}
