#include "adaptive_means.h"

#include "../../../base.h"
#include "../utils/kmeansutil.h"

const std::string AdaptiveMeansID::CLASSTAG = "AdaptiveMeans";

AdaptiveMeansID::AdaptiveMeansID(uint32_t s) : seed(s)
{
	std::stringstream sstream;
	sstream << "AdaptiveMeans" << "_i" << s;
	this->nametag = sstream.str();
}

void AdaptiveMeansID::parseDescriptionString(std::vector< InitDescription* >& descriptions, std::string descriptionstring)
{
	std::istringstream iss(descriptionstring);
	std::vector<unsigned int> sList = Description::readIntList(iss);
	for(idx_type s=0; s<sList.size(); ++s)
		descriptions.push_back(new AdaptiveMeansID(sList[s]));
}


Parameters AdaptiveMeansID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	//std::cout << "adaptive means... k=" << k << ", input = " << input.points.rows() << "x" << input.points.cols() << "." << std::endl;
	
	std::mt19937 gen(seed);
	return  initializer::wrappedAdaptiveMeans(input, k, gen);
	
	//std::cout << "initial = " << initial << std::endl;	
	//return initial;
}

Matrix initializer::adaptiveMeans(commonutil::DataSet const& input, idx_type k, std::mt19937& gen)
{
	idx_type n = input.points.cols();
	idx_type d = input.points.rows();

	initializer::check(k, d, n);
	
	Matrix means(d,k);
	Vector sqNorms;
	for (idx_type i=0; i<k; ++i)
		if (i==0)
			means.col(i) = input.points.col(commonutil::randomIndex(input.weights, gen));
		else
		{
			if (i==1)
				sqNorms = (input.points.colwise()-means.col(0)).colwise().squaredNorm();
			else
			{
				for (idx_type j=0; j<n; ++j)
				{
					fp_type sqn = (input.points.col(j)-means.col(i-1)).squaredNorm();
					if (sqn<sqNorms[j])
						sqNorms[j]=sqn;
				}
			}

			means.col(i) = input.points.col(commonutil::randomIndex(sqNorms, gen));
		}

	return means;
}


Parameters initializer::wrappedAdaptiveMeans(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, bool computeWeightAndCovar)
{
	return kmeansutil::wrapMeans(input, adaptiveMeans(input,k,gen), computeWeightAndCovar);
}
