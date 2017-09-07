#include "gonzalez_means.h"

#include "../utils/kmeansutil.h"

const std::string GonzalezMeansID::CLASSTAG = "GonzalezMeans";

GonzalezMeansID::GonzalezMeansID(uint32_t s) : seed(s)
{
	std::stringstream sstream;
	sstream << "GonzalezMeans" << "_i" << s;
	this->nametag = sstream.str();
}

void GonzalezMeansID::parseDescriptionString(std::vector< InitDescription* >& descriptions, std::string descriptionstring)
{
	std::istringstream iss(descriptionstring);
	std::vector<unsigned int> sList = Description::readIntList(iss);
	for(idx_type s=0; s<sList.size(); ++s)
		descriptions.push_back(new GonzalezMeansID(sList[s]));
}


Parameters GonzalezMeansID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	//std::cout << "adaptive means... k=" << k << ", input = " << input.points.rows() << "x" << input.points.cols() << "." << std::endl;
	
	std::mt19937 gen(seed);
	return  initializer::wrappedGonzalez(input, k, gen, true);
	
	//std::cout << "initial = " << initial << std::endl;	
	//return initial;
}

Matrix initializer::gonzalez(commonutil::DataSet const& input, idx_type k,  std::mt19937& gen)
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
			
			idx_type index;
			sqNorms.maxCoeff(&index);
			means.col(i) = input.points.col(index);
		}
		
		return means;
}



Parameters initializer::wrappedGonzalez(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, bool computeWeightAndCovar)
{	
	return kmeansutil::wrapMeans(input, initializer::gonzalez(input, k, gen), computeWeightAndCovar);
}