#include "uniform_means.h"

#include "../../../base.h"
#include "../utils/kmeansutil.h"

const std::string UniformMeansID::CLASSTAG = "UniformMeans";

UniformMeansID::UniformMeansID(uint32_t s) : seed(s)
{
	std::stringstream sstream;
	sstream << "UniformMeans" << "_i" << s;
	this->nametag = sstream.str();
}

void UniformMeansID::parseDescriptionString(std::vector< InitDescription* >& descriptions, std::string descriptionstring)
{
	std::istringstream iss(descriptionstring);
	std::vector<unsigned int> sList = Description::readIntList(iss);
	for(idx_type s=0; s<sList.size(); ++s)
		descriptions.push_back(new UniformMeansID(sList[s]));
}


Parameters UniformMeansID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	std::mt19937 gen(seed);
	return initializer::wrappedUniformMeans(input, k, gen);
}


Matrix initializer::uniformMeans(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, std::vector<idx_type> const& indices_of_points_to_be_used)
{
	idx_type n = indices_of_points_to_be_used.empty() ? input.points.cols() : indices_of_points_to_be_used.size();
	idx_type d = input.points.rows();

	initializer::check(k, d, n);
	
	Matrix means(d,k);
	if(indices_of_points_to_be_used.empty())
		for (idx_type i=0; i<k; ++i)
			means.col(i) = input.points.col(commonutil::randomIndex(input.weights, gen));	
	else
	{
		std::uniform_int_distribution<> uid(0,n-1);
		for (idx_type i=0; i<k; ++i)
			means.col(i) = input.points.col(indices_of_points_to_be_used.at(uid(gen)));	
	}
	
	return means;
}



Parameters initializer::wrappedUniformMeans(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, bool computeWeightAndCovar)
{
	return kmeansutil::wrapMeans(input, uniformMeans(input,k,gen), computeWeightAndCovar);
}
