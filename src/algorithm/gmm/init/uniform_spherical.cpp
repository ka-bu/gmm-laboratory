#include "uniform_spherical.h"

#include "../../../base.h"

const std::string UniformSphericalID::CLASSTAG = "UniformSpherical";

UniformSphericalID::UniformSphericalID(uint32_t s) : seed(s)
{
	std::stringstream sstream;
	sstream << "UniformSpherical" << "_i" << s;
	this->nametag = sstream.str();
}

void UniformSphericalID::parseDescriptionString(std::vector< InitDescription* >& descriptions, std::string descriptionstring)
{
	std::istringstream iss(descriptionstring);
	std::vector<unsigned int> sList = Description::readIntList(iss);
	for(idx_type s=0; s<sList.size(); ++s)
		descriptions.push_back(new UniformSphericalID(sList[s]));
}


Parameters UniformSphericalID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	std::mt19937 gen(seed);
	return initializer::uniformSphericalGMM(input, k, gen);
}


Parameters initializer::uniformSphericalGMM(commonutil::DataSet const& input, idx_type k, std::mt19937& gen)
{
	idx_type n = input.points.cols();
	idx_type d = input.points.rows();
	
	initializer::check(k, d, n);

	Parameters desc;
	desc.weights = Vector::Constant(k, 1.0/k);
	desc.covariances = std::vector<Matrix>(k, (1.0/(2*d))*Matrix::Identity(d, d));

	desc.means = Matrix(d,k);
	for (size_t i=0; i<k; ++i)
		desc.means.col(i) = input.points.col(commonutil::randomIndex(input.weights, gen));

	for (size_t i=0; i<k; ++i)
	{
		Vector sq = (desc.means.colwise()-desc.means.col(i)).colwise().squaredNorm();

		fp_type min = 0;
		for (size_t j=0; j<k; ++j)
			if (j!=i && sq[j]>0 && (min==0 || sq[j]<min))
				min = sq[j];
		if (min<=0)
			min = 1;

		desc.covariances[i] *= min;
	}

	return desc;
}