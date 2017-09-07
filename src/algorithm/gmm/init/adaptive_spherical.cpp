#include "adaptive_spherical.h"

#include "../../../base.h"


const std::string AdaptiveSphericalID::CLASSTAG = "AdaptiveSpherical";

AdaptiveSphericalID::AdaptiveSphericalID(uint32_t s) : seed(s)
{
	std::stringstream sstream;
	sstream << "AdaptiveSpherical" << "_i" << s;
	this->nametag = sstream.str();
}

Parameters AdaptiveSphericalID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	std::mt19937 gen(seed);
	return initializer::adaptiveSphericalGMM(input, k, gen);
}

void AdaptiveSphericalID::parseDescriptionString(std::vector< InitDescription* >& descriptions, std::string descriptionstring)
{
	std::istringstream iss(descriptionstring);
	std::vector<unsigned int> sList = readIntList(iss);
	for(idx_type s=0; s<sList.size(); ++s)
		descriptions.push_back(new AdaptiveSphericalID(sList[s]));
}

Parameters initializer::adaptiveSphericalGMM(commonutil::DataSet const& input, idx_type k, std::mt19937& gen)
{
	idx_type d = input.points.rows();
	idx_type n = input.points.cols();
	
	initializer::check(k, d, n);

	Matrix samples(d,k+1);	// two samples for the first gaussian, one for each of the others
	Vector sqDist(k);		// minimum distance to nearest sample (for each sample)

	Parameters desc;
	for (idx_type i=0; i<k+1; ++i)
	{
		if (i<2)
			// draw the first two points uniformly
			samples.col(i) = input.points.col(commonutil::randomIndex(input.weights, gen));
		else
		{
			// draw next point w.r.t. current mixture
			Vector densities = gmmutil::adaptiveDensities(input, desc, 1);
			samples.col(i) = input.points.col(commonutil::randomIndex(densities, gen));
		}

//		std::cout << "sampled points = " << samples.col(i).transpose() << std::endl;

		if (i>0)
		{
			desc.weights = Vector::Constant(i, 1.0/i);
			desc.means = samples.block(0,1,d,i);
			desc.covariances = std::vector<Matrix>(i, (1.0/(2*d))*Matrix::Identity(d,d));

			if (i==1)	// after uniformly drawing two samples create the first model with one component
			{
				sqDist[0] = (samples.col(0)-samples.col(1)).squaredNorm();
				if (sqDist[0]<=0)
					sqDist[0] = 1;
			}
			else if (i==2)
			{
				sqDist[0] = (samples.col(1)-samples.col(2)).squaredNorm();
				if (sqDist[0]<=0)
					sqDist[0] = 1;
				sqDist[1] = sqDist[0];
			}
			else if (i>2)
			{
				Vector mins = (samples.block(0,1,d,i-1).colwise()-samples.col(i)).colwise().squaredNorm();

				// set sqDist to the minimum non-zero distance to one of the other samples
				// this is due to the possibility of sampling the same point multiple times
				sqDist[i-1] = 0;
				for (idx_type j=0; j<i-1; ++j)
					if (mins[j]>0 && (sqDist[i-1]==0 || mins[j]<sqDist[i-1]))
						sqDist[i-1] = mins[j];
				if (sqDist[i-1]<=0)
					sqDist[i-1] = 1;	// prevent zero covariances in pathological cases

				for (idx_type j=0; j<i-1; ++j)
					if (mins[j]>0 && mins[j]<sqDist[j])
						sqDist[j] = mins[j];
			}

			for (idx_type j=0; j<i; ++j)
				desc.covariances[j] *= sqDist[j];
		}
	}

	return desc;
}