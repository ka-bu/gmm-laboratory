#ifndef DATAUTIL_H
#define DATAUTIL_H

#include "../../base/parameters.h"
#include "../../base/commonutil.h"

#include <random>

namespace datautil
{
	/**
	 * drawFromSG draws from a given single Gaussian.
	 * drawFromSMPE draws from a given single MPE.
	 * drawFromMPEMM draws from a given mixture of MPEs.
	 */
	Matrix drawFromSG(Vector const& mean, Matrix const& cov, std::mt19937& gen, idx_type n, bool isCholesky = false);
	Matrix drawFromSMPE(Vector const& mean, Matrix const& cov, fp_type const& kappa, std::mt19937& gen, idx_type n, bool isCholesky = false);
	Matrix drawFromMPEMM(Parameters const& desc, std::mt19937& gen, idx_type n, bool isCholesky = false);
	
	struct DataStats
	{
		Vector spreads;
		fp_type maxSpread = 0;
		fp_type totalWeight = 0;
		fp_type minWeight = 0;
	};
	
	void stats(commonutil::DataSet const& ds, DataStats& dataStats);
	
	
	

}

#endif
