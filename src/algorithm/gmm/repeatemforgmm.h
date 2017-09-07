#ifndef REPEATEMFORGMM_H
#define REPEATEMFORGMM_H

#include <iostream>

#include "emforgmm.h"
#include "../base/randomalgorithm.h"
#include "../description/algo_description.h"
#include "../../base/commonutil.h"
#include "utils/gmmutil.h"

/**
 * 
 * 
 */
class RepeatEMforGMM : public RandomAlgorithm
{
public:
	
	RepeatEMforGMM(commonutil::DataSet const& input, bool verbose, uint32_t seed, unsigned int substeps, bool spherical = false);
	RepeatEMforGMM(commonutil::DataSet const& ds, bool v, std::mt19937& gen, unsigned int subst, bool spherical);
	virtual ~RepeatEMforGMM()
	{
		delete emforgmm;
	}
		
	virtual void run(unsigned int numSteps = 1);
	
	void setIndicesOfInputPointsToBeUsed(std::vector<idx_type> indices);
	

private:
	bool spherical = false;
	unsigned int substeps;
	EMforGMM* emforgmm;
	
	Parameters best_desc;
	fp_type best_cost;
	
};

/**
 *@brief: description 
 */
class RepeatEMforGMMAD : public AlgoDescription
{
public:
	static const std::string CLASSTAG;
	RepeatEMforGMMAD(unsigned int st, uint32_t sd, bool spherical, unsigned int substeps);
	virtual Algorithm* create(const commonutil::DataSet& input) const;
	
	static void parseDescriptionString(std::vector<AlgoDescription*>& descriptions, const std::string descriptionstring, const bool first = false);
private:
	uint32_t seed;
	unsigned int substeps;
	bool spherical;
	
};


#endif
