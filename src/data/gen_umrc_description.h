#ifndef GEN_UMRC_DESCRIPTION_H
#define GEN_UMRC_DESCRIPTION_H

#include "data_description.h"

#include <vector>
#include <boost/concept_check.hpp>

/**
 * 
 */
class GenUMRCDescription: public DataDescription
{
public:
	GenUMRCDescription(unsigned int n, unsigned int d, unsigned int k, fp_type weightFloatExp, 
										fp_type separation,
										fp_type sqrtEWFloatExp, fp_type minMinSqrtEW, fp_type maxMinSqrtEW, fp_type minSqrtEWProp, fp_type maxSqrtEWProp, 
										std::string measurementError,   
										uint32_t seed);

	virtual ~GenUMRCDescription()
	{
	}
	
	Parameters retrieve(commonutil::DataSet& input);

private:
	unsigned int size, dim, comp;
	fp_type weightFloatExp;
	fp_type separation;
	fp_type sqrtEWFloatExp, minMinSqrtEW, maxMinSqrtEW, minSqrtEWProp, maxSqrtEWProp;
	std::string measurementError;
	
	uint32_t genSeed;
};

#endif

