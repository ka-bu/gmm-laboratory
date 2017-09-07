#ifndef EMFORGMM_H
#define EMFORGMM_H

#include <iostream>

#include "../base/randomalgorithm.h"
#include "../description/algo_description.h"
#include "../../base/commonutil.h"
#include "utils/gmmutil.h"

/**
 * @brief EM algorithm for GMMs
 */
class EMforGMM : public RandomAlgorithm
{
public:
	
	EMforGMM(commonutil::DataSet const& input, bool verbose, uint32_t seed, bool spherical = false);
	EMforGMM(commonutil::DataSet const& input, bool verbose, std::mt19937& gen, bool spherical = false);
	virtual ~EMforGMM()
	{
	};
		
	void run(unsigned int numSteps = 1) override;

	void setIndicesOfInputPointsToBeUsed(std::vector<idx_type> indices);
	
	/**
	* @brief does a dynamic cast of the given GMMalgorithm to EMforGMM
	* @return NULL if the Algorithm is not a EMforGMM instance
	*/
	static EMforGMM* toEMforGMM(Algorithm* a);

protected:
	void initializeTotalAndMinWeight();
	
	std::vector<idx_type> indices_of_input_points_to_be_used;
	bool spherical = false;
	fp_type minWeight, totalWeight;
	
};

/**
 *@brief: description 
 */
class EMforGMMAD : public AlgoDescription
{
public:
	static const std::string CLASSTAG;
	EMforGMMAD(unsigned int st, uint32_t sd, bool spherical);
	Algorithm* create(const commonutil::DataSet& input) const override;
	
	static void parseDescriptionString(std::vector<AlgoDescription*>& descriptions, const std::string descriptionstring, const bool first = false);
private:
	uint32_t seed;
	bool spherical;
	
};

std::ostream& operator<<(std::ostream&os, const EMforGMMAD& emforgmmad
);

#endif
