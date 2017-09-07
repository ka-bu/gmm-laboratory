#ifndef LLOYDSFORGMM_H
#define LLOYDSFORGMM_H

#include "../base/randomalgorithm.h"
#include "../description/algo_description.h"
#include "../../base/parameters.h"
#include "../../base/commonutil.h"

/**
* @brief Lloyd type algorithm for GMMs
*/
class LloydsForGMM : public RandomAlgorithm
{
public:
	/**
	* @brief Constructor initializing input and configuring the algorithm.
	* @remarks Remember to init the mixture model.
	*/
	LloydsForGMM(commonutil::DataSet const& ds, bool v, uint32_t seed, fp_type ratio = -1., bool spherical = false);
	
	LloydsForGMM(commonutil::DataSet const& ds, bool v, std::mt19937& gen, fp_type ratio = -1., bool spherical = false);

	virtual ~LloydsForGMM()
	{
	}

	virtual void init(Parameters const& gmmdesc);
	virtual void run(unsigned int numSteps = 1);

	/**
	* @brief does a dynamic cast of the given GMMalgorithm to LloydsForGMM
	* @return NULL if the Algorithm is not a LloydsForGMM instance
	*/
	static LloydsForGMM* toLloydsForGMM(Algorithm* a);

private:
	fp_type ratio;
	bool spherical;
};

/**
 *@brief: description 
 */
class LloydsForGMMAD : public AlgoDescription
{
public:
	static const std::string CLASSTAG;
	LloydsForGMMAD(unsigned int st, uint32_t sd, fp_type ratio, bool spherical);
	virtual Algorithm* create(const commonutil::DataSet& input) const;
	
	static void parseDescriptionString(std::vector<AlgoDescription*>& descriptions, const std::string descriptionstring, const bool first = false);
private:
	fp_type ratio;
	bool spherical;
	uint32_t seed;
	
};

#endif
