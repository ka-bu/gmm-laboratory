#ifndef LLOYDS_H
#define LLOYDS_H

#include "../base/randomalgorithm.h"
#include "../description/algo_description.h"

#include "../../base/parameters.h"
#include "../../base/commonutil.h"

/**
 * (Randomized) Lloyds algorithm.
 * 
 * @param ratio determines the probability that the assignment to the clusters is done 
 * randomly. I.e., if ratio is set to 0, then the assignment is deterministic.
 * If ratio is set to 1, then the assignement is always done randomly
 * @brief Lloyds algorithm
 */
class Lloyds : public RandomAlgorithm
{
public:
	/**
	* @brief Constructor initializing input and configuring the algorithm.
	* @remarks Remember to init the mixture model.
	* 
	* @param ratio determines the probability that the assignment to the clusters is done 
	* randomly. I.e., if ratio is set to <= 0, then the assignment is deterministic.
	* If ratio is set to 1, then the assignement is always done randomly
	*/
	Lloyds(commonutil::DataSet const& input, bool verbose, uint32_t seed, fp_type ratio = -1);

	Lloyds(commonutil::DataSet const& input, bool verbose, std::mt19937& gen,  fp_type ratio = -1);
	
	virtual ~Lloyds()
	{
	}
	
	void init(Matrix const& means);
	
	virtual void run(unsigned int numSteps = 1);

	/**
	* @brief does a dynamic cast of the given GMMalgorithm to Lloyds
	* @return NULL if the Algorithm is not a Lloyds instance
	*/
	static Lloyds* toLloyds(Algorithm* a);

private:
	fp_type ratio;
};

/**
 *@brief: description 
 */
class LloydsAD : public AlgoDescription
{
public:
	static const std::string CLASSTAG;
	LloydsAD(unsigned int steps,  fp_type ratio = 0, uint32_t seed = 1);
	virtual Algorithm* create(const commonutil::DataSet& input) const;
	
	static void parseDescriptionString(std::vector<AlgoDescription*>& descriptions, const std::string descriptionstring, const bool first = false);
private:
	uint32_t seed;
	fp_type ratio;
	
};

#endif
