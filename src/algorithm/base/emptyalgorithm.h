#ifndef EMPTY_ALGO_H
#define EMPTY_ALGO_H

#include <iostream>

#include "randomalgorithm.h"
#include "../description/algo_description.h"

/**
 * @brief does nothing
 */
class EmptyAlgorithm : public Algorithm
{
public:
	
	EmptyAlgorithm();
	virtual ~EmptyAlgorithm()
	{
	}
		
	void run(unsigned int numSteps = 1) override;

	
	/**
	* @brief does a dynamic cast of the given GMMalgorithm to EMforGMM
	* @return NULL if the Algorithm is not a EMforGMM instance
	*/
	static EmptyAlgorithm* toEmptyAlgorithm(Algorithm* a);
	
};

/**
 *@brief: description 
 */
class EmptyAlgorithmAD : public AlgoDescription
{
public:
	static const std::string CLASSTAG;
	EmptyAlgorithmAD();
	virtual Algorithm* create(const commonutil::DataSet& input) const;
	
	static void parseDescriptionString(std::vector<AlgoDescription*>& descriptions, const std::string descriptionstring, const bool first = false);
	
};


#endif
