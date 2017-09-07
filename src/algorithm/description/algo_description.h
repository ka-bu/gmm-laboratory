#ifndef ALGO_DESCRIPTION_H
#define ALGO_DESCRIPTION_H

#include "../../base/description.h"
#include "../base/algorithm.h"

class AlgoDescription : public Description
{
public:
	virtual Algorithm* create(commonutil::DataSet const&) const = 0;
	
 	virtual ~AlgoDescription()
	{
	}

	void setNext(AlgoDescription*);
	AlgoDescription* getNext() const;
	unsigned int getSteps() const;
	
	static void parseDescriptionString(std::vector<AlgoDescription*> descriptions, std::string descriptionstring, bool first = false);
	static std::vector< AlgoDescription* > createAlgoDescriptions(std::string filename);
	
protected:
	
	AlgoDescription* next = NULL;
	unsigned int steps = 0;
	
};

#endif