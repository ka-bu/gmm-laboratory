#ifndef INIT_DESCRIPTION_H
#define INIT_DESCRIPTION_H

#include "../../base/description.h"
#include "../base/algorithm.h"

class InitDescription : public Description
{
public:
	virtual Parameters compute(commonutil::DataSet const& input, unsigned int k) const = 0;
		
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring, bool first = false);
	
	static std::vector< InitDescription* > createInitDescriptions(std::string filename);
};

// TBD: move somewhere else...
namespace initializer
{
	/**
	 * throws an exception if k==0 || d==0 || n==0.
	 */
	void check(idx_type k, idx_type d, idx_type n);
}

#endif