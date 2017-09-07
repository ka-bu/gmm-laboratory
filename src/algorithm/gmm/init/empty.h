#ifndef EMPTY_H
#define EMPTY_H

#include "../../../base/parameters.h"
#include "../../description/init_description.h"
#include "adaptive_lloyds_means2gmm.h"

#include <vector>
#include "../utils/gmmutil.h"

namespace initializer
{
	Parameters emptyGMM();
}

/**
 * @brief Descriptor for Initialization
 */
class EmptyID : public InitDescription
{
public:
	const static std::string CLASSTAG;
	EmptyID();
	Parameters compute(commonutil::DataSet const&, unsigned int) const;
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring);
};

#endif
