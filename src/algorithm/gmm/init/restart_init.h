#ifndef RESTART_INIT_H
#define RESTART_INIT_H

#include "../../../base/parameters.h"
#include "../../description/init_description.h"

#include <vector>
#include "../utils/gmmutil.h"

namespace initializer
{
	/**
	 * restarts the EM-Algorithm as many times as defined by subruns (if subruns == 0, it will be set to k). 
	 * Each EM instance is initialized via initmethod and runs for substeps steps.
	 * Returns the resulting GMMDesc with the smallest cost. 
	 */
	Parameters restartInit(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, std::string initmethod, idx_type subruns, idx_type substeps, fp_type factor = 1.);
}

/**
 * @brief Descriptor for Initialization
 */
class RestartInitID : public InitDescription
{
public:
	const static std::string CLASSTAG;
	RestartInitID(std::string initmethod, unsigned int subruns, unsigned int substeps, uint32_t seed, fp_type factor = 1.);
	Parameters compute(commonutil::DataSet const&, unsigned int) const;
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring);

private:
	uint32_t seed;
	std::string initmethod;
	unsigned int subruns;
	unsigned int substeps;
	fp_type factor;
};

#endif
