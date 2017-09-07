#ifndef FILE_DESCRIPTION_H
#define FILE_DESCRIPTION_H

#include "data_description.h"

#include <vector>

/**
 * @brief Descriptor for data sets to be loaded from file
 */
class FileDescription : public DataDescription
{
public:
	FileDescription(std::string const&, unsigned int, unsigned int, unsigned int, bool);

	virtual ~FileDescription()
	{
	}
	
	void add_k(unsigned int);
	
	Parameters retrieve(commonutil::DataSet& input);

private:
	std::string prefix;
	unsigned int firstCoord = 0;
	unsigned int numCoords = 0;
	unsigned int blockLength = 1;
	bool normalize = false;
};

#endif
