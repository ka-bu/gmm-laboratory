#ifndef DATADESCRIPTOR_H
#define DATADESCRIPTOR_H

#include "../base/parameters.h"
#include "../base/description.h"

#include <vector>
#include "../algorithm/gmm/utils/gmmutil.h"

/**
 * @brief Interface for describing data sets to be loaded or generated
 */
class DataDescription : public Description
{
public:
	typedef std::vector<unsigned int>::const_iterator kIterator;

	virtual Parameters retrieve(commonutil::DataSet&) = 0;
	std::string tag();
	kIterator first_k() const;
	kIterator last_k() const;
	
	static std::vector<DataDescription*> createDataDescriptions(std::string const& filename);
	
protected:
	std::vector<unsigned int> kList;	
};

#endif
