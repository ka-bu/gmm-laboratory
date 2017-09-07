#ifndef DATA_IUTIL_H
#define DATA_IUTIL_H

#include "../base.h"
#include "ioutil.h"
#include "../base/commonutil.h"

#include <vector>

/**
 * @brief 
 */
class DataIOUtil : public IOUtil
{
public:
	static const size_t AUDIO_FEATUREVECTOR_ENTRY_SIZE;
	static const size_t AUDIO_FEATUREVECTOR_LENGTH;
	static const std::string DATA_FILE_EXTENSION;
	
	static void printStatistis(commonutil::DataSet const& dataset);
	static void runDataGeneration(std::string const& dataDescriptionFile, std::string const& outputDir = std::string());
	static void loadData(commonutil::DataSet& target, std::string const& pattern, unsigned int firstCoord, unsigned int numCoords, unsigned int blockLength, bool normalize, bool noFrameDrop);
	static idx_type loadDataFromFE(commonutil::DataSet& target, std::string const& filename, unsigned int firstCoord, unsigned int numCoords, unsigned int blockLength, bool normalize, bool noFrameDrop);
	static void loadDataFromCSV(commonutil::DataSet& target, std::string const& filename, unsigned int firstCoord, unsigned int numCoords, unsigned int blockLength, bool normalize);
	
	void store(Vector const& entry);
	void store(commonutil::DataSet const& input);
};

#endif