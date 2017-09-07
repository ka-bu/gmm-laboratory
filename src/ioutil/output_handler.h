#ifndef OUTPUT_HANDLER
#define OUTPUT_HANDLER

#include "diff_outil.h"
#include "dist_outil.h"
#include "data_ioutil.h"
#include "cost_outil.h"
#include "wgt_outil.h"
#include "time_outil.h"


class OutputHandler
{
public: 
	DistOUtil dist;
	DiffOUtil diff;
	WgtOUtil wgt;
	CostOUtil cost;
	CostOUtil cost_diff; // cost of diffalgo in DIFF mode
	DataIOUtil dataset;
	TimeOUtil time;
	DiffOUtil difftruth;
	
	void init(std::string const& outputDir, unsigned int k, bool truth, std::string const& filename_algo1, std::string const& filename_algo2 = "NA", std::string const& filename_both = "NA");
	void close();
	void initTimeOutput(std::string const& outputDir, std::string const& outputFile);
	void closeTimeOutput();
	
private:
	std::string createOutputFileInFolder(std::string const& outputFile, std::string const& outputDir, std::string const& subDir = std::string());	
};


#endif


