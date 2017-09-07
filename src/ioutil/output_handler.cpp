#include "output_handler.h"

#include "diff_outil.h"
#include "dist_outil.h"
#include "cost_outil.h"
#include "wgt_outil.h"

#include "../settings/settings.h"

#include <boost/filesystem.hpp>

void OutputHandler::close()
{
	dist.close();
	diff.close();
	dataset.close();
	wgt.close();
	cost.close();
	cost_diff.close();
	difftruth.close();
}

std::string OutputHandler::createOutputFileInFolder(std::string const& outputFile, std::string const& outputDir, std::string const& subDir)
{
	std::stringstream pathStream;
	pathStream << outputDir;
	if(!subDir.empty())
		pathStream << "/" << subDir;
	boost::filesystem::path outPath(pathStream.str());
	boost::filesystem::create_directories(outPath);
	std::stringstream filenameStream;
	filenameStream << pathStream.str() << "/" << outputFile;
	return filenameStream.str();
}


void OutputHandler::init(std::string const& outputDir, unsigned int k, bool truth, std::string const& filename_algo1, std::string const& filename_algo2, std::string const& filename_both)
{
	if (commonSettings().csv)
		cost.open(createOutputFileInFolder(filename_algo1, outputDir, "csv"));
	if (commonSettings().dataset)
	{
		size_t pos = outputDir.length();
		while (outputDir.compare(--pos,1,"/"));
		dataset.open(createOutputFileInFolder(outputDir.substr(pos+1, outputDir.size()-pos), outputDir, "dataset"));
	}
	if (commonSettings().dist && truth)
		dist.open(createOutputFileInFolder(filename_algo1, outputDir, "dist"));
	if(commonSettings().wgt)
		wgt.open(createOutputFileInFolder(filename_algo1, outputDir, "wgt"));
	if(commonSettings().difftruth)
		difftruth.open(createOutputFileInFolder(filename_algo1, outputDir, "truth"),k);
}

void OutputHandler::initTimeOutput(const std::string& outputDir, const std::string& outputFile)
{
	if(testlabSettings().time)
	{
		time.open(createOutputFileInFolder(outputFile, outputDir));
	}
}

void OutputHandler::closeTimeOutput()
{
	time.close();
}

