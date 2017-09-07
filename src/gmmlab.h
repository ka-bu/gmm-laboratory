#ifndef GMMLAB_H
#define GMMLAB_H

#include "data/data_description.h"
#include "algorithm/description/algo_description.h"
#include "algorithm/description/init_description.h"
#include "ioutil/output_handler.h"

class GMMLab
{
public:
	
	static const std::string COST_GAUSSIAN;
	static const std::string COST_KMEANS;
	
	commonutil::DataSet input;
	Parameters truth;
	std::vector<commonutil::DataSet> sampleSets;
	std::vector<double> runtimes;
	std::vector<Parameters> solutions;
	std::vector<std::string> tags;
	std::vector<fp_type> costs;
	
	GMMLab(bool display, bool computeCosts, std::string outputFile = std::string());
	~GMMLab();
	
	void runDataGenerationOnly(std::vector<DataDescription*> const& dataDescriptors);
	void runDefaultTests(std::vector<DataDescription*>  const&  dataDescriptors, std::vector<InitDescription*> const&  initDescriptors, std::vector<AlgoDescription*> const&  algoDescriptions);
	void runDiffTests(std::vector<DataDescription*> const& dataDescriptors, std::vector<InitDescription*>  const&  initDescriptors, std::vector<AlgoDescription*>  const& algo1Descriptors, std::vector<AlgoDescription*>  const& algo2Descriptors);
	
	void close();
	
private:
	
	bool display = false;
	
	
	bool computeCosts;
	
	bool cost_initialized;
	fp_type min_cost;
	fp_type last_cost;
	double min_cost_runtime = 0;
	std::size_t min_cost_round = 0;

	OutputHandler output;
	std::string outputFile = "";
	
	
	/**
	 * computes cost of the given solution depending on cost measure stored in commonutils.
	 */
	fp_type getCosts(commonutil::DataSet const& data, Parameters const& gmmdesc);
	
	/**
	 * handles the results, i.e. stores outputs, computes differences between solutions and so on.
	 * Remember to set cost_initialized = false before the first time handleSolutions is called during an experiment.
	 */
	void handleSolution(Algorithm* algorithm, std::size_t round, Parameters const& gmmdesc, double time, bool hasChanged, Algorithm* diff_algorithm = NULL);
	
	/**
	 * 
	 */
	double runAlgoSequence(AlgoDescription const* description, Parameters initialSolution, double initialTime,	std::string const& outputDir, std::string const& outputFile);
	
	/**
	 * 
	 */
	void computeDiff(Algorithm const& algo1, Algorithm const& algo2, idx_type k, std::size_t round);
	
	/**
	 * 
	 */
	void diffAlgoSequence(AlgoDescription const* ad1, AlgoDescription* ad2, unsigned int k, Parameters initialSolution, double initialTime, std::string const& outputDir, std::string const& filename_algo1, std::string const& filename_algo2, std::string const& filename_both_algos);
	
	
};

inline double secondsSince(clock_t t)
{
	return double(clock()-t)/CLOCKS_PER_SEC;
};

#endif
