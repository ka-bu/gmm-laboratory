#include "gmmlab.h"

#include "base.h"
#include "algorithm/base/samplingalgorithm.h"
#include "settings/settings.h"
#include "base/similarities.h"


#include <iomanip>
#include "algorithm/kmeans/utils/kmeansutil.h"


const std::string GMMLab::COST_GAUSSIAN = "GAUSSIAN"; // negative log-likelihood
const std::string GMMLab::COST_KMEANS = "KMEANS"; // sum of squared errors


GMMLab::GMMLab(bool display, bool computeCosts, std::string outputFile) : display(display), computeCosts(computeCosts), outputFile(outputFile)
{
	std::cout << std::endl << std::endl
		<< "GMMLab - Run automated tests of gmm algorithms" << std::endl
		<< "------------------------------------------------" << std::endl << std::endl;

//	clog.clear(ios_base::badbit);
//	cout << setiosflags(ios::fixed) << setprecision(3);
	std::cout << std::setiosflags(std::ios::scientific)
		<< std::setiosflags(std::ios::showpoint)
		<< std::setprecision(std::numeric_limits<fp_type>::digits10);
}

GMMLab::~GMMLab()
{
	this->close();
}

fp_type GMMLab::getCosts(commonutil::DataSet const& data, Parameters const& gmmdesc)
{  		
	
	if(commonSettings().costmeasure == COST_GAUSSIAN)
		return gmmutil::nll(input, gmmdesc);
	else if(commonSettings().costmeasure == COST_KMEANS)
		return kmeansutil::kmeanscost(input, gmmdesc.means);
	else
		gmmlab_throw("undefined costmeasure")
}

void GMMLab::close()
{
	output.close();
	output.closeTimeOutput();
}



void GMMLab::handleSolution(Algorithm* algorithm, std::size_t round, Parameters const& gmmdesc, double time, bool hasChanged, Algorithm* diff_algorithm)
{
	
	//std::cout << "round = " << round << std::endl
	//          << "\ngmmdesc = " << gmmdesc << std::endl << std::endl;
	
	// display? then store the samples
	if(display)
	{
		solutions.push_back(gmmdesc);
		std::ostringstream stream;
		stream << "round " << round;
		tags.push_back(stream.str());
		runtimes.push_back(time);

		if(algorithm != NULL)
		{
			SamplingAlgorithm* sa = SamplingAlgorithm::toSamplingAlgorithm(algorithm);
			if (sa != NULL)
				sampleSets.push_back(sa->getSamples());
		}
		
		
		// uncomment the following "if" if the samplesets shall be displayed *before* the new gmm
		if(algorithm == NULL)
		{
			commonutil::DataSet emptyset;
			sampleSets.push_back(emptyset);
		}
	}
	
	if(output.difftruth.is_open() && !truth.empty())
	{
		unsigned int k = gmmdesc.components();
		Vector diffs(3*k+1);
		
		if(algorithm == NULL)
			diffs = Vector::Zero(3*k+1);
		else
		{
			diffs[0] = algorithm->getRuntime(); //fabs(algo1.getRuntime()-algo2.getRuntime());

			for (idx_type i=0; i<k; ++i)
				diffs[1+i] = fabs(algorithm->getGMMDesc().weights[i]-truth.weights[i]);
			for (idx_type i=0; i<k; ++i)
				diffs[1+k+i] = (algorithm->getGMMDesc().means.col(i)-truth.means.col(i)).norm();
			for (idx_type i=0; i<k; ++i)
				diffs[1+2*k+i] = (algorithm->getGMMDesc().covariances[i]-truth.covariances[i]).norm();
		}
		
		output.difftruth.store(round, diffs);
	}

	
	if (output.dist.is_open())
	{
		GMMDis dis;
		Vector x = dis(gmmdesc, truth);
		output.dist.store(round, time, x);
	}
	


	fp_type cost;	
	if (output.cost.is_open() || computeCosts)
	{
	
		// compute costs (only if the solution has changed)
		if(cost_initialized)
		{
			if(!hasChanged)
				cost = last_cost;
			else
				cost = getCosts(input, gmmdesc);
		}
		else
		{
			cost = getCosts(input, gmmdesc);
			//cost_initialized = true;
		}
		
		// output: display
		if(computeCosts)
			costs.push_back(cost);
		
		// output: cost file
		if(output.cost.is_open())
			output.cost.store(round, time, cost);

		// minimum?
		if(!cost_initialized || (min_cost > cost))
		{
			min_cost = cost;
			min_cost_runtime = time;
			min_cost_round = round;
			cost_initialized = true;
		}
		
		last_cost = cost;
	}
	
	if(output.cost_diff.is_open())
	{
		if(diff_algorithm != NULL)
		{
			Parameters gmmdesc = diff_algorithm->getGMMDesc();
			output.cost_diff.store(round, time, getCosts(input, gmmdesc));
		}
		else
		{
			output.cost_diff.store(round, time, cost);
		}
	}
	
	
	if(algorithm!= NULL)
		algorithm->outputStatistics(round, output);
}




double GMMLab::runAlgoSequence(AlgoDescription const* description, Parameters initialSolution, double initialTime,
	std::string const& outputDir, std::string const& outputFile)
{	
	// Initialize output
	output.init(outputDir, initialSolution.components(), !truth.empty(),  outputFile);

	if (commonSettings().dataset)
		output.dataset.store(input);
	
	// compute maximum spread
	Vector spreads = (input.points.rowwise().maxCoeff()
					 -input.points.rowwise().minCoeff());
					 
	if (commonSettings().verbose)
		std::cout << "maximum spread is " << spreads.maxCoeff() << std::endl;


	// Initialisation
	double totaltime = initialTime;
	Parameters last = initialSolution;

	
	// handle initial solution
	cost_initialized = false;
	if(!last.empty())
		handleSolution(NULL, 0, last, totaltime, true);

	// start computation of new algorithm sequence
	std::cout << "starting new algorithm sequence:" << std::endl;
			  
	AlgoDescription const* nextAlgo = description;	
	std:size_t roundOffset = 0;
	std::vector<std::size_t> finalRounds;
	std::vector<double> finalTimes;
	while (nextAlgo!=NULL)
	{
		std::cout << "creating algorithm from descriptor: "
				  << nextAlgo->getNameTag() << " (last change: " << __TIMESTAMP__ << ")"<< std::endl;
				  
		Algorithm* algorithm = nextAlgo->create(input);
		algorithm->init(last);
		
		// run current algorithm
		for (std::size_t i=0; i<nextAlgo->getSteps(); ++i)
		{
// 			if(commonSettings().verbose)
				std::cout << "> round " << i << std::endl;

			algorithm->run();
			if(commonSettings().intermediate)
				handleSolution(algorithm, roundOffset+i+1, algorithm->getGMMDesc(),
					algorithm->getRuntime()+totaltime, algorithm->hasChanged());
		}
		
		if(!commonSettings().intermediate)
			handleSolution(algorithm, roundOffset+nextAlgo->getSteps(), algorithm->getGMMDesc(),
					algorithm->getRuntime()+totaltime, algorithm->hasChanged());
			
		// final verboseness
		algorithm->finalVerboseness();
		
		// Save results and delete algorithm
		last = algorithm->getGMMDesc();
		roundOffset += nextAlgo->getSteps();
		totaltime += algorithm->getRuntime();
		finalRounds.push_back(roundOffset);
		finalTimes.push_back(totaltime);
		std::cout << "   computation took " << algorithm->getRuntime() << " seconds." << std::endl;
// 		std::cout << "computed solution" << std::endl << last << std::endl;


		delete algorithm;

		nextAlgo = nextAlgo->getNext();
	}
	std::cout << std::endl;

	if (output.cost.is_open())
	{
		// Store statistics markers
		output.cost.storeMarkerVector(finalRounds, finalTimes, min_cost_round, min_cost_runtime);

		// Print overall results
		if (!truth.empty())
			std::cout << "         cost of truth = " << getCosts(input, truth) << std::endl;
	}

	if(computeCosts)
	{
		std::cout << "cost of final solution = " << last_cost << std::endl;
		std::cout << "     cheapest solution = " << min_cost << " after " << min_cost_runtime << " seconds" << std::endl;
	}
	
	std::cout << std::endl;

		
	// Finalize output
	output.close();
	
	return totaltime-initialTime;
}







void GMMLab::runDataGenerationOnly(std::vector<DataDescription*> const& dataDescriptors)
{
	idx_type dataIndex = 0;
	
	if(dataDescriptors.size() != 1)
		gmmlab_throw("Define exactly one data descriptions");
	
	// clear input
	input.weights.resize(0);
	input.points.resize(0,0);
	
	std::cout << ">" << std::endl
		  << ">> retrieving new data set from descriptor: "
		  << dataDescriptors.at(dataIndex)->tag() << " (last change: " << __TIMESTAMP__ << ")" << std::endl
		  << ">" << std::endl;

	// load or generate data from descriptor		
	truth = dataDescriptors.at(dataIndex)->retrieve(input);
}

void GMMLab::runDefaultTests(std::vector<DataDescription*>  const&  dataDescriptors, std::vector<InitDescription*> const&  initDescriptors, std::vector<AlgoDescription*> const&  algoDescriptions)
{
	output.initTimeOutput(commonSettings().outputDir,outputFile.empty()?"time":outputFile); 
	
	// loop over all data sets
	for (unsigned int dataIndex=0; dataIndex<dataDescriptors.size(); ++dataIndex)
	{
		// clear input
		input.weights.resize(0);
		input.points.resize(0,0);
	
		std::cout << ">" << std::endl
				  << ">> retrieving new data set from descriptor: "
				  << dataDescriptors.at(dataIndex)->tag() << " (last change: " << __TIMESTAMP__ << ")" << std::endl
				  << ">" << std::endl;

		// load or generate data from descriptor		
		truth = dataDescriptors.at(dataIndex)->retrieve(input);

		if (input.weights.size()==0)
		{
			std::cout << "   skipping empty data set." << std::endl << std::endl;
			continue;
		}
		
		idx_type n = input.points.cols();
		idx_type d = input.points.rows();
		
		
		std::cout << "Using " << d << " dimensions of the " << n << " input points." << std::endl;
		
		// loop over all values of k
		for (DataDescription::kIterator kIter=dataDescriptors[dataIndex]->first_k();
				kIter!=dataDescriptors[dataIndex]->last_k(); ++kIter)
		{
			unsigned int k = *kIter;
			std::cout << std::endl << "> computing solutions with k=" << k << std::endl << std::endl;
										  
			std::stringstream outDirStream;
			outDirStream << commonSettings().outputDir << "/"
						 << dataDescriptors[dataIndex]->tag()
						 << "_k" << k;
				
			// loop over all initial solutions
			for (unsigned int initIndex=0; initIndex<initDescriptors.size(); ++initIndex)
			{
		
				std::cout << "computing initial solution from descriptor: "
						  << initDescriptors[initIndex]->getNameTag() << " (last change: " << __TIMESTAMP__ << ")"<< std::endl;
							  
				clock_t t = clock();
				Parameters initial = initDescriptors[initIndex]->compute(input, k);
				double inittime = secondsSince(t);
				std::cout << "   computed in " << inittime << " seconds." << std::endl << std::endl;
				
// std::cout << "initial solution = " << initial << std::endl;
	
				std::cout << "computing algorithms..." << std::endl;
				
				for (int algoIndex=0; algoIndex<algoDescriptions.size(); ++algoIndex)
				{
					std::stringstream outFileStream;
					if(outputFile.empty())
					{
						outFileStream << initDescriptors.at(initIndex)->getNameTag() << "__";
						AlgoDescription* ad = algoDescriptions.at(algoIndex);
						while (ad!=NULL)
						{
							outFileStream << ad->getNameTag();
							ad = ad->getNext();
							if (ad!=NULL)
								outFileStream << "__";
						}
					}

					try
					{
						double algotime = runAlgoSequence(algoDescriptions.at(algoIndex), initial, inittime,
										outputFile.empty()?outDirStream.str():commonSettings().outputDir, 
										  outputFile.empty()?outFileStream.str():outputFile);
						
std::cout << "algotime= " << algotime << std::endl;
std::cout << "totalAlgotime = " << inittime+algotime << std::endl;

						if (output.time.is_open())
						{
			
							std::stringstream dataset;
							dataset << dataDescriptors.at(dataIndex)->tag() << "_k" << k; //<< ", "<<initDescriptors.at(initIndex)->getNameTag();
							output.time.store(dataset.str(), initDescriptors.at(initIndex)->getNameTag(), algoDescriptions.at(algoIndex)->getNameTag(), algotime);
						}
					}
					catch (...)
					{
						std::cout << "runDefaultTests() - EXCEPTION CAUGHT! Continuing with next algorithm sequence." << std::endl;
						output.close();
					}
				}	
			}
			
			
		}
	}
	output.closeTimeOutput();
	
}



