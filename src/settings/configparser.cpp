#include "configparser.h"
#include "settings.h"

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <exception>

namespace po = boost::program_options;

void configparser::parseCommmonConfiguration(po::options_description& config)
{

	config.add_options()

	// algo
	("intermediate,i", po::value<bool>(&commonSettings().intermediate),"compute intermediate results")	
	("costmeasure", po::value<std::string>(&commonSettings().costmeasure), "cost measure (default: NLL)") 
	
	// input
	("data", po::value<std::string>(&commonSettings().datafile), "data descriptor file")
	("init", po::value<std::string>(&commonSettings().initfile), "init descriptor file")
	("algo", po::value<std::string>(&commonSettings().algofile), "algo descriptor file")
	
	// output
	("output", po::value<std::string>(&commonSettings().outputDir), "output directory")
	("csv", po::value<bool>(&commonSettings().csv), "generate csv file")
	("dataset", po::value<bool>(&commonSettings().dataset), "generate dataset file")
	("wgt", po::value<bool>(&commonSettings().wgt), "generate wgt file")
	("dist", po::value<bool>(&commonSettings().dist), "generate dist file")
	("cls", po::value<bool>(&commonSettings().cls),"append to cls file")
	("difftruth", po::value<bool>(&commonSettings().difftruth), "generate diff file for comparison with truth")

	// output
	("verbose,v", po::value<bool>(&commonSettings().verbose), "print detailed information")
	("printstats", po::value<bool>(&commonSettings().printStatistics), "generate diff file for comparison with truth")
	;
}

void configparser::finalizeCommonConfiguration()
{

}

void configparser::parseGMMLabConfiguration(int argc, char* argv[])
{
	std::cout << "command line: ";

	for (int i = 0; i < argc; ++i)
		std::cout << argv[i] << " ";

	std::cout << std::endl;

	std::string config_file;
	
	// Declare a group of options that will be
	// allowed only on command line
	po::options_description generic("Generic options");
	generic.add_options()
	("version", "print version string")
	("help,h", "produce help message")
	("config", po::value<std::string>(&config_file)->default_value("gmmlab.cfg"),"name of configuration file.")
	;
	
	// Declare a group of options that will be allowed both on command line and in config file
	po::options_description config("Configuration");
	config.add_options()
	("mode,m", po::value<std::string>(&gmmlabSettings().mode), "mode of operation")
	("costs,c",  po::value<bool>(&gmmlabSettings().computeCosts),  "omit cost computation")
	("output-file,o", po::value<std::string>(&gmmlabSettings().outputFile), "output file")
	("fullscreen,f", po::value<bool>(&displaySettings().fullscreen), "use fullscreen instead of window")
	("colored_cluster", po::value<bool>(&displaySettings().coloredCluster), "color points depending on their cluster")
 	("ndisplay", po::value<bool>(&commonSettings().ndisplay), "disable gui")
	;
	parseCommmonConfiguration(config);

	
	// Hidden options, will be allowed both on command line and
	// in config file, but will not be shown to the user.
	po::options_description hidden("Hidden options");
	hidden.add_options()
	;

	po::options_description cmdline_options;
	cmdline_options.add(generic).add(config).add(hidden);

	po::options_description config_file_options;
	config_file_options.add(config).add(hidden);

	po::options_description visible("Allowed options");
	visible.add(generic).add(config);

	po::positional_options_description p;
	p.add("file", 1).add("load", 1);
	
	try
	{
		po::variables_map vm;
		store(po::command_line_parser(argc, argv).
		options(cmdline_options).positional(p).run(), vm);
		notify(vm);

		std::ifstream ifs(config_file.c_str());

		if (ifs)
		{
			store(parse_config_file(ifs, config_file_options), vm);
			notify(vm);
		}
		
		finalizeCommonConfiguration();

		if (vm.count("help"))
		{
			std::cout << visible << "\n";
			exit(0);
		}

		if (vm.count("version"))
		{
			std::cout << "GMM Laboratory, version 0.5\n";
			exit(0);
		}

			
		/*        if (!commonSettings().csv && !gmmlabSettings().lab && !gmmlabSettings().dist
				&& !gmmlabSettings().cls && !gmmlabSettings().gmm
				&& vm.count("output"))
				std::cout << "WARNING: output filename was given without specifying --csv, --lab, --dist, --cls or --gmm" << std::endl;
		*/
		if (commonSettings().verbose)
			std::cout << "Verbose output is ON." << std::endl << std::endl;

	}
	catch (std::exception& e)
	{
		std::cout << e.what() << "\n";
		exit(1);
	}
	
	
}

void configparser::parseTestLabConfiguration(int argc, char* argv[])
{
    std::cout << "command line: ";
    for (int i = 0; i < argc; ++i)
        std::cout << argv[i] << " ";
    std::cout << std::endl;

    std::string config_file, mode;

    // Declare a group of options that will be
    // allowed only on command line
    po::options_description generic("Generic options");
    generic.add_options()
    ("version", "print version string")
    ("help,h", "produce help message")
    ("config", po::value<std::string>(&config_file)->default_value("testlab.cfg"),"name of configuration file.")
    ;
    
    // Declare a group of options that will be allowed both on command line and in config file
    po::options_description config("Configuration");
    config.add_options()
    // testlab-settings
    ("mode,m", po::value<std::string>(&testlabSettings().mode), "mode of operation")
    ("time", po::value<bool>(&testlabSettings().time), "generate time file")
    ;
    parseCommmonConfiguration(config);

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden("Hidden options");
    hidden.add_options()
    ;

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config).add(hidden);

    po::options_description config_file_options;
    config_file_options.add(config).add(hidden);

    po::options_description visible("Allowed options");
    visible.add(generic).add(config);

    po::positional_options_description p;
    p.add("data", 1).add("init", 1).add("algo", 1).add("diff", 1);

    try
    {
        po::variables_map vm;
        store(po::command_line_parser(argc, argv).
              options(cmdline_options).positional(p).run(), vm);
        notify(vm);

        std::ifstream ifs(config_file.c_str());

        if (ifs)
        {
            store(parse_config_file(ifs, config_file_options), vm);
            notify(vm);
        }
        
        finalizeCommonConfiguration();

        if (vm.count("help"))
        {
            std::cout << visible << "\n";
            exit(0);
        }

        if (vm.count("version"))
        {
            std::cout << "Test Laboratory, version 0.1\n";
            exit(0);
        }

	commonSettings().ndisplay = true;

        if (commonSettings().verbose)
        {
            std::cout << "Verbose output is ON." << std::endl << std::endl;
        }
    }
    catch (std::exception& e)
    {
        std::cout << e.what() << "\n";
        exit(1);
    }
    
}
