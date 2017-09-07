#ifndef GMMALGORITHM_H
#define GMMALGORITHM_H

#include "../../base/parameters.h"
#include "../../base/commonutil.h"
#include "../../ioutil/output_handler.h"

/**
* @brief Abstract base class for algorithms.
*/
class Algorithm
{
public:
	/**
	* @brief Constructor initializing input and configuring the algorithm.
	* @remarks Remember to init the mixture model.
	*/
	Algorithm(commonutil::DataSet const&, bool = false);
	
	Algorithm(const Algorithm&);
	Algorithm& operator=(const Algorithm&);
	virtual ~Algorithm()
	{
	}

	bool hasChanged() const;
	Parameters const& getGMMDesc() const;
	double getRuntime() const;

	virtual void setInput(commonutil::DataSet const&);
	virtual void init(Parameters const&);
	virtual void run(unsigned int numSteps = 1) = 0;
	
	virtual void finalVerboseness();
	virtual void outputStatistics(unsigned int round, OutputHandler& outputhandler);
	inline bool const errorHandlingTriggered() const { return this->error_handled; };

protected:
	commonutil::DataSet const* input;
	Parameters desc;
	
	double runtime = 0;
	bool change = false;
	bool verbose = false;
	inline void errorHandling(std::string const& errormsg)
	{ 
		this->error_handled = true; 
		if(this->verbose)
			std::cout << errormsg << std::endl;
		
	};
	
private:
	bool error_handled = false;
};


#endif
