#include "algorithm.h"

#include <algorithm>
#include <ctime>
#include <random>


Algorithm::Algorithm(commonutil::DataSet const& data, bool v) : input(&data), runtime(0),
	change(false), verbose(v)
{
	assert(data.points.cols()==data.weights.size());
}

Algorithm::Algorithm(const Algorithm& rhs) : input(rhs.input), desc(rhs.desc),
	runtime(rhs.runtime), change(rhs.change), verbose(rhs.verbose)
{	
}

Algorithm& Algorithm::operator=(const Algorithm& rhs)
{
	this->input = rhs.input;
	this->desc = rhs.desc;
	this->runtime = rhs.runtime;
	this->change = rhs.change;
	return *this;
}

void Algorithm::setInput(commonutil::DataSet const& data)
{
	this->input = &data;
	if(data.points.cols() == 0)
		gmmlab_throw("DataSet is empty.");
}

void Algorithm::init(Parameters const& dsc)
{
	this->desc = dsc;
	this->change = true;
// 	if(this->desc.empty())
// 		gmmlab_throw("Initial solution is empty.");
}

bool Algorithm::hasChanged() const
{
	return this->change;
}

Parameters const& Algorithm::getGMMDesc() const
{
	return this->desc;
}

double Algorithm::getRuntime() const
{
	return this->runtime;
}

void Algorithm::outputStatistics(unsigned int round, OutputHandler& outputhandler)
{
	// do nothing
}

void Algorithm::finalVerboseness()
{
	// do nothing
}

