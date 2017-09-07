#include "time_outil.h"

TimeOUtil::TimeOUtil(): IOUtil(".time")
{

}


void TimeOUtil::header()
{
	this->filestream << "dataset" << SEPARATOR << "init" << SEPARATOR << "algo" << SEPARATOR << "avgruntime";
}

void TimeOUtil::store(std::string dataset, std::string init, std::string algo, double runtime)
{
	this->filestream << dataset << SEPARATOR << init << SEPARATOR << algo << SEPARATOR << runtime << std::endl;
}


