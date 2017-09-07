#include "dist_outil.h"

DistOUtil::DistOUtil() : IOUtil(".dist")
{
}

void DistOUtil::header()
{
	this->filestream << "round" << SEPARATOR << "runtime" << SEPARATOR << "weight" << SEPARATOR << "mean" << SEPARATOR << "covariance";
}

void DistOUtil::store(size_t round, double time, Vector const& dist)
{
	if(!is_open())
		gmmlab_throw("File not open.");
	if(dist.size()!=3)
		gmmlab_throw("Cannot store the distance.. wrong number of entries.");
	filestream << round << SEPARATOR << time << SEPARATOR << dist(0)  << SEPARATOR << dist(1) << SEPARATOR << dist(2) << std::endl;
}

