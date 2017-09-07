#include "diff_outil.h"

DiffOUtil::DiffOUtil(): IOUtil(".diff")
{

}


void DiffOUtil::header()
{
	this->filestream << "round" << SEPARATOR << "time";
	for (size_t i=1; i<=k; ++i)
		this->filestream << SEPARATOR << "weight" << i;
	for (size_t i=1; i<=k; ++i)
		this->filestream << SEPARATOR << "mean" << i;
	for (size_t i=1; i<=k; ++i)
		this->filestream << SEPARATOR << "covar" << i;
	
	if(k==0)
		gmmlab_throw("You should have specified parameter k when calling DiffOUtil::open()")
}

void DiffOUtil::store(std::size_t round, const Vector& diffs)
{
	if(!is_open())
		gmmlab_throw("File not open.");
	this->filestream << round;
	for (idx_type i=0; i<diffs.size(); ++i)
		this->filestream << SEPARATOR << diffs[i];
	this->filestream << std::endl;
}

void DiffOUtil::storeNA(std::size_t round, idx_type diffsize)
{
	if(!is_open())
		gmmlab_throw("File not open.");
	this->filestream << round;
	for (idx_type i=0; i<diffsize; ++i)
		this->filestream << SEPARATOR << NA;
	this->filestream << std::endl;
}
