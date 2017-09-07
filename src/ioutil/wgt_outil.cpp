#include "wgt_outil.h"

WgtOUtil::WgtOUtil(): IOUtil(".wgt")
{

}

void WgtOUtil::header()
{
	this->filestream << "round";
	for (size_t i=1; i<=2; ++i)
		for (size_t j=1; j<=k; ++j)
			this->filestream << SEPARATOR << "a" << i << "w" << j;
}

void WgtOUtil::store(size_t round, Vector const& w1, Vector const& w2)
{
	if(!is_open())
		gmmlab_throw("File not open.");
	filestream << round;
	for (idx_type i=0; i<w1.size(); ++i)
		filestream << SEPARATOR << w1[i];
	for (idx_type i=0; i<w2.size(); ++i)
		filestream << SEPARATOR << w2[i];
	filestream << std::endl;
}



