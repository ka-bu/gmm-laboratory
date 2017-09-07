#ifndef DIST_OUTIL_H
#define DIST_OUTIL_H

#include "../base.h"
#include "ioutil.h"

#include <vector>

/**
 * @brief 
 */
class DistOUtil : public IOUtil
{
public:	
	DistOUtil();
	virtual void header();
	void store(size_t round, double time, Vector const& dist);

};

#endif