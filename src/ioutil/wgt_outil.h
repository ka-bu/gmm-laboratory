#ifndef WGT_OUTIL_H
#define WGT_OUTIL_H

#include "../base.h"
#include "ioutil.h"

#include <vector>

/**
 * @brief 
 */
class WgtOUtil : public  IOUtil
{
public:
	WgtOUtil();
	virtual void header();
	void store(size_t round, Vector const& w1, Vector const& w2);
};

#endif