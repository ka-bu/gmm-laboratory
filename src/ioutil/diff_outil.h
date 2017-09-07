#ifndef DIFF_OUTIL_H
#define DIFF_OUTIL_H

#include "../base.h"
#include "ioutil.h"

#include <vector>

/**
 * @brief 
 */
class DiffOUtil : public IOUtil
{
public:
	DiffOUtil();
	virtual void header();
	void store(size_t round, Vector const& diffs);
	void storeNA(size_t round, idx_type diffsize);

};

#endif