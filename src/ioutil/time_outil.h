#ifndef TIME_OUTIL_H
#define TIME_OUTIL_H

#include "../base.h"
#include "ioutil.h"

#include <vector>

/**
 * @brief 
 */
class TimeOUtil : public IOUtil
{
public:
	TimeOUtil();
	virtual void header();
	void store(std::string dataset, std::string init, std::string algo, double runtime);

};

#endif