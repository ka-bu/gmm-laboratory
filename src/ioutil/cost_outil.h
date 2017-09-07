#ifndef COST_OUTIL_H
#define COST_OUTIL_H

#include "../base.h"
#include "ioutil.h"

#include <vector>

/**
 * @brief 
 */
class CostOUtil : public IOUtil
{
public:
	CostOUtil();
	virtual void header();
	void store(std::size_t round, double runtime, fp_type costs);
	void storeMarker(std::size_t mainrnd, double maintime, std::size_t minrnd, double mintime);
	void storeMarkerVector(std::vector<std::size_t> endRounds, std::vector<double> endTimes, std::size_t minrnd, double mintime);

};

#endif