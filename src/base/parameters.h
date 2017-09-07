#ifndef GMMDESC_H
#define GMMDESC_H

#include "../base.h"

#include <vector>

/**
 * @brief Data structure for the parameters of an mpe mixture.
 */
struct Parameters
{
public:
	Vector weights;
	Matrix means;
	std::vector<Matrix> covariances;
	Vector kappa;
	
	virtual bool operator==(Parameters const& rhs) const;

	virtual bool empty() const;
	
	idx_type components() const;
};

std::ostream& operator<<(std::ostream& os, Parameters const& desc);

#endif
