#ifndef SIMILARITIES_H
#define SIMILARITIES_H

#include "../base.h"
#include "parameters.h"

class MahalanobisDis
{
public:
	fp_type operator()(Vector const&, Matrix const&, Vector const&, Matrix const&);
};

class GMMDis
{
public:
	Vector operator()(Parameters const&, Parameters const&);
};

class AverageLinkage
{
public:
	fp_type operator()(Vector const& sum1, idx_type count1, Vector const& sum2, idx_type count2);
};

#endif