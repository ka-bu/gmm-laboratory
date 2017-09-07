#ifndef COMMONUTIL_H
#define COMMONUTIL_H

#include "../base.h"

#include <iostream>

namespace commonutil
{	
	struct CompleteData
	{
		Vector point;
		idx_type source;
	};
	
	struct DataSet
	{
		Matrix points;
		Vector weights;
		DataSet(){};
		/**
		 * returns the number of columns. Each column of points contains one data point. 
		 */
		idx_type n() const
		{
			return points.cols();
		}
		/**
		 * returns the number of rows. The i-th row of points contains the i-th dimension (i.e. the i-th component of each data point). 
		 */
		idx_type d() const
		{
			return points.rows();
		}
	};
	
	/**
	 * Erases the i-th element from a vector.
	 * @param v vector
	 * @param i the index
	 */
	void erase(Vector&, idx_type);
	
	/**
	 * Erases the i-th column from a matrix.
	 * @param m matrix
	 * @param i the index
	 */
	void eraseColumn(Matrix&, idx_type);
	
	/**
	 * fills the given vector by real values. The i-th entry is drawn uniformly at random from [min(i), max(i)].
	 */
	void fill(Vector& vector, std::mt19937& gen, Vector min, Vector max);
	
	/**
	 * fills the given vector by real values that are drawn uniformly at random from [min, max].
	 */
	void fill(Vector& vector, std::mt19937& gen, fp_type min, fp_type max);
	
	/**
	 * fills the given vector by real values that are drawn uniformly at random from [min, max].
	 */
	void fill(Matrix& matrix, std::mt19937& gen, fp_type min, fp_type max);
	
	/**
	 * returns true iff there does *not* exist an index i where p1(i)=p2(i). p1 and p2 must have same length.
	 */
	bool componentWiseDifferent(Vector, Vector);
	
	void appendData(DataSet& dataset, DataSet const& additionalData);
	
	
	void get_uniform_subset(DataSet const& dataset, DataSet& samples, idx_type number, std::mt19937& gen);


	template<typename T> const std::string toString(const T&);
	
	// template<typename T> const T fromString(const std::string&);
	
	template<typename RndEngine> const idx_type randomIndex(Vector const& weights, RndEngine& re);
	
	std::vector<idx_type> sort_indices(Vector const& vector, bool increasing = true);
	
	std::vector< idx_type > nth_element_indices(const Vector& vector, idx_type n, bool increasing = true);
	
	void map_indices(std::vector<idx_type>& indices, std::vector<idx_type> const& mapping);

	
}

std::ostream& operator<<(std::ostream&, commonutil::DataSet const&);

template<class T> const std::string commonutil::toString(const T& t)
{
     std::ostringstream stream;
     stream << t;
     return stream.str();
}

// template<class T> const T commonutil::fromString(const std::string& s)
// {
//      std::istringstream stream (s);
//      T t;
//      stream >> t;
//      return t;
// }

template<typename RndEngine> const idx_type commonutil::randomIndex(Vector const& weights, RndEngine& re)
{
	idx_type size = weights.size();
	if (size==0)
		return 0;

	fp_type sum = weights.sum();
	assert(sum>=0);

	std::uniform_real_distribution<> urd(0, sum);
	fp_type r = urd(re);

	std::size_t index = 0;
	fp_type w = weights[0];
	while (r>w&&index<size-1)
	{
		r -= w;
		w = weights[++index];
	}
	return index;
}

#endif
