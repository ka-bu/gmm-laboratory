#ifndef BASE_H
#define BASE_H

#define NDEBUG

#include <cassert>
#include <cstddef>
#include <Eigen/Dense>
#include <random>
#include <stdint.h>
#include <iostream>

typedef double fp_type;
typedef Eigen::MatrixXd Matrix;
typedef Eigen::VectorXd Vector;
typedef Eigen::DenseIndex idx_type;


inline void __gmmlab_throw(const std::string message, const char* filename, const std::size_t lineNumber)
{
	std::ostringstream stream;
	stream << "ERROR: " << message << "\nin file " << filename << "  at line " << lineNumber << "\n";
	throw std::runtime_error(stream.str());
}

inline void __gmmlab_warning(const std::string message, const char* filename, const std::size_t lineNumber)
{
	std::cout << "WARNING: " << message << "\nin file " << filename << "  at line " << lineNumber << std::endl;
}


#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#define gmmlab_throw(message) \
{__gmmlab_throw(message,__FILE__,__LINE__);}

#define gmmlab_warning(message) \
{__gmmlab_warning(message,__FILE__,__LINE__);}

#endif
