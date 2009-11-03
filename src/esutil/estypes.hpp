#ifndef _ESTYPES_HPP
#define _ESTYPES_HPP

// define this to "float" if you want to have single precision
typedef double real;

// define this to "long long" if you need longer integers
typedef int longint;

// OBSOLETE
typedef int integer;

static const real ROUND_ERROR_PREC = 1.0e-14;

// RNG config
#include <boost/random.hpp>
typedef boost::lagged_fibonacci607 RNGType;
// If you REALLY need speed!
//typedef boost::rand48 RNGType;


#endif
