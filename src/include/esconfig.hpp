#ifndef _ESCONFIG_HPP
#define _ESCONFIG_HPP

#include <boost/random.hpp>
#include <limits>
namespace espresso {
  // define to "float" for single precision (i.e. typedef float real;)
  // define to "double" for double precision (i.e. typedef double real;)
  typedef double real;

  static const real infinity = std::numeric_limits< real >::infinity();
  static const real ROUND_ERROR_PREC = std::numeric_limits< real >::epsilon();

  // define this to "long long" if you need longer integers
  typedef int longint;
  
  // RNG config
  typedef boost::lagged_fibonacci607 RNGType;
  // If you REALLY need speed use the line below instead
  //typedef boost::rand48 RNGType;
}

#endif
