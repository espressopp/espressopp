#ifndef _ESCONFIG_HPP
#define _ESCONFIG_HPP

#include <boost/random.hpp>
#include <limits>
namespace espresso {
  // define this to "float" if you want to have single precision
  typedef double real;
  static const real infinity = std::numeric_limits< real >::infinity();
  static const real ROUND_ERROR_PREC = std::numeric_limits< real >::epsilon();

  // define this to "long long" if you need longer integers
  typedef int longint;
  
  
  // RNG config
  typedef boost::lagged_fibonacci607 RNGType;
  // If you REALLY need speed!
  //typedef boost::rand48 RNGType;
}

#endif
