#ifndef _ESCONFIG_HPP
#define _ESCONFIG_HPP

#include <boost/random.hpp>
namespace espresso {
  // define this to "float" if you want to have single precision
  typedef double real;
  
  // define this to "long long" if you need longer integers
  typedef int longint;
  
  static const real ROUND_ERROR_PREC = 1.0e-14;
  
  // RNG config
  typedef boost::lagged_fibonacci607 RNGType;
  // If you REALLY need speed!
  //typedef boost::rand48 RNGType;
}

#endif
