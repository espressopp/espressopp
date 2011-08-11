// ESPP_CLASS
#ifndef _ESUTIL_NORMALVARIATE_HPP
#define _ESUTIL_NORMALVARIATE_HPP
#include <boost/random.hpp>
#include "types.hpp"

namespace espresso {
  namespace esutil {
    using namespace boost;
      /** This class generates normal distributed random
	  variates. The class also keeps a shared pointer to the boost
	  RNG object, so that it is not destroyed if the espresso RNG
	  object dies. */
    class NormalVariate 
      : variate_generator< RNGType&, normal_distribution< real > > 
    {
      typedef normal_distribution< real > DistType;
      typedef variate_generator< RNGType&, DistType > Super;

      /// store the shared pointer to the RNG
      shared_ptr< RNG > rng;

    public:
      NormalVariate(shared_ptr< RNG > _rng,
		    const real mean = 0.0, 
		    const real sigma = 1.0);
      using Super::operator();
      static void registerPython();
    };
  }
}
#endif
