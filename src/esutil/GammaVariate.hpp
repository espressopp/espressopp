// ESPP_CLASS
#ifndef _ESUTIL_GAMMAVARIATE_HPP
#define _ESUTIL_GAMMAVARIATE_HPP
#include <boost/random.hpp>
#include "types.hpp"

namespace espresso {
  namespace esutil {
    using namespace boost;
      /** This class generates gamma distributed random
	  variates. The class also keeps a shared pointer to the boost
	  RNG object, so that it is not destroyed if the espresso RNG
	  object dies. */
    class GammaVariate
      : variate_generator< RNGType&, gamma_distribution< real > >
    {
      typedef gamma_distribution< real > DistType;
      typedef variate_generator< RNGType&, DistType > Super;

      /// store the shared pointer to the RNG
      shared_ptr< RNG > rng;

    public:
      GammaVariate(shared_ptr< RNG > _rng,
		    const int alpha = 1,
		    const real beta = 1.0);
      using Super::operator();
      static void registerPython();
    };
  }
}
#endif
