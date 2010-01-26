#ifndef _ESUTIL_UNIFORMONSPHERE_HPP
#define _ESUTIL_UNIFORMONSPHERE_HPP
#include <boost/random.hpp>
#include "types.hpp"

namespace espresso {
  namespace esutil {
    using namespace boost;
    /** This class generates random vectors that are uniformly
	distributed on a sphere. */
    class UniformOnSphere
      : variate_generator< RNGType, uniform_on_sphere< real > > 
    {
      typedef uniform_on_sphere< real > DistType;
      typedef variate_generator< RNGType, DistType > Super;
      
      /// store the shared pointer to the RNG
      shared_ptr< RNGType > rng;
      
    public:
      UniformOnSphere(shared_ptr< RNGType > _rng, int dim = 3);
      using Super::operator();
    };
  }
}
#endif
