#ifndef _ESUTIL_UNIFORMONSPHERE_HPP
#define _ESUTIL_UNIFORMONSPHERE_HPP
#include <boost/random.hpp>
#include "types.hpp"
#include "Real3D.hpp"

namespace espresso {
  namespace esutil {
    using namespace boost;
    /** This class generates random vectors that are uniformly
	distributed on a sphere. */
    class UniformOnSphere
      : variate_generator< RNGType&, uniform_on_sphere< real, Real3D > > 
    {
      typedef uniform_on_sphere< real, Real3D > DistType;
      typedef variate_generator< RNGType&, DistType > Super;
      
      /// store the shared pointer to the RNG
      shared_ptr< RNG > rng;
      
    public:
      UniformOnSphere(shared_ptr< RNG > _rng);
      using Super::operator();
      static void registerPython();
    };
  }
}
#endif
