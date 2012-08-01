// ESPP_CLASS
#ifndef _ESUTIL_RNG_HPP
#define _ESUTIL_RNG_HPP
#include <boost/random.hpp>
#include "Real3D.hpp"
#include <vector>


#include "types.hpp"

namespace espresso {
  namespace esutil {
    /** A class that allows to choose between different RNGs, but that
	provides a uniform interface for these.
    */
    using namespace boost;

    class RNG {

    public:
      /** Init the RNG, use the given seed. */
      RNG(long _seed = 12345);

      /** Seed the RNG. */
      void seed(long _seed);
    
      /** returns a uniformly distributed random number between 0 and
	  1. */
      real operator()();

      /** returns a uniformly distributed integer random number in the
	  interval [0, N-1]. */
      int operator()(int N);

      /** returns a normal distributed random number, with mean of 0.0
	  and sigma of 1.0. If you need to generate many normal
	  idstributed random numbers with different mean or sigma,
	  use createNormalVariate() to create a variate generator object. */
      real normal();

      /** returns a random 3D vector that is uniformly distributed on a sphere. */
      Real3D uniformOnSphere();

      shared_ptr< RNGType > getBoostRNG();

      static void registerPython();

    private:
      typedef 
      variate_generator< RNGType&, normal_distribution< real > >
      NormalVariate;

      typedef 
      variate_generator< RNGType&, uniform_on_sphere< real, Real3D > > 
      UniformOnSphereVariate;
    
      shared_ptr< RNGType > boostRNG;

      NormalVariate normalVariate;

      UniformOnSphereVariate uniformOnSphereVariate;
    };
  }
}
#endif
