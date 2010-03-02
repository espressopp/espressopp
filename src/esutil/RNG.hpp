#ifndef _ESUTIL_RNG_HPP
#define _ESUTIL_RNG_HPP
#include <boost/random.hpp>
#include <vector>

#include "types.hpp"
#include "NormalVariate.hpp"
#include "UniformOnSphere.hpp"

namespace espresso {
  namespace esutil {
    /** A class that allows to choose between different RNGs, but that
	provides a uniform interface for these.
    */
    using namespace boost;

    class RNG {
      shared_ptr< RNGType > boostRNG;

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

      /** returns a variate generator object that generates random
	  variates with a normal distribution and the given parameters. */
      shared_ptr < NormalVariate > 
      createNormalVariate(const real mean = 0.0,
			  const real sigma = 1.0);

      /** returns a random 3D vector that is uniformly distributed on a sphere. */
      std::vector< real > uniformOnSphere();

      shared_ptr< UniformOnSphere >
      createUniformOnSphere(int dim = 3);

      static void registerPython();

    private:
      // the standard normal variate (mean 0, sigma 1)
      NormalVariate normalVariate;
      UniformOnSphere uniformOnSphereVariate;
    };
  }
}
#endif
