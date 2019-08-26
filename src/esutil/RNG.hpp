/*
  Copyright (C) 2012,2013,2019
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI

  This file is part of ESPResSo++.

  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// ESPP_CLASS
#ifndef _ESUTIL_RNG_HPP
#define _ESUTIL_RNG_HPP
#include <boost/random.hpp>
#include "Real3D.hpp"
#include <vector>


#include "types.hpp"

namespace espressopp {
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

      /** Gets RNG seed. */
      long get_seed();

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

      /** returns a gamma distributed random number with shape parameter
       \alpha and scale parameter 1.  */
      real gamma(unsigned int alpha=1);
      real gammaOf1(){return gamma(1);}; //this would look nicer with an optional argument, but it did not work (for boost)

      /** returns a random 3D vector that is uniformly distributed on a sphere. */
      Real3D uniformOnSphere();

      shared_ptr< RNGType > getBoostRNG();

      void saveState(long long);

      void loadState(const char*);

      static void registerPython();

    private:
      long seed_;

      typedef
      variate_generator< RNGType&, normal_distribution< real > >
      NormalVariate;

      //      typedef
      //      variate_generator< RNGType&, gamma_distribution< real > >
      //      GammaVariate;

      typedef
      variate_generator< RNGType&, uniform_on_sphere< real, Real3D > >
      UniformOnSphereVariate;

      shared_ptr< RNGType > boostRNG;

      NormalVariate normalVariate;

      //GammaVariate gammaVariate;

      UniformOnSphereVariate uniformOnSphereVariate;
    };
  }
}
#endif
