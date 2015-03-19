/*
  Copyright (C) 2012,2013
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
#ifndef _ESUTIL_GAMMAVARIATE_HPP
#define _ESUTIL_GAMMAVARIATE_HPP
#include <boost/random.hpp>
#include "types.hpp"

namespace espressopp {
  namespace esutil {
    using namespace boost;
      /** This class generates gamma distributed random
	  variates. The class also keeps a shared pointer to the boost
	  RNG object, so that it is not destroyed if the espressopp RNG
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
