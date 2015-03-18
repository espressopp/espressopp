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
#ifndef _INTERACTION_ZERO_HPP
#define _INTERACTION_ZERO_HPP

#include "Potential.hpp"

namespace espressopp {
  namespace interaction {
    /** This class provides methods for a zero potential
     * no interactions between particles, mainly used for debugging and testing
    */
    class Zero : public PotentialTemplate< Zero > {

    public:
      static void registerPython();

      Zero() {} ;

      real _computeEnergySqrRaw(real distSqr) const {
        return 0;
      }

      bool _computeForceRaw(Real3D& force,
                            const Real3D& dist,
                            real distSqr) const {
        force = Real3D(0,0,0);
        return true;
      }
    };
    // provide pickle support
    struct Zero_pickle : boost::python::pickle_suite
    {
      static
      boost::python::tuple
      getinitargs(Zero const& pot)
      {
          return boost::python::make_tuple();
      }
    };

  }
}

#endif
