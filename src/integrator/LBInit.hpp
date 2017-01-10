/*
  Copyright (C) 2012-2016
      Max Planck Institute for Polymer Research
  Copyright (C) 2008-2011
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
#ifndef _INTEGRATOR_LBINIT_HPP
#define _INTEGRATOR_LBINIT_HPP

#include "LatticeBoltzmann.hpp"

namespace espressopp {
  namespace integrator {
    /** Abstract base class for arbitrary Init for LB simulations. */
    class LBInit {
    public:
      /* Constructor for the class */
      LBInit(shared_ptr< System > _system,
                         shared_ptr< LatticeBoltzmann > _latticeboltzmann) {
                            latticeboltzmann = _latticeboltzmann;
      }
      /* Destructor for the class */
      virtual ~LBInit () {}

      /* HANDLING INITIAL DENSITIES AND VELOCITIES */
      virtual void createDenVel (real _rho0, Real3D _u0) = 0;

      /* HANDLING EXTERNAL FORCES */
      virtual void setForce (Real3D _force) = 0;
      virtual void addForce (Real3D _force) = 0;

      static void registerPython();

    protected:
      shared_ptr<LatticeBoltzmann> latticeboltzmann;
      real rho0;
      Real3D u0;
    };
  }
}

#endif
