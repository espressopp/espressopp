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

#ifndef _LBINIT_POPWAVE_HPP
#define _LBINIT_POPWAVE_HPP

#include "LBInit.hpp"

namespace espressopp {
  namespace integrator {
    class LBInitPopWave : public LBInit {
      public:
      LBInitPopWave(shared_ptr<System> _system,
                    shared_ptr< LatticeBoltzmann > _latticeboltzmann);

      void createDenVel (real _rho0, Real3D _u0);

      void setForce (Real3D _force);
      void addForce (Real3D _force);

      static void registerPython();
      private:
    };
  }
}

#endif
