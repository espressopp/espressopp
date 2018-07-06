/*
  Copyright (C) 2018
      Max Planck Institute for Polymer Research

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
#ifndef _INTEGRATOR_VELOCITYVERLETRESPA_HPP
#define _INTEGRATOR_VELOCITYVERLETRESPA_HPP

#include "types.hpp"
#include "MDIntegrator.hpp"
#include <boost/signals2.hpp>

namespace espressopp {
  namespace integrator {

    /** Velocity Verlet Integrator */
    class VelocityVerletRESPA : public MDIntegrator {

      public:

        VelocityVerletRESPA(shared_ptr<class espressopp::System> system);

        virtual ~VelocityVerletRESPA();

        void run(int nsteps);

        /** Setter routine for multistep. */
        void setmultistep(int _multistep);
        /** Getter routine for multistep. */
        int getmultistep() {
          return multistep;
        }

        /** Setter routine for timestep. */
        void setTimeStep(real _dt);

        /** Register this class so it can be used from Python. */
        static void registerPython();

      protected:
        bool resortFlag;
        real maxDist;
        real maxCut;
        int multistep;
        real dtlong;

        real integrate1();
        void integrate2(bool slow);
        void initForces();
        void updateForces(bool slow);
        void calcForces(bool slow);
    };
  }
}

#endif
