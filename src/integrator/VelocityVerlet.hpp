/*
  Copyright (C) 2017
      Jakub Krajniak (jkrajniak at gmail.com)
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
#ifndef _INTEGRATOR_VELOCITYVERLET_HPP
#define _INTEGRATOR_VELOCITYVERLET_HPP

#include "types.hpp"
#include "MDIntegrator.hpp"
#include "esutil/Timer.hpp"
#include <boost/signals2.hpp>

namespace espressopp {
  namespace integrator {

    /** Velocity Verlet Integrator */
    class VelocityVerlet : public MDIntegrator {

      public:

        VelocityVerlet(shared_ptr<class espressopp::System> system);

        virtual ~VelocityVerlet();

        void run(int nsteps);
        
        /** Load timings in array to export to Python as a tuple. */
        void loadTimers(real t[10]);

        void resetTimers();

        /** Register this class so it can be used from Python. */
        static void registerPython();

      protected:
        bool resortFlag;  //!< true implies need for resort of particles
        real maxDist;

        real maxCut;

        /** Method updates particle positions and velocities.
            \return maximal square distance a particle has moved.
        */
        real integrate1();

        void integrate2();

        void initForces();

        void updateForces();

        void calcForces();

        void printPositions(bool withGhost);

        void printForces(bool withGhost);

        void setUp();   //!< set up for a new run

        void printTimers();

        esutil::WallTimer timeIntegrate;  //!< used for timing

        // variables that keep time information about different phases
        real timeRun;
        real timeLost;
        real timeForce;
        real timeForceComp[10000];
        real timeComm1;
        real timeComm2;
        real timeInt1;
        real timeInt2;
        real timeResort;

        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
