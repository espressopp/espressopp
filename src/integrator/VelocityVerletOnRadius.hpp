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
#ifndef _INTEGRATOR_VELOCITYVERLETONRADIUS_HPP
#define _INTEGRATOR_VELOCITYVERLETONRADIUS_HPP

#include "types.hpp"
#include "logging.hpp"
#include "Extension.hpp"
#include "ParticleGroup.hpp"
#include "boost/signals2.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
namespace espressopp {
  namespace integrator {

    class VelocityVerletOnRadius : public Extension {
      public:
        VelocityVerletOnRadius(shared_ptr< System > _system, real _radialDampingMass);
        virtual ~VelocityVerletOnRadius() {};

        void setRadialDampingMass(real _radialDampingMass) {
        	radialDampingMass = _radialDampingMass;
        }

        real getRadialDampingMass() {
        	return radialDampingMass;
        }

        /** Register this class so it can be used from Python. */
        static void registerPython();
      private:
        boost::signals2::connection _aftIntP, _aftIntV, _aftInitF;
        void connect();
        void disconnect();
        void integrate1();
        void integrate2();
        void initForces();
        real radialDampingMass;
        /** Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };

  }
}
#endif
