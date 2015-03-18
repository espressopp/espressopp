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
#ifndef _INTEGRATOR_CAPFORCE_HPP
#define _INTEGRATOR_CAPFORCE_HPP

#include "types.hpp"
#include "logging.hpp"
#include "Extension.hpp"
#include "ParticleGroup.hpp"
#include "boost/signals2.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"

namespace espressopp {
  namespace integrator {

    /** CapForce */

    class CapForce : public Extension {

      public:

        CapForce(shared_ptr< System > _system, const Real3D& _CapForce);

        CapForce(shared_ptr< System > _system, real _AbsCapForce);

        CapForce(shared_ptr< System > _system, const Real3D& _CapForce, shared_ptr< ParticleGroup > _particleGroup);

        CapForce(shared_ptr< System > _system, real _AbsCapForce, shared_ptr< ParticleGroup > _particleGroup);

        void setCapForce(Real3D& _CapForce);

        void setAbsCapForce(real _AbsCapForce);

        Real3D& getCapForce();
        real getAbsCapForce();

        void setAdress(bool _adress);
        bool getAdress();
        
        void setParticleGroup(shared_ptr< ParticleGroup > _particleGroup);

        shared_ptr< ParticleGroup > getParticleGroup();

        void applyForceCappingToGroup();
        void applyForceCappingToAll();

        virtual ~CapForce() {};

        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:
        boost::signals2::connection _aftCalcF;
        shared_ptr< ParticleGroup > particleGroup;
        bool allParticles;
        bool absCapping;
        bool adress;
        Real3D capForce;
        real absCapForce;
        void connect();
        void disconnect();

        /** Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
