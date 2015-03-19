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
#ifndef _INTEGRATOR_LANGEVINBAROSTAT_HPP
#define _INTEGRATOR_LANGEVINBAROSTAT_HPP

#include "types.hpp"
#include "logging.hpp"

#include "Extension.hpp"
#include "VelocityVerlet.hpp"

#include "boost/signals2.hpp"

namespace espressopp {
  namespace integrator {

    /** Barostat. In a couple with Langevin thermostat it performs Langevin dynamics in 
     * a Hoover-style extended system. See
     * D. Quigley, M.I.J. Probert, J. Chem. Phys., 120, 2004, p. 11432 */
    class LangevinBarostat : public Extension {

      public:
        LangevinBarostat(shared_ptr< System >, shared_ptr< esutil::RNG >, real);

        void setGammaP(real);
        real getGammaP();

        void setPressure(real);
        real getPressure();
        
        void setMass(real);
        real getMass();

        void setMassByFrequency(real);
        
        virtual ~LangevinBarostat();

        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:
        boost::signals2::connection _runInit, _befIntP, _inIntP, _aftIntV, _aftCalcF;
        
        void frictionBarostat(class Particle&, real);

        // initial parameters
        real gammaP;                // friction coefficient for pressure control
        real mass;                  // fictitious mass
        real externalPressure;      // desired external pressure
        
        real desiredTemperature;    // desired temperature
        
        // system variable
        real momentum;              // momentum variable corresponding to the volume variable
        real momentum_mass;         // momentum variable divided by mass pe/w
        
        
        real pref3;  // prefactor, force term
        real pref4;  // prefactor, for the volume momentum
        real pref5;  // prefactor, for the volume momentum
        real pref6;  // prefactor, for the volume momentum

        shared_ptr< esutil::RNG > rng;  //!< random number generator used for friction term

        void initialize();    // initialize barostat prefactors
        
        void upd_Vp(); // it is for signals at first we modify volume then momentum
        void upd_pV(); // the other way around

        void updVolume();           // scale the volume according to the evolution equations
        void updVolumeMomentum();   // update local momentum which corresponds to the volume variable
        void updForces();               // update forces with an additional term
        void updDisplacement(real&);         // returns pe/W in order to update particle positions
        
        void connect();
        void disconnect();
        
        /** Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
