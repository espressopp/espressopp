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

/* Berendsen thermostat */

#ifndef _INTEGRATOR_BERENDSENTHERMOSTAT_HPP
#define	_INTEGRATOR_BERENDSENTHERMOSTAT_HPP

#include "types.hpp"
#include "logging.hpp"

#include "analysis/Temperature.hpp"
#include "Extension.hpp"
#include "VelocityVerlet.hpp"

#include "boost/signals2.hpp"

namespace espressopp {
  
  using namespace analysis;

  namespace integrator {

    class BerendsenThermostat: public Extension {

      public:
        BerendsenThermostat(shared_ptr< System > system);
        
        void setTau(real);
        real getTau();
        void setTemperature(real);
        real getTemperature();

        ~BerendsenThermostat();
        
        /* Register in Python. */
        static void registerPython();

      private:
        boost::signals2::connection _runInit, _aftIntV;

        real tau;   // time constant 1/(2*gamma), where gamma = friction constant
        real T0;    // external temperature
        
        real pref;  // prefactor for the temperature calculations
        
        void initialize();

        // velocities scaling per time step
        void thermostat();
        void scaleVelocity(Particle&, real);

        void connect();
        void disconnect();
        
        /* Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif	/* _INTEGRATOR_BERENDSENTHERMOSTAT_HPP */

