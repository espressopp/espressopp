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

/* Berendsen barostat */

#ifndef _INTEGRATOR_BERENDSENBAROSTAT_HPP
#define	_INTEGRATOR_BERENDSENBAROSTAT_HPP

#include "types.hpp"
#include "logging.hpp"

#include "analysis/Pressure.hpp"
#include "Extension.hpp"

#include "boost/signals2.hpp"
#include "Int3D.hpp"

namespace espressopp {
  
  using namespace analysis;

  namespace integrator {

    class BerendsenBarostat: public Extension {

      public:
        BerendsenBarostat(shared_ptr< System > system);
        
        void setFixed(Int3D);
        Int3D getFixed();
        void setTau(real);
        real getTau();
        void setPressure(real);
        real getPressure();

        ~BerendsenBarostat();

        void connect();
        void disconnect();
        
        /* Register in Python. */
        static void registerPython();

      private:
        boost::signals2::connection _runInit, _aftIntV;
        
        real tau;   // time constant
        real P0;    // external pressure
        
        real pref;  // prefactor for the pressure calc
        
        Int3D fixed; // fixed directions. If (0,1,1) then Lx=const.
                     // By default (1,1,1). Can not be (0,0,0)
        real exponent; // precalculated exponent
        
        void initialize();

        /* rescale the system size and coord. of particles */
        void barostat();

        /* Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif	/* _INTEGRATOR_BERENDSENBAROSTAT_HPP */

