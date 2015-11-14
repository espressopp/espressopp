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
#ifndef _INTEGRATOR_ISOKINETIC_HPP
#define _INTEGRATOR_ISOKINETIC_HPP

#include "types.hpp"
#include "logging.hpp"

#include "Extension.hpp"
#include "boost/signals2.hpp"

namespace espressopp {
  namespace integrator {
    /** Langevin */

    class Isokinetic : public Extension{

      public:
        Isokinetic(shared_ptr< System > system);

        void setTemperature(real temperature);

        real getTemperature();

        void setCoupling(int coupling);

        int getCoupling();

        ~Isokinetic();

        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:
        boost::signals2::connection _aftIntV;
        
        real temperature;  //!< desired user temperature
        int coupling; // how often to couple to the thermostat
        int couplecount;

        // not yet needed
        // shared_ptr< esutil::RNG > rng;  //!< random number generator used for friction term

        void rescaleVelocities();

        void connect();
        void disconnect();
        
        /** Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
