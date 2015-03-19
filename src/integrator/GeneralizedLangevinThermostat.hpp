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
#ifndef _INTEGRATOR_GENERALIZEDLANGEVINTHERMOSTAT_HPP
#define _INTEGRATOR_GENERALIZEDLANGEVINTHERMOSTAT_HPP

#include "types.hpp"
#include "logging.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"
#include "interaction/Interpolation.hpp"
#include <map>

#include "Extension.hpp"
#include "VelocityVerlet.hpp"


#include "boost/signals2.hpp"


namespace espressopp {
  namespace integrator {

    /** Langevin thermostat */

    class GeneralizedLangevinThermostat : public Extension {

      public:

        GeneralizedLangevinThermostat(shared_ptr<System> system);
        virtual ~GeneralizedLangevinThermostat();

                /** Setter for the filename, will read in the table. */
        void addCoeffs(int itype, const char* _filename, int type);
        const char* getFilename() const { return filename.c_str(); }
        
        void integrate();
        void friction();
        //void setGamma(real gamma);
        //real getGamma();
        //void setTemperature(real temperature);
        //real getTemperature();
        //void setAdress(bool _adress);
        //bool getAdress();
        
        //void initialize();
        /** update of forces to thermalize the system */
        //void thermalizeAdr(); // same as above, for AdResS

        /** very nasty: if we recalculate force when leaving/reentering the integrator,
            a(t) and a((t-dt)+dt) are NOT equal in the vv algorithm. The random
            numbers are drawn twice, resulting in a different variance of the random force.
            This is corrected by additional heat when restarting the integrator here.
            Currently only works for the Langevin thermostat, although probably also others
            are affected.
        */
        //void heatUp();

        /** Opposite to heatUp */
        //void coolDown();

        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:

        boost::signals2::connection _integrate, _friction;

        //void frictionThermo(class Particle&);

        // this connects thermalizeAdr
        //void enableAdress();
        //bool adress;

        void connect();
        void disconnect();

        std::string filename;
        typedef shared_ptr <interaction::Interpolation> Table;
        std::map<int, Table> coeffs; // map type to force

        //static LOG4ESPP_DECL_LOGGER(theLogger);
        //real gamma;        //!< friction coefficient
        //real temperature;  //!< desired user temperature

        //real pref1;  //!< prefactor, reduces complexity of thermalize
        //real pref2;  //!< prefactor, reduces complexity of thermalize

        //real pref2buffer; //!< temporary to save value between heatUp/coolDown

        //shared_ptr< esutil::RNG > rng;  //!< random number generator used for friction term

    };
  }
}

#endif
