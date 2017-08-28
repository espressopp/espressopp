/*
  Copyright (C) 2017
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
#ifndef _INTEGRATOR_LANGEVINTHERMOSTATONRADIUS_HPP
#define _INTEGRATOR_LANGEVINTHERMOSTATONRADIUS_HPP

#include "types.hpp"
#include "logging.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"

#include "Extension.hpp"
#include "VelocityVerlet.hpp"


#include "boost/signals2.hpp"
#include "boost/unordered_set.hpp"

namespace espressopp {
    namespace integrator {
	
	/** Langevin thermostat */

	class LangevinThermostatOnRadius : public Extension {
	    
	public:
	    
	    LangevinThermostatOnRadius(shared_ptr<System> system, real _dampingmass);
	    virtual ~LangevinThermostatOnRadius();
	    
	    void setGamma(real gamma);
	    real getGamma();
	    
	    void setTemperature(real temperature);
	    real getTemperature();
	    
	    void initialize();
	    
	    /** update of forces to thermalize the system */
	    void thermalize();
	    
	    /** Add pid to exclusionlist */
	    void addExclpid(int pid) { exclusions.insert(pid); }
	    
	    /** very nasty: if we recalculate force when leaving/reentering the integrator,
		a(t) and a((t-dt)+dt) are NOT equal in the vv algorithm. The random
		numbers are drawn twice, resulting in a different variance of the random force.
		This is corrected by additional heat when restarting the integrator here.
		Currently only works for the Langevin thermostat, although probably also others
		are affected.
	    */
	    void heatUp();
	    
	    /** Opposite to heatUp */
	    void coolDown();
	    
	    /** Register this class so it can be used from Python. */
	    static void registerPython();
	    
	private:
	    
	    boost::signals2::connection _initialize, _heatUp, _coolDown,
		_thermalize;
	    
	    void frictionThermo(class Particle&);
	    
	    /** pid eclusion list */
	    std::set<longint> exclusions;
	    
	    void connect();
	    void disconnect();
	    
	    real gamma;        //!< friction coefficient
	    real temperature;  //!< desired user temperature
	    
	    real dampingmass; //!< radial damping mass
	    real massf; //!< square root of radial damping mass
	    
	    real pref1;  //!< prefactor, reduces complexity of thermalize
	    real pref2;  //!< prefactor, reduces complexity of thermalize
	    
	    real pref2buffer; //!< temporary to save value between heatUp/coolDown
	    
	    shared_ptr< esutil::RNG > rng;  //!< random number generator used for friction term
	    
	};
    }
}

#endif
