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

#include "python.hpp"
#include "LangevinThermostatOnRadius.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"

namespace espressopp {

    namespace integrator {
	
	using namespace espressopp::iterator;
	
	
	LangevinThermostatOnRadius::LangevinThermostatOnRadius(shared_ptr<System> system, real _dampingmass)
	    :Extension(system) {
	    
	    type = Extension::Thermostat;
	    
	    gamma  = 0.0;
	    temperature = 0.0;
	    
	    dampingmass = _dampingmass;
	    massf = sqrt(dampingmass);
	    
	    exclusions.clear();
	    
	    if (!system->rng) {
		throw std::runtime_error("system has no RNG");
	    }
	    
	    rng = system->rng;
	    
	    LOG4ESPP_INFO(theLogger, "Langevin constructed");
	    
	    
	}
	
	void LangevinThermostatOnRadius::setGamma(real _gamma) {
	    gamma = _gamma;
	}
	
	real LangevinThermostatOnRadius::getGamma() {
	    return gamma;
	}
	
	void LangevinThermostatOnRadius::setTemperature(real _temperature) {
	    temperature = _temperature;
	}
	
	real LangevinThermostatOnRadius::getTemperature() {
	    return temperature;
	}
	
	LangevinThermostatOnRadius::~LangevinThermostatOnRadius() {
	    disconnect();
	}
	
	
	void LangevinThermostatOnRadius::disconnect() {
	    
	    _initialize.disconnect();
	    _heatUp.disconnect();
	    _coolDown.disconnect();
	    _thermalize.disconnect();
	    
	}
	
	void LangevinThermostatOnRadius::connect() {
	    
	    // connect to initialization inside run()
	    _initialize = integrator->runInit.connect(
                boost::bind(&LangevinThermostatOnRadius::initialize, this));
	    
	    _heatUp = integrator->recalc1.connect(
                boost::bind(&LangevinThermostatOnRadius::heatUp, this));
	    
	    _coolDown = integrator->recalc2.connect(
                boost::bind(&LangevinThermostatOnRadius::coolDown, this));
	    
	    _thermalize = integrator->aftCalcF.connect(
                boost::bind(&LangevinThermostatOnRadius::thermalize, this));
	}
	
	
	void LangevinThermostatOnRadius::thermalize() {
	    LOG4ESPP_DEBUG(theLogger, "thermalize");
	    
	    System& system = getSystemRef();
	    
	    CellList cells = system.storage->getRealCells();
	    
	    for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
		
		if(exclusions.count(cit->id()) == 0) {
		    frictionThermo(*cit);
		}
		
	    }
	}
	
	void LangevinThermostatOnRadius::frictionThermo(Particle& p) {
	    
	    // get a random value for a radius
	    
	    real ranval = (*rng)() - 0.5;
	    
	    p.fradius() += pref1 * p.vradius() * dampingmass +
		pref2 * ranval * massf;
	    
	    LOG4ESPP_TRACE(theLogger, "new force of p = " << p.force());
	}
	
	void LangevinThermostatOnRadius::initialize() { 

            // calculate the prefactors
	    
	    real timestep = integrator->getTimeStep();
	    
	    LOG4ESPP_INFO(theLogger, "init, timestep = " << timestep <<
			  ", gamma = " << gamma << 
			  ", temperature = " << temperature);
	    
	    pref1 = -gamma;
	    pref2 = sqrt(24.0 * temperature * gamma / timestep);
	    
	}
	
	/** very nasty: if we recalculate force when leaving/reentering the integrator,
	    a(t) and a((t-dt)+dt) are NOT equal in the vv algorithm. The random
	    numbers are drawn twice, resulting in a different variance of the random force.
	    This is corrected by additional heat when restarting the integrator here.
	    Currently only works for the Langevin thermostat, although probably also others
	    are affected.
	*/
	
	void LangevinThermostatOnRadius::heatUp() {
	    LOG4ESPP_INFO(theLogger, "heatUp");
	    
	    pref2buffer = pref2;
	    pref2       *= sqrt(3.0);
	}
	
	/** Opposite to heatUp */
	
	void LangevinThermostatOnRadius::coolDown() {
	    LOG4ESPP_INFO(theLogger, "coolDown");
	    
	    pref2 = pref2buffer;
	}
	
	/****************************************************
	 ** REGISTRATION WITH PYTHON
	 ****************************************************/
	
	void LangevinThermostatOnRadius::registerPython() {
	    
	    
	    using namespace espressopp::python;
	    
	    
	    class_<LangevinThermostatOnRadius, shared_ptr<LangevinThermostatOnRadius>, bases<Extension> >
		("integrator_LangevinThermostatOnRadius", init<shared_ptr<System>, real >())
		.def("connect", &LangevinThermostatOnRadius::connect)
		.def("disconnect", &LangevinThermostatOnRadius::disconnect)
		.def("addExclpid", &LangevinThermostatOnRadius::addExclpid)
		.add_property("gamma", &LangevinThermostatOnRadius::getGamma, &LangevinThermostatOnRadius::setGamma)
		.add_property("temperature", &LangevinThermostatOnRadius::getTemperature, &LangevinThermostatOnRadius::setTemperature)
		;
	    
	    
	}
	
    }
}
