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

#include "python.hpp"
#include "LangevinThermostatHybrid.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"

namespace espressopp {

  namespace integrator {

    using namespace espressopp::iterator;


    LangevinThermostatHybrid::LangevinThermostatHybrid(shared_ptr<System> system,shared_ptr<FixedTupleListAdress> _fixedtupleList)
    :Extension(system),fixedtupleList(_fixedtupleList) {

      type = Extension::Thermostat;

      gamma  = 0.0;
      gammahy  = 0.0;
      gammacg  = 0.0;
      temperature = 0.0;
      
      if (!system->rng) {
        throw std::runtime_error("system has no RNG");
      }

      rng = system->rng;

      LOG4ESPP_INFO(theLogger, "Langevin constructed");


    }

    void LangevinThermostatHybrid::setGamma(real _gamma)
    {
      gamma = _gamma;
    }

    real LangevinThermostatHybrid::getGamma()
    {
      return gamma;
    }

    void LangevinThermostatHybrid::setGammaHybrid(real _gammahy)
    {
      gammahy = _gammahy;
    }

    real LangevinThermostatHybrid::getGammaHybrid()
    {
      return gammahy;
    }

    void LangevinThermostatHybrid::setGammaCG(real _gammacg)
    {
      gammacg = _gammacg;
    }

    real LangevinThermostatHybrid::getGammaCG()
    {
      return gammacg;
    }

    void LangevinThermostatHybrid::setTemperature(real _temperature)
    {
      temperature = _temperature;
    }

    real LangevinThermostatHybrid::getTemperature()
    {
      return temperature;
    }

    LangevinThermostatHybrid::~LangevinThermostatHybrid() {
        disconnect();
    }


    void LangevinThermostatHybrid::disconnect() {

        _initialize.disconnect();
        _heatUp.disconnect();
        _coolDown.disconnect();
        _thermalizeAdr.disconnect();

    }

    void LangevinThermostatHybrid::connect() {

        // connect to initialization inside run()
        _initialize = integrator->runInit.connect(
                boost::bind(&LangevinThermostatHybrid::initialize, this));

        _heatUp = integrator->recalc1.connect(
                boost::bind(&LangevinThermostatHybrid::heatUp, this));

        _coolDown = integrator->recalc2.connect(
                boost::bind(&LangevinThermostatHybrid::coolDown, this));

        _thermalizeAdr = integrator->aftCalcF.connect(
                boost::bind(&LangevinThermostatHybrid::thermalizeAdr, this));
    }

    void LangevinThermostatHybrid::thermalizeAdr()
    {
      LOG4ESPP_DEBUG(theLogger, "thermalizeAdr");

      System& system = getSystemRef();

      // thermalize CG particles
      /*CellList cells = system.storage->getRealCells();
      for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
        frictionThermo(*cit);
      }*/

      // In force-Adress, it doesn't make sense to thermalize both CG and AT particles, since CG particles get velocities of AT particles anyway.
      // TODO May not be compatible with H-Adress
      
      // thermalize AT particles
      CellList cells = system.storage->getRealCells();
      for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
        Particle &vp = *cit;
        real weight = vp.lambda();
        FixedTupleListAdress::iterator it3;
        it3 = fixedtupleList->find(&vp);
        if (it3 != fixedtupleList->end()) {

            std::vector<Particle*> atList;
            atList = it3->second;
            for (std::vector<Particle*>::iterator it2 = atList.begin();
                                 it2 != atList.end(); ++it2) {
                Particle &at = **it2;
                frictionThermo(at,weight);
            }
        }
      }
      
    }

    void LangevinThermostatHybrid::frictionThermo(Particle& p,real weight)
    {
      real massf = sqrt(p.mass());

      // get a random value for each vector component

      Real3D ranval((*rng)() - 0.5, (*rng)() - 0.5, (*rng)() - 0.5);

      if (weight<1.0 && weight>0.0) { //hybrid region
        p.force() += pref1hy * p.velocity() * p.mass() +
                   pref2hy * ranval * massf;
      } else if (weight==1.0) {
        p.force() += pref1 * p.velocity() * p.mass() +
                   pref2 * ranval * massf;
      } else {
        p.force() += pref1cg * p.velocity() * p.mass() +
                   pref2cg * ranval * massf;
      }

      LOG4ESPP_TRACE(theLogger, "new force of p = " << p.force());
    }

    void LangevinThermostatHybrid::initialize()
    { // calculate the prefactors

        real timestep = integrator->getTimeStep();

      LOG4ESPP_INFO(theLogger, "init, timestep = " << timestep <<
		    ", gamma = " << gamma << 
		    ", gammahy = " << gammahy << 
		    ", gammacg = " << gammacg << 
		    ", temperature = " << temperature);

      pref1 = -gamma;
      pref2 = sqrt(24.0 * temperature * gamma / timestep);

      pref1hy = -gammahy;
      pref2hy = sqrt(24.0 * temperature * gammahy / timestep);

      pref1cg = -gammacg;
      pref2cg = sqrt(24.0 * temperature * gammacg / timestep);

    }

    /** very nasty: if we recalculate force when leaving/reentering the integrator,
	a(t) and a((t-dt)+dt) are NOT equal in the vv algorithm. The random
	numbers are drawn twice, resulting in a different variance of the random force.
	This is corrected by additional heat when restarting the integrator here.
	Currently only works for the Langevin thermostat, although probably also others
	are affected.
    */

    void LangevinThermostatHybrid::heatUp()
    {
      LOG4ESPP_INFO(theLogger, "heatUp");

      pref2buffer = pref2;
      pref2       *= sqrt(3.0);
      pref2bufferhy = pref2hy;
      pref2hy       *= sqrt(3.0);
      pref2buffercg = pref2cg;
      pref2cg       *= sqrt(3.0);
    }

    /** Opposite to heatUp */

    void LangevinThermostatHybrid::coolDown()
    {
      LOG4ESPP_INFO(theLogger, "coolDown");

      pref2 = pref2buffer;
      pref2hy = pref2bufferhy;
      pref2cg = pref2buffercg;

    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void LangevinThermostatHybrid::registerPython() {


      using namespace espressopp::python;


      class_<LangevinThermostatHybrid, shared_ptr<LangevinThermostatHybrid>, bases<Extension> >
        ("integrator_LangevinThermostatHybrid", init<shared_ptr<System>,shared_ptr<FixedTupleListAdress> >())
        .def("connect", &LangevinThermostatHybrid::connect)
        .def("disconnect", &LangevinThermostatHybrid::disconnect)
        .add_property("gamma", &LangevinThermostatHybrid::getGamma, &LangevinThermostatHybrid::setGamma)
        .add_property("gammahy", &LangevinThermostatHybrid::getGammaHybrid, &LangevinThermostatHybrid::setGammaHybrid)
        .add_property("gammacg", &LangevinThermostatHybrid::getGammaCG, &LangevinThermostatHybrid::setGammaCG)
        .add_property("temperature", &LangevinThermostatHybrid::getTemperature, &LangevinThermostatHybrid::setTemperature)
        ;


    }

  }
}

