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
#include "LangevinThermostat1D.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"

namespace espressopp {

  namespace integrator {

    using namespace espressopp::iterator;


    LangevinThermostat1D::LangevinThermostat1D(shared_ptr<System> system)
    :Extension(system) {

      type = Extension::Thermostat;

      gamma       = 0.0;
      temperature = 0.0;
      direction   = 0; // default is x direction

      if (!system->rng) {
        throw std::runtime_error("system has no RNG");
      }

      rng = system->rng;

      LOG4ESPP_INFO(theLogger, "Langevin constructed");


    }

    void LangevinThermostat1D::setGamma(real _gamma)
    {
      gamma = _gamma;
    }

    real LangevinThermostat1D::getGamma()
    {
      return gamma;
    }

    void LangevinThermostat1D::setAdress(bool _adress){
        adress = _adress;
    }

    bool LangevinThermostat1D::getAdress(){
        return adress;
    }

    void LangevinThermostat1D::setTemperature(real _temperature)
    {
      temperature = _temperature;
    }

    real LangevinThermostat1D::getTemperature()
    {
      return temperature;
    }

    void LangevinThermostat1D::setDirection(int _direction)
    {
      direction = _direction;
    }

    int LangevinThermostat1D::getDirection()
    {
      return direction;
    }

    LangevinThermostat1D::~LangevinThermostat1D() {
        disconnect();
    }


    void LangevinThermostat1D::disconnect() {

        _initialize.disconnect();
        _heatUp.disconnect();
        _coolDown.disconnect();
        _thermalize.disconnect();
        _thermalizeAdr.disconnect();

    }

    void LangevinThermostat1D::connect() {

        // connect to initialization inside run()
        _initialize = integrator->runInit.connect(
                boost::bind(&LangevinThermostat1D::initialize, this));

        _heatUp = integrator->recalc1.connect(
                boost::bind(&LangevinThermostat1D::heatUp, this));

        _coolDown = integrator->recalc2.connect(
                boost::bind(&LangevinThermostat1D::coolDown, this));

        if (adress) {
            _thermalizeAdr = integrator->aftCalcF.connect(
                boost::bind(&LangevinThermostat1D::thermalizeAdr, this));
        }
        else {
            _thermalize = integrator->aftCalcF.connect(
                boost::bind(&LangevinThermostat1D::thermalize, this));
        }
    }


    void LangevinThermostat1D::thermalize()
    {
      LOG4ESPP_DEBUG(theLogger, "thermalize");

      System& system = getSystemRef();
      
      CellList cells = system.storage->getRealCells();

      for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
        frictionThermo(*cit);
      }
    }

    // for AdResS
    void LangevinThermostat1D::thermalizeAdr()
    {
      LOG4ESPP_DEBUG(theLogger, "thermalize");

      System& system = getSystemRef();

      // thermalize CG particles
      CellList cells = system.storage->getRealCells();
      for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
        frictionThermo(*cit);
      }

      // TODO: It doesn't make that much sense to thermalize both CG and AT particles, since CG particles get velocities of AT particles anyway.
      
      // thermalize AT particles
      ParticleList& adrATparticles = system.storage->getAdrATParticles();
      for (std::vector<Particle>::iterator it = adrATparticles.begin();
              it != adrATparticles.end(); it++) {
        frictionThermo(*it);
      }
      
    }

    void LangevinThermostat1D::frictionThermo(Particle& p)
    {
      real massf = sqrt(p.mass());
      real ranval = (*rng)() - 0.5;

      p.force()[direction] += pref1 * p.velocity()[direction] * p.mass() +
                   pref2 * ranval * massf;

      LOG4ESPP_TRACE(theLogger, "new force of p = " << p.force());
    }

    void LangevinThermostat1D::initialize()
    { // calculate the prefactors

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

    void LangevinThermostat1D::heatUp()
    {
      LOG4ESPP_INFO(theLogger, "heatUp");

      pref2buffer = pref2;
      pref2       *= sqrt(3.0);
    }

    /** Opposite to heatUp */

    void LangevinThermostat1D::coolDown()
    {
      LOG4ESPP_INFO(theLogger, "coolDown");

      pref2 = pref2buffer;
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void LangevinThermostat1D::registerPython() {


      using namespace espressopp::python;


      class_<LangevinThermostat1D, shared_ptr<LangevinThermostat1D>, bases<Extension> >
        ("integrator_LangevinThermostat1D", init<shared_ptr<System> >())
        .def("connect", &LangevinThermostat1D::connect)
        .def("disconnect", &LangevinThermostat1D::disconnect)
        .add_property("adress", &LangevinThermostat1D::getAdress, &LangevinThermostat1D::setAdress)
        .add_property("gamma", &LangevinThermostat1D::getGamma, &LangevinThermostat1D::setGamma)
        .add_property("direction", &LangevinThermostat1D::getDirection, &LangevinThermostat1D::setDirection)
        .add_property("temperature", &LangevinThermostat1D::getTemperature, &LangevinThermostat1D::setTemperature)
        ;


    }

  }
}

