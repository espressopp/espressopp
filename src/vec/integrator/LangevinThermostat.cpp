/*
  Copyright (C) 2020-2021
      Max Planck Institute for Polymer Research & JGU Mainz
  Copyright (C) 2012,2013,2014,2015,2016
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

#include "LangevinThermostat.hpp"
#include "vec/Vectorization.hpp"
#include "vec/iterator/ParticleArrayIterator.hpp"

#include "python.hpp"
#include "types.hpp"
#include "System.hpp"
#include "SystemAccess.hpp"
#include "storage/Storage.hpp"
#include "esutil/RNG.hpp"

namespace espressopp { namespace vec {

  namespace integrator {

    LangevinThermostat::LangevinThermostat(std::shared_ptr<System> system)
      : Extension(system)
    {
      type = Extension::Thermostat;
      gamma  = 0.0;
      temperature = 0.0;
      adress = false;
      exclusions.clear();

      if (!getSystem()->rng) {
        throw std::runtime_error("system has no RNG");
      }
      rng = getSystem()->rng;

      LOG4ESPP_INFO(theLogger, "Langevin constructed");
    }

    void LangevinThermostat::setGamma(real _gamma)
    {
      gamma = _gamma;
    }

    real LangevinThermostat::getGamma()
    {
      return gamma;
    }

    void LangevinThermostat::setAdress(bool _adress)
    {
      if(_adress){
        throw std::runtime_error("LangevinThermostat::setAdress Not implemented for _adress=true");
      }
      adress = _adress;
    }

    bool LangevinThermostat::getAdress()
    {
      return adress;
    }

    void LangevinThermostat::setTemperature(real _temperature)
    {
      temperature = _temperature;
    }

    real LangevinThermostat::getTemperature()
    {
      return temperature;
    }

    LangevinThermostat::~LangevinThermostat()
    {
      disconnect();
    }

    void LangevinThermostat::disconnect()
    {
      _initialize.disconnect();
      _heatUp.disconnect();
      _coolDown.disconnect();
      _thermalize.disconnect();
      _thermalizeAdr.disconnect();
    }

    void LangevinThermostat::connect()
    {
      // connect to initialization inside run()
      _initialize = integrator->runInit.connect(
              std::bind(&LangevinThermostat::initialize, this));

      _heatUp = integrator->recalc1.connect(
              std::bind(&LangevinThermostat::heatUp, this));

      _coolDown = integrator->recalc2.connect(
              std::bind(&LangevinThermostat::coolDown, this));

      if (adress) {
          _thermalizeAdr = integrator->aftCalcF.connect(
              std::bind(&LangevinThermostat::thermalizeAdr, this));
      }
      else {
          _thermalize = integrator->aftCalcF.connect(
              std::bind(&LangevinThermostat::thermalize, this));
      }
    }


    void LangevinThermostat::thermalize()
    {
      LOG4ESPP_DEBUG(theLogger, "thermalize");

      if(!exclusions.empty()) throw std::runtime_error(
        "LangevinThermostat::thermalize exclusions not implemented");

      auto& particles = getSystem()->vectorization->particles;
      for(iterator::ParticleArrayIterator pit(particles, true); pit.isValid(); ++pit)
      {
        const real mass  = pit.mass();
        const real massf = sqrt(mass);

        const real ranval_x = (*rng)() - 0.5;
        const real ranval_y = (*rng)() - 0.5;
        const real ranval_z = (*rng)() - 0.5;

        pit.f_x() += pref1 * pit.v_x() * mass + pref2 * ranval_x * massf;
        pit.f_y() += pref1 * pit.v_y() * mass + pref2 * ranval_y * massf;
        pit.f_z() += pref1 * pit.v_z() * mass + pref2 * ranval_z * massf;
      }
    }

    // for AdResS
    void LangevinThermostat::thermalizeAdr()
    {
      throw std::runtime_error("LangevinThermostat::thermalizeAdr Function not implemented");

    #if 0
      LOG4ESPP_DEBUG(theLogger, "thermalize");

      System& system = getSystemRef();

      // thermalize AT particles
      ParticleList& adrATparticles = system.storage->getAdrATParticles();
      for (auto it = adrATparticles.begin(); it != adrATparticles.end(); it++) {
        if(exclusions.count((*it).id()) == 0)
        {
          frictionThermo(*it);
        }
      }
    #endif
    }

    void LangevinThermostat::initialize()
    {
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
    void LangevinThermostat::heatUp()
    {
      LOG4ESPP_INFO(theLogger, "heatUp");

      pref2buffer = pref2;
      pref2       *= sqrt(3.0);
    }

    /** Opposite to heatUp */

    void LangevinThermostat::coolDown()
    {
      LOG4ESPP_INFO(theLogger, "coolDown");

      pref2 = pref2buffer;
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void LangevinThermostat::registerPython()
    {
      using namespace espressopp::python;

      class_<LangevinThermostat, std::shared_ptr<LangevinThermostat>, bases<Extension> >
        ("vec_integrator_LangevinThermostat", init<std::shared_ptr<System>>())
        .def("connect", &LangevinThermostat::connect)
        .def("disconnect", &LangevinThermostat::disconnect)
        .def("addExclpid", &LangevinThermostat::addExclpid)
        .add_property("adress", &LangevinThermostat::getAdress, &LangevinThermostat::setAdress)
        .add_property("gamma", &LangevinThermostat::getGamma, &LangevinThermostat::setGamma)
        .add_property("temperature", &LangevinThermostat::getTemperature, &LangevinThermostat::setTemperature)
        ;
    }
  }
}}

