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

#include "python.hpp"
#include "types.hpp"
#include "System.hpp"
#include "SystemAccess.hpp"
#include "storage/Storage.hpp"
#include "esutil/RNG.hpp"

namespace espressopp { namespace vec {

  namespace integrator {

    LangevinThermostat::LangevinThermostat(shared_ptr<Vectorization> vectorization)
      : Extension(vectorization)
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
              boost::bind(&LangevinThermostat::initialize, this));

      _heatUp = integrator->recalc1.connect(
              boost::bind(&LangevinThermostat::heatUp, this));

      _coolDown = integrator->recalc2.connect(
              boost::bind(&LangevinThermostat::coolDown, this));

      if (adress) {
          _thermalizeAdr = integrator->aftCalcF.connect(
              boost::bind(&LangevinThermostat::thermalizeAdr, this));
      }
      else {
          _thermalize = integrator->aftCalcF.connect(
              boost::bind(&LangevinThermostat::thermalize, this));
      }
    }


    void LangevinThermostat::thermalize()
    {
      LOG4ESPP_DEBUG(theLogger, "thermalize");

      if(!exclusions.empty()) throw std::runtime_error(
        "LangevinThermostat::thermalize exclusions not implemented");

      {
        auto& particles                    = vectorization->particles;
        const auto& realCells              = particles.realCells();
        const size_t* __restrict cellRange = particles.cellRange().data();
        const size_t* __restrict sizes     = particles.sizes().data();

        if(particles.mode_aos())
        {
          for(size_t ircell=0; ircell<realCells.size(); ircell++)
          {
            const size_t rcell = particles.realCells()[ircell];
            const size_t start = cellRange[rcell];
            const size_t size  = sizes[rcell];

            using espressopp::vec::Real4D;
            const Real4D* __restrict v = &(particles.velocity[start]);
            const real* __restrict mass = &(particles.mass[start]);

            Real4D* __restrict f = &(particles.force[start]);

            /// NOTE: Not vectorizable due to rng function call
            // #pragma vector always
            // #pragma vector aligned
            // #pragma ivdep
            for(size_t ip=0; ip<size; ip++)
            {
              const real massf = sqrt(mass[ip]);

              const real ranval_x = (*rng)() - 0.5;
              const real ranval_y = (*rng)() - 0.5;
              const real ranval_z = (*rng)() - 0.5;

              f[ip].x += pref1 * v[ip].x * mass[ip] + pref2 * ranval_x * massf;
              f[ip].y += pref1 * v[ip].y * mass[ip] + pref2 * ranval_y * massf;
              f[ip].z += pref1 * v[ip].z * mass[ip] + pref2 * ranval_z * massf;
            }
          }
        }
        else
        {
          for(size_t ircell=0; ircell<realCells.size(); ircell++)
          {
            const size_t rcell = particles.realCells()[ircell];
            const size_t start = cellRange[rcell];
            const size_t size  = sizes[rcell];

            const real* __restrict v_x = &(particles.v_x[start]);
            const real* __restrict v_y = &(particles.v_y[start]);
            const real* __restrict v_z = &(particles.v_z[start]);
            const real* __restrict mass = &(particles.mass[start]);

            real* __restrict f_x = &(particles.f_x[start]);
            real* __restrict f_y = &(particles.f_y[start]);
            real* __restrict f_z = &(particles.f_z[start]);

            /// NOTE: Not vectorizable due to rng function call
            // #pragma vector always
            // #pragma vector aligned
            // #pragma ivdep
            for(size_t ip=0; ip<size; ip++)
            {
              const real massf = sqrt(mass[ip]);

              const real ranval_x = (*rng)() - 0.5;
              const real ranval_y = (*rng)() - 0.5;
              const real ranval_z = (*rng)() - 0.5;

              f_x[ip] += pref1 * v_x[ip] * mass[ip] + pref2 * ranval_x * massf;
              f_y[ip] += pref1 * v_y[ip] * mass[ip] + pref2 * ranval_y * massf;
              f_z[ip] += pref1 * v_z[ip] * mass[ip] + pref2 * ranval_z * massf;
            }
          }
        }
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
      for (std::vector<Particle>::iterator it = adrATparticles.begin();
              it != adrATparticles.end(); it++) {

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

      class_<LangevinThermostat, shared_ptr<LangevinThermostat>, bases<Extension> >
        ("vec_integrator_LangevinThermostat", init<shared_ptr<Vectorization>>())
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

