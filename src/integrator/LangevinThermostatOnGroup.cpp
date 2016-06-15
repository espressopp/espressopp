/*
  Copyright (C) 2016
      Jakub Krajniak (c) (jkrjaniak at gmail.com)

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
#include "LangevinThermostatOnGroup.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"

namespace espressopp {
namespace integrator {

using namespace espressopp::iterator;

LangevinThermostatOnGroup::LangevinThermostatOnGroup(
    shared_ptr<System> system, shared_ptr<ParticleGroup> _pg)
        : Extension(system), particle_group(_pg) {
  type = Extension::Thermostat;

  gamma = 0.0;
  temperature = 0.0;

  if (!system->rng) {
    throw std::runtime_error("system has no RNG");
  }

  rng = system->rng;

  LOG4ESPP_INFO(theLogger, "Langevin constructed");
}

void LangevinThermostatOnGroup::setGamma(real _gamma) {
  gamma = _gamma;
}

real LangevinThermostatOnGroup::getGamma() {
  return gamma;
}

void LangevinThermostatOnGroup::setTemperature(real _temperature) {
  temperature = _temperature;
}

real LangevinThermostatOnGroup::getTemperature() {
  return temperature;
}

LangevinThermostatOnGroup::~LangevinThermostatOnGroup() {
  disconnect();
}

void LangevinThermostatOnGroup::disconnect() {
  _initialize.disconnect();
  _heatUp.disconnect();
  _coolDown.disconnect();
  _thermalize.disconnect();
}

void LangevinThermostatOnGroup::connect() {
  // connect to initialization inside run()
  _initialize = integrator->runInit.connect(
      boost::bind(&LangevinThermostatOnGroup::initialize, this));

  _heatUp = integrator->recalc1.connect(
      boost::bind(&LangevinThermostatOnGroup::heatUp, this));

  _coolDown = integrator->recalc2.connect(
      boost::bind(&LangevinThermostatOnGroup::coolDown, this));

  _thermalize = integrator->aftCalcF.connect(
      boost::bind(&LangevinThermostatOnGroup::thermalize, this));
}


void LangevinThermostatOnGroup::thermalize() {
  LOG4ESPP_DEBUG(theLogger, "thermalize");

  System &system = getSystemRef();

  for (ParticleGroup::iterator it = particle_group->begin(); it != particle_group->end(); it++) {
    frictionThermo(**it);
  }
}

void LangevinThermostatOnGroup::frictionThermo(Particle &p) {
  real massf = sqrt(p.mass());

  // get a random value for each vector component

  Real3D ranval((*rng)() - 0.5, (*rng)() - 0.5, (*rng)() - 0.5);

  p.force() += pref1 * p.velocity() * p.mass() +
      pref2 * ranval * massf;

  LOG4ESPP_TRACE(theLogger, "new force of p = " << p.force());
}

void LangevinThermostatOnGroup::initialize() {  // calculate the prefactors
  real timestep = integrator->getTimeStep();

  pref1 = -gamma;
  pref2 = sqrt(24.0 * temperature * gamma / timestep);

  LOG4ESPP_INFO(theLogger, "init, timestep = " << timestep <<
      ", gamma = " << gamma <<
      ", temperature = " << temperature << " pref2=" << pref2);
}

/** very nasty: if we recalculate force when leaving/reentering the integrator,
a(t) and a((t-dt)+dt) are NOT equal in the vv algorithm. The random
numbers are drawn twice, resulting in a different variance of the random force.
This is corrected by additional heat when restarting the integrator here.
Currently only works for the Langevin thermostat, although probably also others
are affected.
*/
void LangevinThermostatOnGroup::heatUp() {
  LOG4ESPP_INFO(theLogger, "heatUp");

  pref2buffer = pref2;
  pref2 *= sqrt(3.0);
}

/** Opposite to heatUp */

void LangevinThermostatOnGroup::coolDown() {
  LOG4ESPP_INFO(theLogger, "coolDown");

  pref2 = pref2buffer;
}

/****************************************************
** REGISTRATION WITH PYTHON
****************************************************/
void LangevinThermostatOnGroup::registerPython() {
  using namespace espressopp::python;

  class_<LangevinThermostatOnGroup, shared_ptr<LangevinThermostatOnGroup>, bases<Extension> >
      ("integrator_LangevinThermostatOnGroup",
           init<shared_ptr<System>, shared_ptr<ParticleGroup> >())
      .def("connect", &LangevinThermostatOnGroup::connect)
      .def("disconnect", &LangevinThermostatOnGroup::disconnect)
      .add_property("gamma",
                    &LangevinThermostatOnGroup::getGamma,
                    &LangevinThermostatOnGroup::setGamma)
      .add_property("temperature",
                    &LangevinThermostatOnGroup::getTemperature,
                    &LangevinThermostatOnGroup::setTemperature);
}

}  // end namespace integrator
}  // end namespace espressopp
