#include "System.hpp"
#include "Storage.hpp"
#include "Langevin.hpp"
#include "esutil/RNG.hpp"

#include "types.hpp"

using namespace espresso;
using namespace espresso::integrator;

LOG4ESPP_LOGGER(Langevin::theLogger, "Langevin");

Langevin::Langevin(shared_ptr<System> _system) :

  rng(_system->rng)

{
  system = _system;
  gamma  = 0.0;
  temperature = 0.0;

  LOG4ESPP_INFO(theLogger, "Langevin constructed");
}

void Langevin::setGamma(real _gamma)
{
  gamma = _gamma;
}

real Langevin::getGamma()
{
  return gamma;
}

void Langevin::setTemperature(real _temperature)
{
  temperature = _temperature;
}

real Langevin::getTemperature()
{
  return temperature;
}

Langevin::~Langevin()
{
}

void Langevin::thermalize()
{
  LOG4ESPP_INFO(theLogger, "thermalize")

  std::vector< Cell >& localCells = system.lock()->storage->getLocalCells();

  for (size_t c = 0; c < localCells.size(); c++) {
    Cell* localCell = &localCells[c];
    for (size_t index = 0; index < localCell->particles.size(); index++) {
      Particle* particle  = &localCell->particles[index];
      frictionThermo(particle);
    }
  }
}

void Langevin::frictionThermo(Particle *p)
{
  // double massf = sqrt(PMASS(*p));
 
  double massf = 1.0;

  for (int j = 0 ; j < 3 ; j++) {
    p->f.f[j] += pref1*p->m.v[j] + pref2*(rng()-0.5)*massf;
  }

  LOG4ESPP_DEBUG(theLogger, "new force of p = " << p->f.f[0] << " "
                             << p->f.f[1] << " " << p->f.f[2]);
}

void Langevin::init(real timestep)

{ // calculate the prefactors

  LOG4ESPP_INFO(theLogger, "init, timestep = " << timestep <<
                           ", gamma = " << gamma << 
                           ", temperature = " << temperature);

  pref1 = -gamma * timestep;
  pref2 = sqrt(24.0 * temperature * gamma / timestep);

}

/** very nasty: if we recalculate force when leaving/reentering the integrator,
    a(t) and a((t-dt)+dt) are NOT equal in the vv algorithm. The random
    numbers are drawn twice, resulting in a different variance of the random force.
    This is corrected by additional heat when restarting the integrator here.
    Currently only works for the Langevin thermostat, although probably also others
    are affected.
*/

void Langevin::heatUp()
{
  LOG4ESPP_INFO(theLogger, "heatUp");

  pref2buffer = pref2;
  pref2       *= sqrt(3.0);
}

/** Opposite to heatUp */

void Langevin::coolDown()
{
  LOG4ESPP_INFO(theLogger, "coolDown");

  pref2 = pref2buffer;
}

