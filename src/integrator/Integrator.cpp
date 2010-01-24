
#include "Integrator.hpp"
#include "System.hpp"

using namespace espresso;
using namespace integrator;

LOG4ESPP_LOGGER(Integrator::theLogger, "Integrator");

Integrator::Integrator(shared_ptr<System> _system)
{
  system = _system;
  LOG4ESPP_INFO(theLogger, "construct Integrator");
}

Integrator::~Integrator()
{
  LOG4ESPP_INFO(theLogger, "free Integrator");
}

void Integrator::setTimeStep(double _dt)
{
  dt = _dt;
}
