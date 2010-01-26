
#include "MDIntegrator.hpp"
#include "System.hpp"

using namespace espresso;
using namespace integrator;

LOG4ESPP_LOGGER(MDIntegrator::theLogger, "MDIntegrator");

MDIntegrator::MDIntegrator(shared_ptr<System> _system)
{
  system = _system;
  LOG4ESPP_INFO(theLogger, "construct Integrator");
}

MDIntegrator::~MDIntegrator()
{
  LOG4ESPP_INFO(theLogger, "free Integrator");
}

void MDIntegrator::setTimeStep(real _dt)
{
  dt = _dt;
}
