#include "FENE.hpp"
#include <cmath>
#include <python.hpp>

using namespace espresso;
using namespace espresso::potential;

/* ---------------------------------------------------------------------- */

LOG4ESPP_LOGGER(FENE::theLogger, "potential.FENE");

/* ---------------------------------------------------------------------- */

FENE::FENE() {}
       
FENE::~FENE() {}

void FENE::set(real _K, real _r0, real _rMax) {
  computer.K = _K;
  computer.r0 = _r0;
  computer.rMax = _rMax;
}

real FENE::getK() const { return computer.K; }
real FENE::getR0() const { return computer.r0; }
real FENE::getRMax() const { return computer.rMax; }
real FENE::getCutoffSqr() const { return -1.0; }

real FENE::computeEnergySqr (const real distSqr) const {
  return computer.computeEnergySqr(distSqr);
}

Real3D FENE::computeForce (const Real3D &dist) const {
  return computer.computeForce(dist);
}

real FENE::BasicComputer::computeEnergy(const real r) const {
  real energy = -0.5 * pow(rMax, 2) * K * log(1 - pow((r - r0) / rMax, 2));
  return energy;
}

real FENE::BasicComputer::computeEnergySqr(const real dist) const {
  return computeEnergy(sqrt(dist));
}

Real3D FENE::BasicComputer::computeForce (const Real3D &dist) const {
  Real3D f = 0.0;
  real r = sqrt(dist.sqr());

  real ffactor = -K * (r - r0) / (1 - pow((r - r0) / rMax, 2)) / r;
  f = ffactor * dist;

  return f;
}

pairs::EnergyComputer*
FENE::createEnergyComputer(const pairs::EnergyComputer &templ) const
{ return new pairs::SquareDistEnergyComputerFacade<FENE::BasicComputer>(templ, computer); }

pairs::ForceComputer*
FENE::createForceComputer(const pairs::ForceComputer &templ) const
{ return new pairs::VectorForceComputerFacade<FENE::BasicComputer>(templ, computer); }


//////////////////////////////////////////////////
  // REGISTRATION WITH PYTHON
  //////////////////////////////////////////////////

void
FENE::registerPython() {
  using namespace boost::python;
  
  class_< FENE, bases< CentralPotential > >("potential_FENE")
    .def("set", &FENE::set)
    .def("getK", &FENE::getK)
    .def("getR0", &FENE::getR0)
    .def("getRMax", &FENE::getRMax)
    ;
}



