#include "FENE.hpp"
#include <cmath>
#include <python.hpp>

using namespace espresso;
using namespace espresso::interaction;

/* ---------------------------------------------------------------------- */

LOG4ESPP_LOGGER(FENE::theLogger, "interaction.FENE");

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

real FENE::computeEnergy (const Real3D &dist) const {
  return computer.computeEnergySqr(dist.sqr());
}
      
real FENE::computeEnergy (const real dist) const {
  return computer.computeEnergySqr(dist*dist);
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
  
  Real3D (FENE::*computeForceOverload)(const Real3D&) const =
    &FENE::computeForce;
  real (FENE::*computeEnergyOverload1)(const Real3D &) const =
    &FENE::computeEnergy;
  real (FENE::*computeEnergyOverload2)(const real) const =
    &FENE::computeEnergy;
  
  class_<FENE, boost::shared_ptr<FENE>, bases<Interaction> >("interaction_FENE", init<>())
    .def("set", &FENE::set)
    .def("getK", &FENE::getK)
    .def("getR0", &FENE::getR0)
    .def("getRMax", &FENE::getRMax)
    .def("computeForce", computeForceOverload)
    .def("computeEnergy", computeEnergyOverload1)
    .def("computeEnergy", computeEnergyOverload2);
}



