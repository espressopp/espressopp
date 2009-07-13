//method to compute the potential energy

// this macro removes all LOG statements with level < DEBUG at source level

#define LOG4ESPP_LEVEL_DEBUG

#include "LennardJones.hpp"
#include <python.hpp>

using namespace espresso;
using namespace espresso::potential;

/* ---------------------------------------------------------------------- */

LOG4ESPP_LOGGER(LennardJones::theLogger, "_espresso.potential.LennardJones");

/* ---------------------------------------------------------------------- */

LennardJones::LennardJones() {}

LennardJones::~LennardJones() {}

void LennardJones::set(real _epsilon, real _sigma, real _cutoff) {

  LOG4ESPP_DEBUG(theLogger, "set epsilon = " << _epsilon << 
                            ", sigma = " << _sigma << ", cutoff = " << _cutoff);

  cutoff = _cutoff;
  computer.epsilon = _epsilon;
  computer.sigma = _sigma;
  computer.cutoffSqr = cutoff*cutoff;
}

real LennardJones::getCutoff() const { return cutoff; }
real LennardJones::getEpsilon() const { return computer.epsilon; }
real LennardJones::getSigma() const { return computer.sigma; }
real LennardJones::getCutoffSqr() const { return computer.cutoffSqr; }

real LennardJones::computeEnergySqr(const real distSqr) const {
  return computer.computeEnergySqr(distSqr);
}

Real3D LennardJones::computeForce(const Real3D &dist) const {
  return computer.computeForce(dist);
}

real LennardJones::BasicComputer::computeEnergySqr(const real distSqr) const {
  if (distSqr < cutoffSqr) {
    real frac2 = sigma*sigma / distSqr;
    real frac6 = frac2 * frac2 * frac2;
    real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
    return energy;
  } else {
    return 0.0;
  }
}
    
Real3D LennardJones::BasicComputer::computeForce(const Real3D &dist) const {
  Real3D f = 0.0;
  real   frac2;
  real   frac6;
      
  real distSqr = dist.sqr();
      
  if (distSqr < cutoffSqr) {
    frac2 = sigma*sigma / distSqr;
    frac6 = frac2 * frac2 * frac2;
    real ffactor = 48.0 * epsilon * (frac6*frac6 - 0.5 * frac6) / distSqr;

    LOG4ESPP_TRACE(LennardJones::theLogger, "computeForce, distSqr = " << distSqr <<
		   ", ffactor = " << ffactor);

    f = dist * ffactor;
  }
      
  return f;
}

pairs::EnergyComputer*
LennardJones::createEnergyComputer(const pairs::EnergyComputer &templ) const
{ return new pairs::SquareDistEnergyComputerFacade<LennardJones::BasicComputer>(templ, computer); }

pairs::ForceComputer*
LennardJones::createForceComputer(const pairs::ForceComputer &templ) const
{ return new pairs::VectorForceComputerFacade<LennardJones::BasicComputer>(templ, computer); }
   
//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void 
LennardJones::registerPython() {
  using namespace espresso::python;

  class_< LennardJones, bases< CentralPotential > >
    ("potential_LennardJones")
    .def("set", &LennardJones::set)
    .def("getCutoff", &LennardJones::getCutoff)
    .def("getEpsilon", &LennardJones::getEpsilon)
    .def("getSigma", &LennardJones::getSigma)
    ;
}
