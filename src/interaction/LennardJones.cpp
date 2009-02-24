//method to compute the potential energy

#include "LennardJones.hpp"
#include <mpi.hpp>
#include <pmi.hpp>
#include <python.hpp>
#include <iostream>

using namespace espresso::interaction;

/* ---------------------------------------------------------------------- */

LOG4ESPP_LOGGER(LennardJones::theLogger, "interaction.LennardJones");

/* ---------------------------------------------------------------------- */

PMI_REGISTER_CLASS("espresso::interaction::LennardJones", espresso::interaction::LennardJones);

LennardJones::LennardJones() {
  epsilon = 1.0;
  sigma = 1.0;
  cutoff = 2.5;
  cutoffSqr = 2.5*2.5;
}

LennardJones::~LennardJones() {}

PMI_DEFINE_SETTER(espresso::interaction, LennardJones, setCutoff, real, _cutoff) {
  cutoff = _cutoff; 
  cutoffSqr = cutoff * cutoff;
}
real LennardJones::getCutoff() const { return cutoff; }
real LennardJones::getCutoffSqr() const { return cutoffSqr; }
    
PMI_DEFINE_SETTER(espresso::interaction, LennardJones, setEpsilon, real, _epsilon) {
  epsilon = _epsilon;
}
real LennardJones::getEpsilon() const { return epsilon; }

PMI_DEFINE_SETTER(espresso::interaction, LennardJones, setSigma, real, _sigma) {
  sigma = _sigma;
}
real LennardJones::getSigma() const { return sigma; }


real LennardJones::computeEnergy (const Real3D &dist,
				  const const_reference p1,
				  const const_reference p2) const {
  return computeEnergy(dist);
}

real LennardJones::computeEnergy (const Real3D &dist) const {
  return computeEnergySqr(dist.sqr());
}

real LennardJones::computeEnergy(const real dist) const {
  return computeEnergySqr(dist*dist);
}
    
real LennardJones::computeEnergySqr (const real distSqr) const {
  if (distSqr < cutoffSqr) {
    real frac2 = sigma*sigma / distSqr;
    real frac6 = frac2 * frac2 * frac2;
    real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
    return energy;
  } else return 0.0;
}
    
Real3D LennardJones::computeForce (const Real3D &dist,
				   const const_reference p1,
				   const const_reference p2) const {
  return computeForce(dist);
}

Real3D LennardJones::computeForce (const Real3D &dist) const {
  Real3D f = 0.0;
  real   frac2;
  real   frac6;
      
  real distSqr = dist.sqr();
      
  if (distSqr < cutoffSqr) {
    frac2 = sigma / distSqr;
    frac6 = frac2 * frac2 * frac2;
    real ffactor = 48.0 * epsilon * (frac6*frac6 - 0.5 * frac6) * frac2;
	
    LOG4ESPP_DEBUG(theLogger, "computeForce, distSqr = " << distSqr <<
		   ", ffactor = " << ffactor);
	
    f = dist * ffactor;
  } 
      
  return f;
}
    


#ifdef HAVE_PYTHON  
//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void 
LennardJones::registerPython() {
  using namespace boost::python;

  // create thin wrappers around overloaded member functions
  Real3D (LennardJones::*computeForceOverload)(const Real3D&) const =
    &LennardJones::computeForce;
  real (LennardJones::*computeEnergyOverload1)(const Real3D &) const =
    &LennardJones::computeEnergy;
  real (LennardJones::*computeEnergyOverload2)(const real) const =
    &LennardJones::computeEnergy;
    
  class_<LennardJones>("interaction_LennardJones", init<>())
    .def("getCutoff", &LennardJones::getCutoff)
    .def("setCutoff", &LennardJones::setCutoff)
    .def("getEpsilon", &LennardJones::getEpsilon)
    .def("setEpsilon", &LennardJones::setEpsilon)
    .def("getSigma", &LennardJones::getSigma)
    .def("setSigma", &LennardJones::setSigma)
    .def("computeForce", computeForceOverload)
    .def("computeEnergy", computeEnergyOverload1)
    .def("computeEnergy", computeEnergyOverload2);
}
#endif
