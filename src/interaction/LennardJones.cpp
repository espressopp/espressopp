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

PMI_REGISTER_CLASS(LennardJones);

LennardJones::LennardJones() {
  // set the defaults here
  setLocal(1.0, 1.0, 2.5);
}

LennardJones::~LennardJones() {}

void LennardJones::set(real _epsilon, real _sigma, real _cutoff) {
#ifndef HAVE_MPI
  setLocal(_epsilon, _sigma, _cutoff);
#else
  real v[3] = { _epsilon, _sigma, _cutoff };
  pmiObject.invoke<&LennardJones::setWorker>();
  boost::mpi::communicator world;
  boost::mpi::broadcast(world, v, 3, pmi::getControllerMPIRank());
  setLocal(v[0], v[1], v[2]);
  }

void LennardJones::setWorker() {
  real v[3];
  boost::mpi::communicator world;
  boost::mpi::broadcast(world, v, 3, pmi::getControllerMPIRank());
  setLocal(v[0], v[1], v[2]);
}

PMI_REGISTER_METHOD(LennardJones, setWorker);
#endif

void LennardJones::setLocal(real _epsilon, real _sigma, real _cutoff) {
  epsilon = _epsilon;
  sigma = _sigma;
  cutoff = _cutoff;
  cutoffSqr = cutoff*cutoff;
}

real LennardJones::getCutoff() const { return cutoff; }
real LennardJones::getCutoffSqr() const { return cutoffSqr; }
real LennardJones::getEpsilon() const { return epsilon; }
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
    .def("set", &LennardJones::set)
    .def("getCutoff", &LennardJones::getCutoff)
    .def("getEpsilon", &LennardJones::getEpsilon)
    .def("getSigma", &LennardJones::getSigma)
    .def("computeForce", computeForceOverload)
    .def("computeEnergy", computeEnergyOverload1)
    .def("computeEnergy", computeEnergyOverload2);
}
#endif
