#include "FENE.hpp"
#include <cmath>
#include <python.hpp>

using namespace espresso::interaction;

/* ---------------------------------------------------------------------- */

LOG4ESPP_LOGGER(FENE::theLogger, "interaction.FENE");

/* ---------------------------------------------------------------------- */

PMI_REGISTER_CLASS(FENE);

FENE::FENE() { setLocal(1.0, 0.0, 1.0); }
       
FENE::~FENE() {}

void FENE::set(real _K, real _r0, real _rMax)
#ifndef HAVE_MPI
{ setLocal(_K, _r0, _rMax); }
#else
{
  real v[3] = { _K, _r0, _rMax };
  pmiObject.invoke<&FENE::setWorker>();
  boost::mpi::communicator world;
  boost::mpi::broadcast(world, v, 3, pmi::getControllerMPIRank());
  setLocal(v[0], v[1], v[2]);
}

void FENE::setWorker() {
  real v[3];
  boost::mpi::communicator world;
  boost::mpi::broadcast(world, v, 3, pmi::getControllerMPIRank());
  setLocal(v[0], v[1], v[2]);
}

PMI_REGISTER_METHOD(FENE, setWorker);
#endif

void FENE::setLocal(real _K, real _r0, real _rMax) {
  K = _K;
  r0 = _r0;
  rMax = _rMax;
}

real FENE::getK() const { return K; }
real FENE::getR0() const { return r0; }
real FENE::getRMax() const { return rMax; }
      
real FENE::getCutoff() const { return -1; }
real FENE::getCutoffSqr() const { return -1; }

real FENE::computeEnergy (const Real3D &dist,
			  const const_reference p1,
			  const const_reference p2) const {
   return computeEnergy(dist);
}

real FENE::computeEnergy (const Real3D &dist) const {
   return computeEnergySqr(dist.sqr());
}
      
real FENE::computeEnergy (const real dist) const {
   real energy = -0.5 * pow(rMax, 2) * K * log(1 - pow((dist - r0) / rMax, 2));
   return energy;
}

real FENE::computeEnergySqr (real distSqr) const {
    return computeEnergy(sqrt(distSqr));
}
      
Real3D FENE::computeForce (const Real3D &dist,
			   const const_reference p1,
			   const const_reference p2) const {
    return computeForce(dist);
}

Real3D FENE::computeForce (const Real3D &dist) const {

    Real3D f = 0.0;
    real r = sqrt(dist.sqr());

    real ffactor = -K * (r - r0) / (1 - pow((r - r0) / rMax, 2)) / r;
    f = ffactor * dist;

    return f;
}

#ifdef HAVE_PYTHON

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

    class_<FENE>("interaction_FENE", init<>())
      .def("set", &FENE::set)
      .def("getK", &FENE::getK)
      .def("getR0", &FENE::getR0)
      .def("getRMax", &FENE::getRMax)
      .def("computeForce", computeForceOverload)
      .def("computeEnergy", computeEnergyOverload1)
      .def("computeEnergy", computeEnergyOverload2);
    ;
  }
#endif



