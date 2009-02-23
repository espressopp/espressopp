#include "FENE.hpp"
#include <cmath>
#include <mpi.hpp>
#include <python.hpp>

using namespace espresso::interaction;

/* ---------------------------------------------------------------------- */

LOG4ESPP_LOGGER(FENE::theLogger, "interaction.FENE");

/* ---------------------------------------------------------------------- */

PMI_REGISTER_CLASS("espresso::interaction::FENE", espresso::interaction::FENE);

FENE::FENE() {
  K = 1.0;
  r0 = 0.0;
  rMax = 1.0;
}
       
FENE::~FENE() {}

PMI_DEFINE_SETTER(espresso::interaction, FENE, setK, real, _K) {
  K = _K;
}
real FENE::getK() const { return K; }

PMI_DEFINE_SETTER(espresso::interaction, FENE, setR0, real, _r0) {
  r0 = _r0;
}
real FENE::getR0() const { return r0; }

PMI_DEFINE_SETTER(espresso::interaction, FENE, setRMax, real, _rMax) {
  rMax = _rMax;
}
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
      .def("setK", &FENE::setK)
      .def("getK", &FENE::getK)
      .def("setR0", &FENE::setR0)
      .def("getR0", &FENE::getR0)
      .def("setRMax", &FENE::setRMax)
      .def("getRMax", &FENE::getRMax)
      .def("computeForce", computeForceOverload)
      .def("computeEnergy", computeEnergyOverload1)
      .def("computeEnergy", computeEnergyOverload2);
    ;
  }
#endif



