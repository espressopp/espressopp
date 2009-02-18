#include "FENE.hpp"
#include <cmath>
#include <mpi.hpp>
#include <python.hpp>

using namespace espresso::interaction;

/* ---------------------------------------------------------------------- */

LOG4ESPP_LOGGER(FENE::theLogger, "interaction.FENE");

/* ---------------------------------------------------------------------- */

FENE::FENE() {}
       
FENE::~FENE() {}

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
      
      
void FENE::setK(real _K) {

#ifdef HAVE_MPI
      // invoke setK for the Workers
      pmiObject.invoke<&FENE::setKWorker>();
      // broadcast _K
      boost::mpi::communicator world;
      boost::mpi::broadcast(world, _K, pmi::getControllerMPIRank());
#endif
      // invoke setK in SPMD mode
      setKLocal(_K);
}

#ifdef HAVE_MPI
void FENE::setKWorker() {

      real _K;
      // broadcast _K
      boost::mpi::communicator world;
      boost::mpi::broadcast(world, _K, pmi::getControllerMPIRank());
      // invoke setK in SPMD mode
      setKLocal(_K);
}
#endif

void FENE::setr0(real _r0) { 

#ifdef HAVE_MPI
      // invoke setr0 for the Workers
      pmiObject.invoke<&FENE::setr0Worker>();
      // broadcast _r0
      boost::mpi::communicator world;
      boost::mpi::broadcast(world, _r0, pmi::getControllerMPIRank());
#endif
      // invoke setr0 in SPMD mode
      setr0Local(_r0);
}

#ifdef HAVE_MPI
void FENE::setr0Worker() {

      real _r0;
      // broadcast _r0
      boost::mpi::communicator world;
      boost::mpi::broadcast(world, _r0, pmi::getControllerMPIRank());
      // invoke setr0 in SPMD mode
      setr0Local(_r0);
}
#endif

void FENE::setrMax(real _rMax) { 

#ifdef HAVE_MPI
      // invoke setrMax for the Workers
      pmiObject.invoke<&FENE::setrMaxWorker>();
      // broadcast _rMax
      boost::mpi::communicator world;
      boost::mpi::broadcast(world, _rMax, pmi::getControllerMPIRank());
#endif
      // invoke setrMax in SPMD mode
      setrMaxLocal(_rMax);
}

#ifdef HAVE_MPI
void FENE::setrMaxWorker() {

      real _rMax;
      // broadcast _rMax
      boost::mpi::communicator world;
      boost::mpi::broadcast(world, _rMax, pmi::getControllerMPIRank());
      // invoke setrMax in SPMD mode
      setrMaxLocal(_rMax);
}
#endif

void FENE::setKLocal(real _K) { 
  LOG4ESPP_INFO(theLogger, "setK (local) : " << _K);
  K = _K; 
}

void FENE::setr0Local(real _r0) { r0 = _r0; }
void FENE::setrMaxLocal(real _rMax) { rMax = _rMax; }
      
real FENE::getCutoff() const { return -1; }
real FENE::getCutoffSqr() const { return -1; }

//////////////////////////////////////////////////
// REGISTRATION WITH PMI
//////////////////////////////////////////////////

PMI_REGISTER_CLASS("espresso::interaction::FENE", espresso::interaction::FENE);
PMI_REGISTER_METHOD("setKWorker", espresso::interaction::FENE, setKWorker);
PMI_REGISTER_METHOD("setr0Worker", espresso::interaction::FENE, setr0Worker);
PMI_REGISTER_METHOD("setrMaxWorker", espresso::interaction::FENE, setrMaxWorker);

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
      .def("setr0", &FENE::setr0)
      .def("setrMax", &FENE::setrMax)
      .def("computeForce", computeForceOverload)
      .def("computeEnergy", computeEnergyOverload1)
      .def("computeEnergy", computeEnergyOverload2);
    ;
  }
#endif



