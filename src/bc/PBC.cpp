// this macro removes all LOG statements with level < WARN at source level
#define LOG4ESPP_LEVEL_WARN

#include "PBC.hpp"
#include <python.hpp>
#include <iostream>

using namespace espresso::bc;

/* ---------------------------------------------------------------------- */

LOG4ESPP_LOGGER(PBC::theLogger, "bc.PBC");

/* ---------------------------------------------------------------------- */

static inline real dround(real x) { return floor(x + 0.5); }

PMI_REGISTER_CLASS(PBC);

PBC::PBC() {}

PBC::PBC(real _length) {
  setLocal(_length);
}

PBC::~PBC() {}

void PBC::set(real _length)
#ifndef HAVE_MPI
{  setLocal(_length); }
#else
{
  pmiObject.invoke<&PBC::setWorker>();
  boost::mpi::communicator world;
  boost::mpi::broadcast(world, _length, pmi::getControllerMPIRank());
  setLocal(_length);
}

void PBC::setWorker() {
  real v;
  boost::mpi::communicator world;
  boost::mpi::broadcast(world, v, pmi::getControllerMPIRank());
  setLocal(v);
}

PMI_REGISTER_METHOD(PBC, setWorker);
#endif

void PBC::setLocal(real _length) {
  length = _length;
  lengthInverse = 1.0 / length;
}

/** Routine delivers the distance vector between two positions.
    \sa bc::BC::getDist */
Real3D PBC::getDist(const Real3D& pos1, const Real3D& pos2) const {

  real xij;
  real yij;
  real zij;

  xij = pos1[0] - pos2[0];
  yij = pos1[1] - pos2[1];
  zij = pos1[2] - pos2[2];

  xij -= dround(xij*lengthInverse)*length;
  yij -= dround(yij*lengthInverse)*length;
  zij -= dround(zij*lengthInverse)*length;

  return Real3D(xij, yij, zij);
}

#ifdef HAVE_PYTHON
//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void
PBC::registerPython() {
  using namespace boost::python;

  class_<PBC>("bc_PBC", init<>())
    .def("set", &PBC::set)
    .def("getDist", &PBC::getDist);
}
#endif
