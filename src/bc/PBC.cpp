// this macro removes all LOG statements with level < WARN at source level
#define LOG4ESPP_LEVEL_WARN

#include "PBC.hpp"
#include <cmath>
#include <python.hpp>
#include <iostream>

using namespace espresso;
using namespace espresso::bc;

/* ---------------------------------------------------------------------- */

LOG4ESPP_LOGGER(PBC::theLogger, "bc.PBC");

/* ---------------------------------------------------------------------- */

static inline real dround(real x) { return floor(x + 0.5); }

PBC::PBC() {}

PBC::PBC(real _length) {
  set(_length);
}

PBC::~PBC() {}

void PBC::set(real _length) {
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
