#define LOG4ESPP_LEVEL_INFO

#include "PeriodicBC.hpp"
#include <cmath>

using namespace espresso;
using namespace espresso::bc;

/* ---------------------------------------------------------------------- */

LOG4ESPP_LOGGER(PeriodicBC::theLogger, "bc.PeriodicBC");

/* ---------------------------------------------------------------------- */

PeriodicBC::PeriodicBC(real _length[3]) {
  LOG4ESPP_DEBUG(theLogger, "constructor, length = " << length);
  setLength(_length);
}

PeriodicBC::~PeriodicBC() {}

void 
PeriodicBC::setLength(real _length[3]) {
  for(int k = 0; k < 3; k++) {
    length[k] = _length[k];
    lengthInverse[k] = 1.0 / length[k];
  }
  LOG4ESPP_INFO(theLogger, "set length = " << length);
}

void 
PeriodicBC::getLength(real &_length[3]) const { 
  for(int k = 0; k < 3; k++)
    _length[k] = length[k];
}
  
void 
PeriodicBC::foldThis(real &pos[3]) const {
  pos[0] -= floor(pos[0] * lengthInverse[0]) * length[0];
  pos[1] -= floor(pos[1] * lengthInverse[1]) * length[1];
  pos[2] -= floor(pos[2] * lengthInverse[2]) * length[2];
}

void 
PeriodicBC::getDist(real dist[3],
                    real &distSqr,
                    const real pos1[3],
                    const real pos2[3]) const {
  for(int k = 0; k < 3; k++)
    dist[k] = pos1[k] - pos2[k];

  dist[0] -= round(dist[0] * lengthInverse[0]) * length[0];
  dist[1] -= round(dist[1] * lengthInverse[1]) * length[1];
  dist[2] -= round(dist[2] * lengthInverse[2]) * length[2];
}

void
PeriodicBC::getRandomPos(real &res[3]) {
  for(int k = 0; k < 3; k++)
    res[k] = length[k];

  res[0] *= drand48();
  res[1] *= drand48();
  res[2] *= drand48();
}
