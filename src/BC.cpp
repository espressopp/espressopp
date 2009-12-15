#include "BC.hpp"
#include <cmath>

using namespace espresso;
using namespace espresso::bc;

/* Constructor for cubic box */
BC::BC(const real _boxL[3]) {
  setBoxL(_boxL);
}

/* Setter method for the box length */
void BC::setBoxL(const real _boxL[3]) {
  for(int i = 0; i < 3; i++) {
    boxL[i]    = _boxL[i];
    invBoxL[i] = 1.0/_boxL[i];
  }
}

/* Returns minimum image vector between two particles */
void BC::getMinimumImageVector(real dist[3],
                    real &distSqr,
                    const real pos1[3],
                    const real pos2[3]) const {
  for(int k = 0; k < 3; k++)
    dist[k] = pos1[k] - pos2[k];

  dist[0] -= round(dist[0] * invBoxL[0]) * boxL[0];
  dist[1] -= round(dist[1] * invBoxL[1]) * boxL[1];
  dist[2] -= round(dist[2] * invBoxL[2]) * boxL[2];
}

/* Fold an individual coordinate in the specified direction */
void BC::foldCoordinate(real pos[3], int imageBox[3], int dir) {
  int tmp = static_cast<int>(floor(pos[dir]*getInvBoxL(dir)));

  imageBox[dir] += tmp;
  pos[dir] -= tmp*getBoxL(dir);

  if(pos[dir] < 0 || pos[dir] >= getBoxL(dir)) {
    /* slow but safe */
    if (fabs(pos[dir]*getInvBoxL(dir)) >= INT_MAX/2) {
# warning ERRORHANDLING MISSING
#if 0
      char *errtext = runtime_error(128 + TCL_INTEGER_SPACE);
      ERROR_SPRINTF(errtext,"{086 particle coordinate out of range, pos = %f, image box = %d} ", pos[dir], image_box[dir]);
#endif
      imageBox[dir] = 0;
      pos[dir] = 0;
    }
  }
}

/* Fold coordinates */
void BC::foldPosition(real pos[3], int imageBox[3]) {
  for(int i = 0; i < 3; ++i)
    foldCoordinate(pos, imageBox, i);
}

/* Unfold coordinates */
void BC::unfoldPosition(real pos[3], int imageBox[3]) {
  for(int i = 0; i < 3; ++i) {
    pos[i] = pos[i] + imageBox[i]*getBoxL(i);
    imageBox[i] = 0;
  }
}

/* Wrap a particle to the central image */
void BC::foldThis(real pos[3]) const {
  pos[0] -= floor(pos[0] * invBoxL[0]) * boxL[0];
  pos[1] -= floor(pos[1] * invBoxL[1]) * boxL[1];
  pos[2] -= floor(pos[2] * invBoxL[2]) * boxL[2];
}

/* Wrap a copy of the coordinates of a particle to the central image */
/* It might be better to pass in pos[3] and a second array which will
   store the folded coordinates. This would avoid the warning that a
   local variable is being returned */
real* BC::fold(const real pos[3]) const {
  real res[3] = {pos[0], pos[1], pos[2]};
  foldThis(res);
  return res;
}

/* Get random position in the central image box */
void BC::getRandomPos(real res[3]) const {
  for(int k = 0; k < 3; k++)
    res[k] = boxL[k];

  res[0] *= drand48();
  res[1] *= drand48();
  res[2] *= drand48();
}
