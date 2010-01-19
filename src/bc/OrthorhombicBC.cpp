#include "OrthorhombicBC.hpp"
#include <cmath>

namespace espresso {
  namespace bc {
    /* Constructor */
    OrthorhombicBC::OrthorhombicBC(const real _boxL[3]) {
      setBoxL(_boxL);
    }

    /* Setter method for the box length */
    void OrthorhombicBC::setBoxL(const real _boxL[3]) {
      for(int i = 0; i < 3; i++) {
	boxL[i]    = _boxL[i];
	invBoxL[i] = 1.0/_boxL[i];
      }
    }

    /* Returns minimum image vector between two particles */
    void OrthorhombicBC::getMinimumImageVector(real dist[3],
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
    void OrthorhombicBC::foldCoordinate(real pos[3], int imageBox[3], int dir) {
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
    void OrthorhombicBC::foldPosition(real pos[3], int imageBox[3]) {
      for (int i = 0; i < 3; ++i)
	foldCoordinate(pos, imageBox, i);
    }

    /* Unfold coordinates */
    void OrthorhombicBC::unfoldPosition(real pos[3], int imageBox[3]) {
      for (int i = 0; i < 3; ++i) {
	pos[i] = pos[i] + imageBox[i]*getBoxL(i);
	imageBox[i] = 0;
      }
    }

    /* Get random position in the central image box */
    void OrthorhombicBC::getRandomPos(real res[3]) const {
      for(int k = 0; k < 3; k++)
	res[k] = boxL[k];

      // TODO: Use real RNG
      res[0] *= drand48();
      res[1] *= drand48();
      res[2] *= drand48();
    }
  }
}
