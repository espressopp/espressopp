#include "python.hpp"
#include "OrthorhombicBC.hpp"
#include <cmath>
#include "Real3D.hpp"
#include "Int3D.hpp"
#include "esutil/RNG.hpp"

namespace espresso {
  namespace bc {
    /* Constructor */
    OrthorhombicBC::
    OrthorhombicBC(shared_ptr< esutil::RNG > _rng,
		   const Real3D& _boxL) 
      : BC(_rng)
    { setBoxL(_boxL); }
  
    /* Setter method for the box length */
    void OrthorhombicBC::setBoxL(const Real3D& _boxL) {
      boxL = _boxL;
      for (int i = 0; i < 3; i++) {
	    invBoxL[i] = 1.0/boxL[i];
	    boxL2[i] = 0.5*boxL[i];
      }
      onBoxDimensionsChanged();
    }
    void OrthorhombicBC::scaleVolume(real s) {
  	  boxL *= s;
  	  boxL2 *= s;
  	  invBoxL /= s;
  	  onBoxDimensionsChanged();
    }

    /* Returns the minimum image vector between two positions */
    void 
    OrthorhombicBC::
    getMinimumImageVector(Real3D& dist,
			  const Real3D& pos1,
			  const Real3D& pos2) const {
      dist = pos1;
      dist -= pos2;

      dist[0] -= round(dist[0] * invBoxL[0]) * boxL[0];
      dist[1] -= round(dist[1] * invBoxL[1]) * boxL[1];
      dist[2] -= round(dist[2] * invBoxL[2]) * boxL[2];
    }

    /* Returns the minimum image vector between two positions */
    void 
    OrthorhombicBC::
    getMinimumImageVectorBox(Real3D& dist,
                             const Real3D& pos1,
                             const Real3D& pos2) const {
      dist = pos1;
      dist -= pos2;

      if (dist[0] < -boxL2[0]) dist[0] += boxL[0];
      else if (dist[0] > boxL2[0]) dist[0] -= boxL[0];
      if (dist[1] < -boxL2[1]) dist[1] += boxL[1];
      else if (dist[1] > boxL2[1]) dist[1] -= boxL[1];
      if (dist[2] < -boxL2[2]) dist[2] += boxL[2];
      else if (dist[2] > boxL2[2]) dist[2] -= boxL[2];
    }

    /* Fold back a nearby position into box */

    void
    OrthorhombicBC::getMinimumDistance(Real3D& dist) const {

      if (dist[0] < -0.5 * boxL[0]) dist[0] += boxL[0];
      else if (dist[0] > 0.5 * boxL[0]) dist[0] -= boxL[0];
      if (dist[1] < -0.5 * boxL[1]) dist[1] += boxL[1];
      else if (dist[1] > 0.5 * boxL[1]) dist[1] -= boxL[1];
      if (dist[2] < -0.5 * boxL[2]) dist[2] += boxL[2];
      else if (dist[2] > 0.5 * boxL[2]) dist[2] -= boxL[2];
    }

    /* Returns the minimum image vector between two positions */
    void
    OrthorhombicBC::
    getMinimumImageVectorX(real dist[3],
                          const real pos1[3],
                          const real pos2[3]) const {

      dist[0] = pos1[0];
      dist[1] = pos1[1];
      dist[2] = pos1[2];

      dist[0] -= pos2[0];
      dist[1] -= pos2[1];
      dist[2] -= pos2[2];

      dist[0] -= round(dist[0] * invBoxL[0]) * boxL[0];
      dist[1] -= round(dist[1] * invBoxL[1]) * boxL[1];
      dist[2] -= round(dist[2] * invBoxL[2]) * boxL[2];
    }

    /* Fold an individual coordinate in the specified direction */
    void 
    OrthorhombicBC::
    foldCoordinate(Real3D& pos, Int3D& imageBox, int dir) const {
      int tmp = static_cast<int>(floor(pos[dir]*invBoxL[dir]));

      imageBox[dir] += tmp;
      pos[dir] -= tmp*boxL[dir];

      if(pos[dir] < 0 || pos[dir] >= boxL[dir]) {
	/* slow but safe */
	if (fabs(pos[dir]*invBoxL[dir]) >= INT_MAX/2) {
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

    /* Unfold an individual coordinate in the specified direction */
    void 
    OrthorhombicBC::
    unfoldCoordinate(Real3D& pos, Int3D& imageBox, int dir) const {
      pos[dir] += imageBox[dir]*boxL[dir];
      imageBox[dir] = 0;
    }

    /* Get random position in the central image box */
    void 
    OrthorhombicBC::
    getRandomPos(Real3D& res) const {
      for(int k = 0; k < 3; k++)
	res[k] = boxL[k];
      
      res[0] *= (*rng)();
      res[1] *= (*rng)();
      res[2] *= (*rng)();
    }
    
    void 
    OrthorhombicBC::
    registerPython() {
      using namespace espresso::python;
      class_<OrthorhombicBC, bases< BC >, boost::noncopyable >
	("bc_OrthorhombicBC", init< shared_ptr< esutil::RNG >, Real3D& >())
	.add_property("boxL", &OrthorhombicBC::getBoxL, &OrthorhombicBC::setBoxL)
      ;
    }
  }
}
