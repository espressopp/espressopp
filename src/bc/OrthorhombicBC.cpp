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
		   const ConstReal3DRef _boxL) 
      : BC(_rng)
    { setBoxL(_boxL); }
  
    /* Setter method for the box length */
    void OrthorhombicBC::setBoxL(const ConstReal3DRef _boxL) {
      boxL = _boxL;
      for (int i = 0; i < 3; i++)
	invBoxL[i] = 1.0/boxL[i];
    }

    /* Returns the minimum image vector between two positions */
    void 
    OrthorhombicBC::
    getMinimumImageVector(Real3DRef dist,
			  const ConstReal3DRef pos1,
			  const ConstReal3DRef pos2) const {
      dist = pos1;
      dist -= pos2;

      dist[0] -= round(dist[0] * invBoxL[0]) * boxL[0];
      dist[1] -= round(dist[1] * invBoxL[1]) * boxL[1];
      dist[2] -= round(dist[2] * invBoxL[2]) * boxL[2];
    }

    /* Fold an individual coordinate in the specified direction */
    void 
    OrthorhombicBC::
    foldCoordinate(Real3DRef pos, Int3DRef imageBox, int dir) const {
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
    unfoldCoordinate(Real3DRef pos, Int3DRef imageBox, int dir) const {
      pos[dir] += imageBox[dir]*boxL[dir];
      imageBox[dir] = 0;
    }

    /* Get random position in the central image box */
    void 
    OrthorhombicBC::
    getRandomPos(Real3DRef res) const {
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
      class_<OrthorhombicBC, bases< BC > >
	("bc_OrthorhombicBC", init< shared_ptr< esutil::RNG >, Real3DRef >())
	.add_property("boxL", &OrthorhombicBC::getBoxL, &OrthorhombicBC::setBoxL)
      ;
    }
  }
}
