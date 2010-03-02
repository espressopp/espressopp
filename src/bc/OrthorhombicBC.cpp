#include "OrthorhombicBC.hpp"
#include <python.hpp>
#include <cmath>
#include "Real3D.hpp"
#include "esutil/RNG.hpp"
#include "System.hpp"
#include "stdio.h"

namespace espresso {
  namespace bc {
    /* Constructor */
    OrthorhombicBC::
    OrthorhombicBC(shared_ptr< System > _system, 
		   const ConstReal3DRef _boxL)
    { setBoxL(_boxL); 
      system = _system;
      if (!system.use_count()) {
        printf("ERROR in OrthorhomicBC::Orthorhombic : no system !\n");
      } else {
        printf("OrthorhombicBC::OrthorhombicBC() : system.use_count()=%i\n",system.use_count());
        try {
          esutil::RNG &rng = *(system.lock()->rng);
          printf("OrthorhombicBC::OrthorhombicBC() : system.lock()->rng()=%f\n",rng());
        }
        catch (boost::bad_weak_ptr) {
          printf("WARNING in OrthorhombicBC::Orthorhombic() : system has no RNG !\n");
        }
      }
    }
  

    /* Setter method for the box length */
    void OrthorhombicBC::setBoxL(const ConstReal3DRef _boxL) {
      boxL = _boxL;
      for (int i = 0; i < 3; i++)
	invBoxL[i] = 1.0/boxL[i];
    }

    /* Returns minimum image vector between two particles */
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
    foldCoordinate(Real3DRef pos, int imageBox[3], int dir) {
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
    void 
    OrthorhombicBC::
    foldPosition(Real3DRef pos, int imageBox[3]) {
      for (int i = 0; i < 3; ++i)
	foldCoordinate(pos, imageBox, i);
    }

    /* Unfold coordinates */
    void 
    OrthorhombicBC::
    unfoldPosition(Real3DRef pos, int imageBox[3]) {
      for (int i = 0; i < 3; ++i) {
	pos[i] = pos[i] + imageBox[i]*getBoxL(i);
	imageBox[i] = 0;
      }
    }

    /* Get random position in the central image box */
    void 
    OrthorhombicBC::
    getRandomPos(Real3DRef res) const {
      for(int k = 0; k < 3; k++)
	res[k] = boxL[k];
      
      if (!system.use_count()) {
        printf("ERROR in OrthorhomicBC::getRandomPos() : BC lost system !\n");
      } else {
        try {
          esutil::RNG &rng = *(system.lock()->rng);
          res[0] *= rng();
          res[1] *= rng();
          res[2] *= rng();
        }
        catch (boost::bad_weak_ptr) {
          printf("ERROR in OrthorhombicBC::getRandomPos() : system lost RNG !\n");
        }
      }
    }

    void 
    OrthorhombicBC::
    registerPython() {
      using namespace espresso::python;
      class_<OrthorhombicBC, bases< BC > >
	("bc_OrthorhombicBC", init< shared_ptr< System >, Real3DRef >())
      ;
    }
  }
}
