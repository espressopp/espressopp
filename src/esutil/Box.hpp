#ifndef _ESUTIL_BOX_HPP
#define _ESUTIL_BOX_HPP

#include "Vector3D.hpp"

namespace espresso {
  namespace esutil {
    /** Represents a box with sides parallel to the coordinates
	axes by two corner points.
     */
    template<class VectorClass>
    class Box {
    public:
      /// Default constructor
      Box(): left(0.0), right(0.0) {}
      /** Constructor from two diagonally located corners.
       */
      Box(const VectorClass &corner1, const VectorClass &corner2) {
	for (size_t i = 0; i < VectorClass::dimension; ++i) {
	  if (corner1[i] < corner2[i]) {
	    left[i]  = corner1[i];
	    right[i] = corner2[i];
	  }
	  else {
	    left[i]  = corner2[i];
	    right[i] = corner1[i];
	  }
	}
      }

      /// get lower coordinates
      const VectorClass &getLeft()  const { return left; }
      /// get upper coordinates
      const VectorClass &getRight() const { return right; }

      VectorClass getExtend() const { return right - left; }

    protected:
      VectorClass left, right;
    };
  }
}
#endif
