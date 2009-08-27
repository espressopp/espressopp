#ifndef _ESUTIL_BOX_HPP
#define _ESUTIL_BOX_HPP

#include "Vector3D.hpp"

namespace espresso {
  namespace esutil {
    /** Represents a box with sides parallel to the coordinates
	axes by two corner points.
     */
    template<class T>
    class Box {
    public:
      /// Default constructor
      Box() {}
      /** Constructor from two diagonally located corners.
       */
      Box(const typename espresso::esutil::Vector3D<T> &corner1,
	  const typename espresso::esutil::Vector3D<T> &corner2) {
	for (size_t i = 0; i < 3; ++i) {
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
      const typename espresso::esutil::Vector3D<T> &getLeft()  const { return left; }
      /// get upper coordinates
      const typename espresso::esutil::Vector3D<T> &getRight() const { return right; }

      espresso::esutil::Vector3D<T> getExtend() const { return right - left; }

    protected:
      espresso::esutil::Vector3D<T> left, right;
    };
  }
}
#endif
