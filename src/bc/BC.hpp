#ifndef _BC_BC_HPP
#define _BC_BC_HPP

#include <boost/python/tuple.hpp>
#include "types.hpp"
#include "log4espp.hpp"

namespace espresso {
  namespace bc {
    /** Abstract base class for boundary conditions. */
    class BC {
    public:
      BC(shared_ptr< esutil::RNG > _rng) : rng(_rng) {}

      /** Virtual destructor for boundary conditions. */
      virtual ~BC() {}

      /** Getter for box dimensions */
      virtual Real3D getBoxL() const = 0;

      /** Getter for the RNG. */
      virtual shared_ptr< esutil::RNG > getRng() { return rng; }
      /** Setter for RNG. */
      virtual void setRng(shared_ptr< esutil::RNG > _rng) { rng = _rng; }

      /** Computes the minimum image distance vector between two
	  positions.

	  \param dist is the distance vector (pos2 - pos1)
	  \param pos1, pos2 are the particle positions 
      */
      virtual void
      getMinimumImageVector(Real3DRef dist,
			    const ConstReal3DRef pos1,
			    const ConstReal3DRef pos2) const = 0;

      virtual Real3D 
      getMinimumImageVector(const ConstReal3DRef pos1,
			    const ConstReal3DRef pos2) const;

      /** fold a coordinate to the primary simulation box.
	  \param pos         the position...
	  \param imageBox    and the box
	  \param dir         the coordinate to fold: dir = 0,1,2 for x, y and z coordinate.

	  Both pos and image_box are I/O,
	  i. e. a previously folded position will be folded correctly.
      */
      virtual void 
      foldCoordinate(Real3DRef pos, Int3DRef imageBox, int dir) const = 0;

      virtual void 
      unfoldCoordinate(Real3DRef pos, Int3DRef imageBox, int dir) const = 0;

      /** Fold the coordinates to the primary simulation box.
	  \param pos the position...
	  \param imageBox and the box

	  Both pos and image_box are I/O,
	  i. e. a previously folded position will be folded correctly.
      */
      virtual void 
      foldPosition(Real3DRef pos, Int3DRef imageBox) const;

      virtual void
      foldPosition(Real3DRef pos) const;

      /** Get the folded position as a Python tuple. */
      virtual class boost::python::tuple
      getFoldedPosition(ConstReal3DRef pos, 
			ConstInt3DRef imageBox) const;

      virtual class boost::python::tuple 
      getFoldedPosition(ConstReal3DRef pos) const;

      /** unfold coordinates to physical position.
	  \param pos the position...
	  \param imageBox and the box
	
	  Both pos and image_box are I/O, i.e. image_box will be (0,0,0)
	  afterwards.
      */
      virtual void 
      unfoldPosition(Real3DRef pos, Int3DRef imageBox) const;

      virtual Real3D
      getUnfoldedPosition(ConstReal3DRef pos, ConstInt3DRef imageBox) const;

      /** Get a random position within the central simulation box. The
	  positions are assigned with each coordinate on [0, boxL]. */
      virtual void
      getRandomPos(Real3DRef res) const = 0;

      virtual Real3D
      getRandomPos() const;

      static void registerPython();

     protected:
      shared_ptr< esutil::RNG > rng;

      static LOG4ESPP_DECL_LOGGER(logger);

    }; 
  }
}

#endif
