#ifndef _BC_HPP
#define _BC_HPP

#include "types.hpp"
#include "Real3DPtr.hpp"
#include "Real3D.hpp"

namespace espresso {
    /** Abstract base class for boundary conditions. */
    class BC {
    public:
      typedef shared_ptr< BC > SelfPtr;

      /** Virtual destructor for boundary conditions. */
      virtual ~BC() {};

      /** Method to set the length of the side of the simulation box */
      virtual void
      setBoxL(Real3DPtr _boxL) = 0;
      
      /** Getters for box dimensions */
      virtual const Real3D &getBoxL() const = 0;
      virtual const Real3D &getInvBoxL() const = 0;
      virtual real getBoxL(int i)      const = 0;
      virtual real getInvBoxL(int i)   const = 0;

      /** Computes the minimum image distance vector between two
          positions. This routine must be implemented by derived
          classes (once the code stabilizes).

          \param dist is the distance vector (pos2 - pos1)
          \param distSqr is the square of the magnitude
                 of the separation vector
          \param pos1, pos2 are the particle positions 
      */
      virtual void
      getMinimumImageVector(Real3DPtr dist,
                            real &distSqr,
                            const Real3DPtr pos1,
                            const Real3DPtr pos2) const = 0;

      /** fold a coordinate to the primary simulation box.
	  \param pos         the position...
	  \param imageBox    and the box
	  \param dir         the coordinate to fold: dir = 0,1,2 for x, y and z coordinate.

	  Both pos and image_box are I/O,
	  i. e. a previously folded position will be folded correctly.
      */
      virtual void foldCoordinate(Real3DPtr pos, int imageBox[3], int dir) = 0;

      /** fold particle coordinates to the primary simulation box.
	  \param pos the position...
	  \param imageBox and the box

	  Both pos and image_box are I/O,
	  i. e. a previously folded position will be folded correctly.
      */
      virtual void foldPosition(Real3DPtr pos, int imageBox[3]) = 0;

      /** unfold coordinates to physical position.
	  \param pos the position...
	  \param imageBox and the box
	
	  Both pos and image_box are I/O, i.e. image_box will be (0,0,0)
	  afterwards.
      */
      virtual void unfoldPosition(Real3DPtr pos, int imageBox[3]) = 0;

      /** Get a random position within the central simulation box. The
          positions are assigned with each coordinate on [0, boxL]. */
      virtual void
      getRandomPos(Real3DPtr res) const = 0;

      virtual Real3D
      getRandomPos() const {
	Real3D res;
	getRandomPos(res);
	return res;
      }
      
      static void registerPython();

    }; 
}

#endif
