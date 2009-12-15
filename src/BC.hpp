#ifndef _BC_BC_HPP
#define _BC_BC_HPP

#include "types.hpp"

namespace espresso {
  namespace bc {
    
    /** This will become an abstract class. */
    class BC {
    private:
      real boxL[3];
      real invBoxL[3];

    public:
      typedef shared_ptr< BC > SelfPtr;

      /** Virtual destructor for boundary conditions. */
      virtual
      ~BC() {};

      /** Constructor for cubic box */
      BC(const real _boxL[3]);

      /** Method to set the length of the side of the cubic simulation cell */
      virtual void
      setBoxL(const real _boxL[3]);

      /** Getters for box dimensions */
      const real *getBoxL()    const { return boxL; }
      const real *getInvBoxL() const { return invBoxL; }
      real getBoxL(int i)      const { return boxL[i]; }
      real getInvBoxL(int i)   const { return invBoxL[i]; }

      /** Computes the minimum image distance vector between two
          positions. This routine must be implemented by derived
          classes (once the code stabilizes).

          \param dist is the distance vector (pos2 - pos1)
          \param distSqr is the square of the magnitude
                 of the separation vector
          \param pos1, pos2 are the particle positions 
      */
      virtual void
      getMinimumImageVector(real dist[3],
                            real &distSqr,
                            const real pos1[3],
                            const real pos2[3]) const;

      /** fold a coordinate to primary simulation box.
	  \param pos         the position...
	  \param imageBox    and the box
	  \param dir         the coordinate to fold: dir = 0,1,2 for x, y and z coordinate.

	  Both pos and image_box are I/O,
	  i. e. a previously folded position will be folded correctly.
      */
      void foldCoordinate(real pos[3], int imageBox[3], int dir);

      /** fold particle coordinates to primary simulation box.
	  \param pos the position...
	  \param imageBox and the box

	  Both pos and image_box are I/O,
	  i. e. a previously folded position will be folded correctly.
      */
      void foldPosition(real pos[3], int imageBox[3]);

      /** unfold coordinates to physical position.
	  \param pos the position...
	  \param imageBox and the box
	
	  Both pos and image_box are I/O, i.e. image_box will be (0,0,0)
	  afterwards.
      */
      void unfoldPosition(real pos[3], int imageBox[3]);


      /** Folds the position \p pos into the central image. 
          This routine must be implemented by derived classes

	  \param pos is the position to be folded
      */
      virtual void
      foldThis(real pos[3]) const;

      /** Returns the central image of the position \p pos.

	  \param pos is the position to be folded
	  \return the folded position */
      virtual real*
      fold(const real pos[3]) const;
      
      /** Get a random position within the central simulation box. The
          positions are assigned with each coordinate on [0, boxL]. */
      virtual void
      getRandomPos(real res[3]) const;

    }; 
  }
}

#endif
