#ifndef _BC_BC_HPP
#define _BC_BC_HPP

#include "types.hpp"

namespace espresso {
  namespace bc {
    
    /** Abstract class for boundary condtions. */
    class BC {
    public:
      typedef shared_ptr< BC > SelfPtr;

      /** Virtual destructor for boundary conditions. */
      virtual ~BC() {}

      /** Delivers the minimum image distance vector between two
          positions. This routine must be implemented by derived
          classes.

          \param dist is the distance vector (pos2 - pos1)
          \param distSqr is the square of the magnitude
                 of the separation vector
          \param pos1, pos2 are the particle positions 
      */
      virtual void
      getMinimumImageVector(real dist[3],
                            real &distSqr,
                            const real pos1[3],
                            const real pos2[3]) const = 0;

      /** Folds the position \p pos into the central image. 
          This routine must be implemented by derived classes

	  \param pos is the position to be folded
      */
      virtual void
      foldThis(real pos[3]) const = 0;

      /** Returns the central image of the position \p pos .

	  \param pos is the position to be folded
	  \return the folded position */
      virtual void
      fold(const real pos[3]) const;

    };
  }
}

#endif
