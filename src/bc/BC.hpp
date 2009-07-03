#ifndef _BC_BC_HPP
#define _BC_BC_HPP

#include <boost/shared_ptr.hpp>
#include "types.hpp"

namespace espresso {
  namespace bc {
    
    /** Abstract class for boundary condtions. */
    class BC {
    public:
      typedef boost::shared_ptr< BC > SelfPtr;

      /** Virtual destructor for boundary conditions. */
      virtual ~BC() {}

      /** Delivers the minimum image distance vector between two
          positions. This routine must be implemented by derived
          classes.

          \param pos1, pos2 are the two positions 
          \returns the distance vector (pos2 - pos1)
      */
      virtual Real3D 
      getDist(const Real3D& pos1, 
	      const Real3D& pos2) const = 0;

      
      /** Folds the position \p pos into the central image. 
          This routine must be implemented by derived classes

	  \param pos is the position to be folded
      */
      virtual void
      foldThis(Real3D& pos) const = 0;

      /** Returns the central image of the position \p pos .

	  \param pos is the position to be folded
	  \return the folded position */
      virtual Real3D
      fold(const Real3D& pos) const;

    public:

      static void registerPython();

    };
  }
}

#endif
