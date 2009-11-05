#ifndef _BC_PERIODICBC_HPP
#define _BC_PERIODICBC_HPP

#include "logging.hpp"
#include "BC.hpp"

namespace espresso {
  namespace bc {
    /** Class for periodic boundary conditions in all three dimensions. */

    class PeriodicBC : public BC {
    public:
      typedef shared_ptr< PeriodicBC > SelfPtr;

    private:
      real length[3];
      real lengthInverse[3];

      static LOG4ESPP_DECL_LOGGER(theLogger);

    public:

      /** Constructor for cubic box */
      PeriodicBC(real _length[3]);
      /** Destructor for periodic boundary conditions */
      virtual ~PeriodicBC();

      /** Method to set the length of the side of the cubic simulation cell */
      virtual void setLength(real _length[3]);
      
      /** Method that returns the length of the side of the cubic simulation cell */
      virtual void getLength(real &length[3]) const;

      // PMI and Python Visible
      /** Fold the position \p pos into the central image. */
      virtual void foldThis(real &pos[3]) const;

      /** Compute the minimum image distance (pos2 - pos1) */
      virtual void getMinimumImageVector(real dist[3],
                                         real &distSqr,
                                         const real pos1[3],
                                         const real pos2[3]) const;

      /** Get a random position within the central simulation box. The
          positions are assigned with each coordinate on [0, length]. */
      virtual void getRandomPos(real &pos[3]);
      
    };
  }
}

#endif
