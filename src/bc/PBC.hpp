#ifndef _BC_PBC_HPP
#define _BC_PBC_HPP

#include <cmath>
#include <logging.hpp>
#include <bc/BC.hpp>

namespace espresso {
  namespace bc {
    /** Class for periodic boundary conditions in all three dimensions. */

    class PBC : public BC {

    private:
      Real3D length;
      Real3D lengthInverse;

      static LOG4ESPP_DECL_LOGGER(theLogger);

    public:

      static void registerPython();

      /** Constructor for cubic box */
      PBC();
      PBC(Real3D _length);
      /** Destructor for periodic boundary conditions */
      virtual ~PBC();

      /** Method to set the length of the side of the cubic simulation cell */
      virtual void set(Real3D length);
      
      /** Method that returns the length of the side of the cubic simulation cell */
      virtual Real3D getLength(void) const;

      // PMI and Python Visible
      /** Fold the position \p pos into the central image. */
      virtual void foldThis(Real3D& pos) const;

      /** Fold the position \pos into the central image. 
       \return the folded position */
      virtual Real3D fold(const Real3D& pos) const;

      /** Compute the minimum image distance (pos2 - pos1) */
      virtual Real3D getDist(const Real3D& pos1, const Real3D& pos2) const;

      /** Get a random position within the central simulation box. The
          positions are assigned with each coordinate on [0, length] */
      virtual Real3D randomPos(void);
      
    };
  }
}

#endif
