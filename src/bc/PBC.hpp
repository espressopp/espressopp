#ifndef _BC_PBC_HPP
#define _BC_PBC_HPP

#include <cmath>
#include <logging.hpp>
#include <bc/BC.hpp>

namespace espresso {
  namespace bc {
    /** Class for periodic boundary conditions in all three dimensions. */
    //currently we assume a cube
    //names like PBC will need to change to cubic or orthorhombic

    class PBC : public BC {

    private:
      real length;
      real lengthInverse;

      static LOG4ESPP_DECL_LOGGER(theLogger);

    public:

      static void registerPython();

      /** Constructor for cubic box */
      PBC();
      PBC(real _length);
      /** Destructor for periodic boundary conditions */
      virtual ~PBC();

      /** Method to set the length of the side of the cubic simulation cell */
      virtual void set(real _length);
      
      /** Method that returns the length of the side of the cubic simulation cell */
      virtual real getLength(void) const;

      //PMI and Python Visible
      /** Method to compute the minimum image distance */
      virtual Real3D getDist(const Real3D& pos1, const Real3D& pos2) const;

      /** Method to get a random position within the central simulation box. The
          positions are assigned with each coordinate on [0, length] */
      virtual Real3D randomPos(void);
      
    };
  }
}

#endif
