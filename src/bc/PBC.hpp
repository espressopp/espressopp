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

      virtual void set(real _length);

      //PMI and Python Visible
      /** Method to compute the minimum image distance */
      virtual Real3D getDist(const Real3D& pos1, const Real3D& pos2) const;
      
    };
  }
}

#endif
