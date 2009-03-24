#ifndef _BASE_REAL3D_HPP
#define _BASE_REAL3D_HPP

#include "estypes.hpp"
#include "esutil/Vector3D.hpp"

namespace espresso {
  namespace base {

    /// make Real3D globally accessible for convenience
    typedef esutil::Vector3D<real> Real3D;

    void registerPythonReal3D();
  }
}

#endif
