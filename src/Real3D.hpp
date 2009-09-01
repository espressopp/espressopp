#ifndef _REAL3D_HPP
#define _REAL3D_HPP

#include "estypes.hpp"
#include "esutil/Vector3D.hpp"
#include "esutil/Box.hpp"

namespace espresso {
    /// make Real3D globally accessible for convenience
    typedef esutil::Vector3D<real> Real3D;
    /// make Box globally accessible for convenience
    typedef esutil::Box<real> RealBox;

    void registerPythonReal3D();
}

#endif
