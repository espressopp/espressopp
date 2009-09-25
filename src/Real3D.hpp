#ifndef _REAL3D_HPP
#define _REAL3D_HPP

#include "estypes.hpp"
#include "esutil/Vector3D.hpp"
#include "esutil/Box.hpp"
#include "esutil/Grid.hpp"

namespace espresso {
    /// make Real3D globally accessible for convenience
    typedef esutil::Vector3D<real> Real3D;
    /// make Real3DBox globally accessible for convenience
    typedef esutil::Box<Real3D> Real3DBox;
    /// make Real3DGrid globally accessible for convenience
    typedef esutil::Grid<Real3D> Real3DGrid;

    void registerPythonReal3D();
}

#endif
