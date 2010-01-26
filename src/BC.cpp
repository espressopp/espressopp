#include "BC.hpp"
#include <python.hpp>
#include <cmath>
#include "Real3D.hpp"
#include "Real3DRef.hpp"

namespace espresso {

  Real3D
  BC::getRandomPos() const {
    Real3D res;
    getRandomPos(res);
    return res;
  }

  //////////////////////////////////////////////////
  // REGISTRATION WITH PYTHON
  //////////////////////////////////////////////////
  void
  BC::registerPython() {
    using namespace espresso::python;
    
    // also register the abstract class BC to make virtual functions available
    // be careful: boost::noncopyable must be used for abstract classes with pure routines
    // no_init must be used as the abstract class BC has no constructor
    
    class_<BC, boost::noncopyable >("bc_BC", no_init)
      .def("setBoxL", &BC::setBoxL)
      ;
  }
}
