#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>

#include "BC.hpp"

using namespace espresso;
using namespace espresso::bc;

Real3D 
BC::fold(const Real3D &pos) const {
  Real3D res(pos);
  foldThis(res);
  return res;
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void
BC::registerPython() {
  using namespace boost::python;

  // also register the abstract class BC to make virtual functions available
  // be careful: boost::noncopyable must be used for abstract classes with pure routines
  // no_init must be used as the abstract class BC has no constructor

  class_<BC, boost::shared_ptr<BC>, boost::noncopyable >("bc_BC", no_init)
    .def("foldThis", &BC::foldThis)
    .def("fold", &BC::fold)
    .def("getDist", &BC::getDist)
    ;
}

