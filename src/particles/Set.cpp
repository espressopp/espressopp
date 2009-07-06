#include <boost/python.hpp>

#include "particles/Set.hpp"
#include "particles/Computer.hpp"

using namespace espresso::particles;

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void
Set::registerPython() {

  using namespace boost::python;

  void (Set::*foreach_nonconst)(Computer& computer) = &Set::foreach;

  // also register the abstract class Set to make virtual functions available
  // be careful: boost::noncopyable must be used for abstract classes with pure routines
  // no_init must be used as the abstract class Set has no constructor

  class_< Set, boost::noncopyable >("particles_Set", no_init)
    .def("foreach", foreach_nonconst)
    ;
}

