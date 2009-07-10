#include "python.hpp"

#include "particles/Set.hpp"
#include "particles/Computer.hpp"

using namespace espresso::particles;

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void
Set::registerPython() {
  using namespace espresso::python;

  void (Set::*foreach_nonconst)(Computer& computer) = &Set::foreach;

  class_< Set, boost::noncopyable >
    ("particles_Set", no_init)
    .def("foreach", pure_virtual(foreach_nonconst))
    ;
}

