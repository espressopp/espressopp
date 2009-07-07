#include "python.hpp"

#include "pairs/Set.hpp"
#include "pairs/Computer.hpp"

using namespace espresso::pairs;

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void
Set::registerPython() {

  using namespace espresso::python;

  void (Set::*foreach_nonconst)(Computer& computer) = &Set::foreach;

  // also register the abstract class Set to make virtual functions available
  // be careful: boost::noncopyable must be used for abstract classes with pure routines
  // no_init must be used as the abstract class Set has no constructor

  class_< Set, boost::noncopyable >("pairs_Set", no_init)
    .def("foreach", foreach_nonconst);
  ;
}


