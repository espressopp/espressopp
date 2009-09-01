#include "python.hpp"
#include "boost/bind.hpp"

#include "pairs/Set.hpp"
#include "pairs/Computer.hpp"

using namespace espresso::pairs;

bool 
Set::foreachPair(Computer &computer) {
  computer.prepare(getLeftStorage(), getRightStorage());
  bool res = foreachPairApply(computer);
  computer.finalize();
  return res;
}

bool
Set::foreachPair(Computer::SelfPtr computer) {
  return foreachPair(*computer);
}

bool
Set::enclForeachPairWithin(Computer &computer) {
  return foreachPair(computer);
}


//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void
Set::registerPython() {
  using namespace espresso::python;

  bool (Set::*pyForeach)(Computer::SelfPtr computer) 
    = &Set::foreachPair;

  class_< Set, boost::noncopyable >("pairs_Set", no_init)
    .def("foreach", pyForeach)
    ;

}


