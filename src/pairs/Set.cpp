#include "python.hpp"
#include "boost/bind.hpp"

#include "pairs/Set.hpp"
#include "pairs/Computer.hpp"

using namespace espresso::pairs;

Set::Set(storage::Storage::SelfPtr _storage1,
	 storage::Storage::SelfPtr _storage2)
  : storage1(_storage1), storage2(_storage2) {}

bool 
Set::foreach(Computer &computer) {
  computer.prepare(storage1, storage2);
  bool res = foreachApply(computer);
  computer.finalize();
  return res;
}

bool
Set::foreach(Computer::SelfPtr computer) {
  return foreach(*computer);
}


//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void
Set::registerPython() {
  using namespace espresso::python;

  bool (Set::*pyForeach)(Computer::SelfPtr computer) 
    = &Set::foreach;

  class_< Set, boost::noncopyable >("pairs_Set", no_init)
    .def("foreach", pyForeach)
    ;

}


