#include "python.hpp"
#include "boost/bind.hpp"

#include "pairs/Set.hpp"
#include "pairs/Computer.hpp"

using namespace espresso::pairs;
#include <iostream>
Set::Set(particles::Storage::SelfPtr _storage1,
	 particles::Storage::SelfPtr _storage2)
  : storage1(_storage1), storage2(_storage2) {}

void 
Set::foreach(Computer &computer) {
  computer.prepare(storage1, storage2);
  try {
    foreachApply(computer);
  } catch (ForeachBreak) {
    computer.finalize();
    throw;
  }
  computer.finalize();
}

void
Set::foreach(Computer::SelfPtr computer) {
  foreach(*computer);
}


//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void
Set::registerPython() {
  using namespace espresso::python;

  void (Set::*foreach1)(Computer::SelfPtr computer) 
    = &Set::foreach;

  // also register the abstract class Set to make virtual functions available
  // be careful: boost::noncopyable must be used for abstract classes with pure routines
  // no_init must be used as the abstract class Set has no constructor

  class_< Set, boost::noncopyable >("pairs_Set", no_init)
    .def("foreach", foreach1)
    .def("bump", &Set::bump)
    ;

}


