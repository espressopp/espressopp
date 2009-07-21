#include "particles/Set.hpp"
#include "particles/Computer.hpp"
#include "particles/ParticleHandle.hpp"
#include "particles/Storage.hpp"

#include <boost/bind.hpp>
#include <python.hpp>

// Implementation of a generic isMember

using namespace espresso::particles;

namespace {
  class Found {};
  void findMember(const ParticleHandle p1, 
		  const ParticleHandle p2) {
    if (p1 == p2) throw Found();
  }
}

bool
Set::contains(ParticleHandle p) {
  try {
    foreach(boost::bind(findMember, p, _1));
  } catch (Found) {
    return true;
  }
  return false;
}

bool                                                                           
Set::contains(ParticleId pid) {
  return contains(getStorage()->getParticleHandle(pid));
}

void
Set::foreach(Computer &computer) {
  computer.prepare(getStorage());
  try {
    foreach(boost::bind(&Computer::apply, &computer, _1));
  } catch (ForeachBreak) {
    computer.finalize();
    throw;
  } 
  computer.finalize();
}
  

void
Set::foreach(const Computer::SelfPtr computer)
{ foreach(*computer); }


//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void
Set::registerPython() {
  using namespace espresso::python;
  
  void (Set::*pyForeach)(Computer::SelfPtr computer) 
    = &Set::foreach;

  bool (Set::*pyContains)(ParticleId pid)
    = &Set::contains;

  class_< Set, boost::noncopyable >
    ("particles_Set", no_init)
    .def("foreach", pyForeach)
    .def("__contains__", pyContains)
    ;
}
