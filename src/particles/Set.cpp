#include "particles/Set.hpp"
#include "particles/Computer.hpp"
#include "storage/ParticleHandle.hpp"
#include "storage/Storage.hpp"

#include <boost/bind.hpp>
#include <python.hpp>

using namespace espresso::particles;
using namespace espresso::storage;

namespace {

  class Found {};
  struct FindMember: public Computer {
    ParticleHandle p1;

    FindMember(ParticleHandle _p1) : p1(_p1) 
    {}

    void prepare(const Storage::SelfPtr storage) {}

    void apply(const ParticleHandle p2) {
      if (p1 == p2) throw Found();
    }
  };
}

bool
Set::contains(ParticleHandle p) {
  FindMember findMember(p);
  try {
    foreach(findMember);
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
    foreachApply(computer);
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
