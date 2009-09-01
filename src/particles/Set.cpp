#include "particles/Set.hpp"
#include "particles/Computer.hpp"
#include "storage/ParticleHandle.hpp"
#include "storage/Storage.hpp"

#include <boost/bind.hpp>
#include <python.hpp>

using namespace espresso::particles;
using namespace espresso::storage;

namespace {

  struct FindMember: public Computer {
    ParticleHandle p1;

    FindMember(ParticleHandle _p1) : p1(_p1) 
    {}

    void prepare(const Storage::SelfPtr storage) {}

    bool apply(const ParticleHandle p2) {
      return p1 != p2;
    }
  };
}

bool
Set::contains(ParticleHandle p) {
  FindMember findMember(p);
  return !foreach(findMember);
}

bool                                                                           
Set::contains(ParticleId pid) {
  return contains(getStorage()->getParticleHandle(pid));
}

bool
Set::foreach(Computer &computer) {
  computer.prepare(getStorage());
  bool res = foreachApply(computer);
  computer.finalize();
  return res;
}
  
bool
Set::foreach(const Computer::SelfPtr computer)
{ return foreach(*computer); }

bool
Set::enclForeachIn(Computer &computer, const RealBox &)
{ return foreach(computer); }


//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void
Set::registerPython() {
  using namespace espresso::python;
  
  bool (Set::*pyForeach)(Computer::SelfPtr computer) 
    = &Set::foreach;

  bool (Set::*pyContains)(ParticleId pid)
    = &Set::contains;

  class_< Set, boost::noncopyable >
    ("particles_Set", no_init)
    .def("foreach", pyForeach)
    .def("__contains__", pyContains)
    ;
}
