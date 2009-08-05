#include "All.hpp"

#include <exception>
#include <cmath>

#include "particles/Computer.hpp"
#include "python.hpp"

using namespace espresso;
using namespace espresso::bc;
using namespace espresso::pairs;
//using namespace espresso::particles;

using namespace boost;

All::All(bc::BC::SelfPtr _bc,
         particles::Set::SelfPtr _set,
         Property< Real3D >::SelfPtr _posProperty ) 
  : Set(_set->getStorage(), _set->getStorage()),
    set(_set), bc(_bc), posProperty(_posProperty) {}

namespace {
  class SameId {};

  struct Traverser2 : particles::Computer {
    bool cont;
    pairs::Computer &computer;
    bc::BC &bc;
    Property< Real3D > &posProperty;

    storage::ParticleHandle p1;
    storage::PropertyHandle< Real3D > pos;

    Traverser2(pairs::Computer &_computer, 
	       bc::BC &_bc, 
	       Property< Real3D > &_posProperty) 
      : cont(true), 
	computer(_computer), bc(_bc), posProperty(_posProperty)
    {}

    void setP1(storage::ParticleHandle _p1) { p1 = _p1; }

    void prepare(storage::Storage::SelfPtr storage) {
      pos = posProperty.getHandle(storage);
    }

    bool apply(storage::ParticleHandle p2) {
      if (p1 == p2) return false;
      Real3D pos1 = pos[p1];
      Real3D pos2 = pos[p2];
      Real3D dist = bc.getDist(pos1, pos2);

//       const Real3D length(21.0);
//       const Real3D length_inv(1.0/length[0], 1.0/length[1], 1.0/length[2]);
//       Real3D dist = pos1 - pos2;
//       dist[0] -= round(dist[0] * length_inv[0]) * length[0];
//       dist[1] -= round(dist[1] * length_inv[1]) * length[1];
//       dist[2] -= round(dist[2] * length_inv[2]) * length[2];
      
      cont = computer.apply(dist, p1, p2);
      return cont;
    }
  };

  struct Traverser1 : particles::Computer {
    Traverser2 traverser2;
    particles::Set &set;

    Traverser1(Traverser2 &_traverser2, particles::Set &_set) 
      : traverser2(_traverser2), set(_set)
    {}

    void prepare(storage::Storage::SelfPtr set) {}

    bool apply(storage::ParticleHandle p1) {
      traverser2.setP1(p1);
      if (!set.foreach(traverser2) && !traverser2.cont) 
	return false;
      return true;
    }
  };

}

bool All::foreachApply(Computer &computer) {
  computer.prepare(set->getStorage(), set->getStorage());
  Traverser2 traverser2(computer, *bc, *posProperty);
  Traverser1 traverser1(traverser2, *set);
  bool res = set->foreach(traverser1);
  computer.finalize();
  return res;
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void 
All::registerPython() {
  using namespace espresso::python;
  class_< All, bases< Set > >
    ("pairs_All", 
     init< 
     bc::BC::SelfPtr, 
     particles::Set::SelfPtr, 
     Property< Real3D >::SelfPtr 
     >())
    ;
}
