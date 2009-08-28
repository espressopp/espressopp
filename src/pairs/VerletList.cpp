#include "VerletList.hpp"
#include "python.hpp"

#include <stdexcept>
#include <algorithm>

using namespace espresso::pairs;
using namespace espresso::storage;
using namespace std;

/* Destructor */
VerletList::~VerletList() {}

/* Constructor */
VerletList::VerletList(bc::BC::SelfPtr _bc,
           Storage::SelfPtr _storage,
	   Property< Real3D >::SelfPtr _posProperty,
		       espresso::real _skin)
  : Set(_storage, _storage),
    bc(_bc),
    posProperty1(_posProperty),
    posProperty2(_posProperty),
    skin(_skin)
    {}

/* Constructor */
VerletList::VerletList(bc::BC::SelfPtr _bc, 
           Storage::SelfPtr _storage1, 
           Storage::SelfPtr _storage2, 
           Property< Real3D >::SelfPtr _posProperty1,
	   Property< Real3D >::SelfPtr _posProperty2,
	   espresso::real _skin) 
  : Set(_storage1, _storage2),
    bc(_bc),
    posProperty1(_posProperty1),
    posProperty2(_posProperty2),
    skin(_skin)
    {}

void VerletList::update() {}

espresso::real VerletList::getSkin() { return skin; }

void VerletList::setSkin(espresso::real _skin) { skin = _skin; }

size_t VerletList::size() const {
   return id_list.size();
}

bool VerletList::foreachApply(Computer &computer) {
   vector<Tuple>::const_iterator it;
   storage::PropertyHandle<Real3D> pos1 = posProperty1->getHandle(storage1);
   storage::PropertyHandle<Real3D> pos2 = posProperty2->getHandle(storage2);

   for(it = id_list.begin(); it != id_list.end(); it++) {
     ParticleHandle p1 = storage1->getParticleHandle(it->first);
     ParticleHandle p2 = storage2->getParticleHandle(it->second);
     Real3D dist = bc->getDist(pos1[p1], pos2[p2]);
    
     if(!computer.apply(dist, p1, p2)) return false;
   }

   return true;
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void
VerletList::registerPython() {
  using namespace espresso::python;

  class_< VerletList, bases< Set > >
    ("pairs_VerletList", "Put pairs::VerletList docstring here.", 
     init< bc::BC::SelfPtr, storage::Storage::SelfPtr, Property< Real3D >::SelfPtr, espresso::real >())
    .def(init< bc::BC::SelfPtr,
	 storage::Storage::SelfPtr,
         storage::Storage::SelfPtr,
         Property< Real3D >::SelfPtr,
         Property< Real3D >::SelfPtr,
         espresso::real >())
    .def("size", &VerletList::size)
    ;
}
