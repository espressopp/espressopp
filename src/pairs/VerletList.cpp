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
		       real _skin)
  : bc(_bc),
    storage1(_storage),
    storage2(_storage),
    skin(_skin)
    {}

/* Constructor */
VerletList::VerletList(bc::BC::SelfPtr _bc, 
		       Storage::SelfPtr _storage1, 
		       Storage::SelfPtr _storage2, 
		       real _skin) 
  : bc(_bc),
    storage1(_storage1),
    storage2(_storage2),
    skin(_skin)
    {}

void VerletList::update() {}

espresso::real VerletList::getSkin() { return skin; }

void VerletList::setSkin(espresso::real _skin) { skin = _skin; }

bool VerletList::foreachPairApply(Computer &computer) {
   vector<Tuple>::const_iterator it;
   storage::PropertyHandle<Real3D> pos1 = getLeftStorage()->getPositionPropertyHandle();
   storage::PropertyHandle<Real3D> pos2 = getRightStorage()->getPositionPropertyHandle();

   for(it = id_list.begin(); it != id_list.end(); it++) {
     ParticleHandle p1 = storage1->getParticleHandle(it->first);
     ParticleHandle p2 = storage2->getParticleHandle(it->second);
     Real3D dist = bc->getDist(pos1[p1], pos2[p2]);
    
     if(!computer.apply(dist, p1, p2)) return false;
   }
   return true;
}

Storage::SelfPtr VerletList::getLeftStorage()
{ return storage1; }

Storage::SelfPtr VerletList::getRightStorage()
{ return storage2; }

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void
VerletList::registerPython() {
  using namespace espresso::python;

  class_< VerletList, bases< Set > >
    ("pairs_VerletList", "Put pairs::VerletList docstring here.", 
     init< bc::BC::SelfPtr,
     storage::Storage::SelfPtr,
     espresso::real >())
    .def(init<bc::BC::SelfPtr,
	 storage::Storage::SelfPtr,
         storage::Storage::SelfPtr,
	 real >())
    ;
}
