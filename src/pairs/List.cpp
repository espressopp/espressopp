#include "types.hpp"
#include "List.hpp"
#include "python.hpp"

#include <stdexcept>
#include <algorithm>

using namespace espresso::pairs;
using namespace espresso::storage;
using namespace std;

/* Destructor */
List::~List()
{}

/* Constructor */
List::List(bc::BC::SelfPtr _bc, 
           Storage::SelfPtr _storage1, 
           Storage::SelfPtr _storage2) 
  : bc(_bc),
    storage1(_storage1),
    storage2(_storage2)
    {}

/* Constructor */
List::List(bc::BC::SelfPtr _bc,
           Storage::SelfPtr _storage)
  : bc(_bc),
    storage1(_storage),
    storage2(_storage)
    {}

size_t List::size() const {
   return id_list.size();
}

bool List::findPair(ParticleId id1, ParticleId id2) const {
   Tuple T(id1, id2);
   vector<Tuple>::const_iterator it = find(id_list.begin(), id_list.end(), T);
   return it != id_list.end();
}

void List::addPair(ParticleId id1, ParticleId id2) {
   Tuple T(id1, id2);
   //TODO: include id's in error messages
   //TODO: next line is throwing confusing message
   if(id1 < 0 || id2 < 0) {
      throw runtime_error("pairs::List::addPair: ParticleId less than zero.");
   }
   //TODO: next test should be only if one storage is used
   if(id1 == id2) {
      throw runtime_error("pairs::List::addPair: Particles id's are equal.");
   }
   //TODO: check if id1/2 d2 is greater than number of particles in storage1/2
   //TODO: if two storages make sure id1 is from storage1 and id2 is from storage2
   id_list.push_back(Tuple(id1, id2));
}

void List::deletePair(ParticleId id1, ParticleId id2) {
   Tuple T(id1, id2);
   vector<Tuple>::iterator it = find(id_list.begin(), id_list.end(), T);
   if (it == id_list.end()) {
     throw runtime_error("pairs::List::deletePair: Tuple not found.");
     //TODO: write id1 and id2 in error message
   }
   id_list.erase(it);
}

bool List::foreachPairApply(Computer &computer) {
   vector<Tuple>::const_iterator it;
   storage::PropertyHandle<Real3D> pos1 = storage1->getPositionPropertyHandle();
   storage::PropertyHandle<Real3D> pos2 = storage2->getPositionPropertyHandle();

   for(it = id_list.begin(); it != id_list.end(); it++) {
     ParticleHandle p1 = storage1->getParticleHandle(it->first);
     ParticleHandle p2 = storage2->getParticleHandle(it->second);
     Real3D dist = bc->getDist(pos1[p1], pos2[p2]);
    
     if(!computer.apply(dist, p1, p2)) return false;
   }

   return true;
}

Storage::SelfPtr List::getLeftStorage()
{ return storage1; }

Storage::SelfPtr List::getRightStorage()
{ return storage2; }


//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void
List::registerPython() {
  using namespace espresso::python;

  class_< List, bases< Set > >
    ("pairs_List", "Put pairs::List docstring here.", 
     init< bc::BC::SelfPtr, storage::Storage::SelfPtr >())
    .def(init< bc::BC::SelfPtr,
	 storage::Storage::SelfPtr,
         storage::Storage::SelfPtr >()) 
    .def("addPair", &List::addPair)
    .def("deletePair", &List::deletePair)
    .def("size", &List::size)
    .def("findPair", &List::findPair)
    ;
}
