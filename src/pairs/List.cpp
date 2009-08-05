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
//TODO: is it right to take _storage1, _posProperty1
List::List(bc::BC::SelfPtr _bc, 
           Storage::SelfPtr _storage1, 
           Storage::SelfPtr _storage2, 
           Property< Real3D >::SelfPtr _posProperty1,
           Property< Real3D >::SelfPtr _posProperty2) 
  : Set(_storage1, _storage2),
   storage(_storage1),
   bc(_bc),
   posProperty(_posProperty1)
   {}

List::List(bc::BC::SelfPtr _bc,
           Storage::SelfPtr _storage,
           Property< Real3D >::SelfPtr _posProperty)
  : Set(_storage, _storage),
    storage(_storage),
    bc(_bc),
    posProperty(_posProperty)
    {}

size_t List::size() const {
   return id_list.size();
}

bool List::findPair(ParticleId id1, ParticleId id2) const {
   Tuple T (id1, id2);
   vector<Tuple>::const_iterator it = find (id_list.begin(), id_list.end(), T);
   return it != id_list.end();
}

void List::addPair(ParticleId id1, ParticleId id2) {
   Tuple T (id1, id2);
   id_list.push_back(Tuple(id1, id2));
}

void List::deletePair(ParticleId id1, ParticleId id2) {  
   Tuple T(id1, id2);
   vector<Tuple>::iterator it = find (id_list.begin(), id_list.end(), T);
   if (it == id_list.end()) {
     throw runtime_error("deletePair: tuple not found");
   }
   id_list.erase(it);
}

void List::foreachApply(Computer &computer) {
  /*computer.prepare(storage, storage);
  storage->foreach(computer);
  computer.finalize();*/
}
/*
void List::foreach(Computer& pairComputer) {
  vector<Tuple>::const_iterator it;
  PropertyHandle<Real3D> pos = *posProperty;

  for (it = id_list.begin(); it != id_list.end(); it++) {
    Handle pref1 = storage->getParticleHandle(it->first);
    Handle pref2 = storage->getParticleHandle(it->second);
    Real3D dist = bc->getDist(pos[pref1], pos[pref2]);
    
    pairComputer(dist, pref1, pref2);
  }
}
*/


//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void
List::registerPython() {
  using namespace espresso::python;

  class_< List, bases< Set > >
    ("pairs_List", 
     init< bc::BC::SelfPtr, storage::Storage::SelfPtr, Property< Real3D >::SelfPtr >())
    .def("addPair", &List::addPair)
    ;
}
