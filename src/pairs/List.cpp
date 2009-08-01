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
           Storage::SelfPtr _storage, 
           Property< Real3D >::SelfPtr _posProperty) :

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

/*--------------------------------------------------------------------------
List::foreach(PairComputer)
-------------------------------------------------------------------------- */

template<class Computer, class Handle>
void List::foreach(Computer& pairComputer) const {
  vector<Tuple>::const_iterator it;
  PropertyHandle<Real3D> pos = *posProperty;

  for (it = id_list.begin(); it != id_list.end(); it++) {
    Handle pref1 = storage->getParticleHandle(it->first);
    Handle pref2 = storage->getParticleHandle(it->second);
    Real3D dist = bc->getDist(pos[pref1], pos[pref2]);
    
    pairComputer(dist, pref1, pref2);
  }
}

void List::foreach(Computer& pairComputer) {
  foreach<Computer, ParticleHandle>(pairComputer);
}

/**
void List::foreach(ConstComputer& pairComputer) const {
  foreach<ConstComputer, ConstParticleHandle>(pairComputer);
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
