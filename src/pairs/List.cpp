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
//TODO: right now only take _storage1, _posProperty1
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

bool List::foreachApply(Computer &pairComputer) {
// computer.prepare(storage, storage);
// computer.finalize();

  bool cont;
  vector<Tuple>::const_iterator it;
  storage::PropertyHandle<Real3D> pos = posProperty->getHandle(storage);

  for (it = id_list.begin(); it != id_list.end(); it++) {
    ParticleHandle p1 = storage->getParticleHandle(it->first);
    ParticleHandle p2 = storage->getParticleHandle(it->second);
    Real3D dist = bc->getDist(pos[p1], pos[p2]);
    
    cont = pairComputer.apply(dist, p1, p2);
  }
  return cont;
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void
List::registerPython() {
  using namespace espresso::python;

  class_< List, bases< Set > >
    ("pairs_List", "Put pairs::List docstring here.", 
     init< bc::BC::SelfPtr, storage::Storage::SelfPtr, Property< Real3D >::SelfPtr >())
    .def(init< bc::BC::SelfPtr,
	 storage::Storage::SelfPtr,
         storage::Storage::SelfPtr,
         Property< Real3D >::SelfPtr,
         Property< Real3D >::SelfPtr >()) 
    .def("addPair", &List::addPair)
    .def("deletePair", &List::deletePair)
    .def("size", &List::size)
    .def("findPair", &List::findPair)
    ;
}
