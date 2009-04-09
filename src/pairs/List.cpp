#include "List.hpp"
#include <stdexcept>
#include <types.hpp>
#include <algorithm>

using namespace espresso::pairs;
using namespace espresso::particles;
using namespace std;

/* Destructor */

List::~List() 

{}

/* Constructor for this class  */

List::List(bc::BC& _bc, 
           Storage& _storage, 
           boost::shared_ptr< Property<Real3D> > _coordinates) :

   storage(_storage),
   bc(_bc),
   coordinates(_coordinates) 

   {}

size_t List::size() const {
   return id_list.size();
}

bool List::findPair(ParticleId id1, ParticleId id2) const {
   Tuple T (id1, id2);
   std::vector<Tuple>::const_iterator it = find (id_list.begin(), id_list.end(), T);
   return it != id_list.end();
}

void List::addPair(ParticleId id1, ParticleId id2) {
   Tuple T (id1, id2);
   id_list.push_back(Tuple(id1, id2));
}

void List::deletePair(ParticleId id1, ParticleId id2) {  
   Tuple T(id1, id2);
   std::vector<Tuple>::iterator it = find (id_list.begin(), id_list.end(), T);
   if (it == id_list.end()) {
     throw std::runtime_error("deletePair: tuple not found");
   }
   id_list.erase(it);
}

/*--------------------------------------------------------------------------
List::foreach(PairComputer)
-------------------------------------------------------------------------- */

template<class Computer, class Handle>
void List::foreach(Computer& pairComputer) const {
  vector<Tuple>::const_iterator it;
  ConstPropertyHandle<Real3D> pos = *coordinates;

  for (it = id_list.begin(); it != id_list.end(); it++) {
    Handle pref1 = storage.getParticleHandle(it->first);
    Handle pref2 = storage.getParticleHandle(it->second);
    Real3D dist = bc.getDist(pos[pref1], pos[pref2]);
    
    pairComputer(dist, pref1, pref2);
  }
}

void List::foreach(Computer& pairComputer) const {
  foreach<Computer, ParticleHandle>(pairComputer);
}

void List::foreach(ConstComputer& pairComputer) const {
  foreach<ConstComputer, ConstParticleHandle>(pairComputer);
}

