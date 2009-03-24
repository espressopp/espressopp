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
           PropertyId _coordinates) :

   storage(_storage),
   bc(_bc),
   coordinates(_coordinates) 

   {}

size_t List::size() {
   return id_list.size();
}

bool List::findPair(ParticleId id1, ParticleId id2) {
   Tuple T (id1, id2);
   std::vector<Tuple>::iterator it = find (id_list.begin(), id_list.end(), T);
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

void List::foreach(Computer& pairComputer) {

  vector<Tuple>::iterator it;

  PropertyReference<Real3D> pos =
    storage.getPropertyReference<Real3D>(coordinates);

  for (it = id_list.begin(); it != id_list.end(); it++) {
    ParticleReference pref1 = storage.getParticleReference(it->first);
    ParticleReference pref2 = storage.getParticleReference(it->second);

    Real3D dist = bc.getDist(pos[pref1], pos[pref2]);
    
    pairComputer(dist, pref1, pref2);
  }
}


/*--------------------------------------------------------------------------
List::foreach(ConstParticlePairComputer)
-------------------------------------------------------------------------- */

void List::foreach(ConstComputer& pairComputer) const {

  vector<Tuple>::const_iterator it;

  PropertyReference<Real3D> pos =
    storage.getPropertyReference<Real3D>(coordinates);
  
  for (it = id_list.begin(); it != id_list.end(); it++) {
    
    ConstParticleReference pref1 = storage.getParticleReference(it->first);
    ConstParticleReference pref2 = storage.getParticleReference(it->second);
    
    Real3D dist = bc.getDist(pos[pref1], pos[pref2]);
    
    pairComputer(dist, pref1, pref2);
  }
}

