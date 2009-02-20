#include "List.hpp"
#include <stdexcept>
#include <types.hpp>

using namespace espresso::pairs;
using namespace espresso::particlestorage;
using namespace std;

/* Destructor */

List::~List() 

{}

/* Constructor for this class  */

List::List(espresso::bc::BC& _bc, 
           ParticleStorage& _storage, 
           size_t _coordinates) :

   storage(_storage),
   bc(_bc),
   coordinates(_coordinates) 

   {}

int List::size() {
   return id_list.size();
}

bool List::findPair(size_t id1, size_t id2) {

   Tuple T (id1, id2);

   std::vector<Tuple>::iterator it = find (id_list.begin(), id_list.end(), T);

   return it != id_list.end();

}

void List::addPair(size_t id1, size_t id2) 

{
   Tuple T (id1, id2);

   id_list.push_back(Tuple(id1, id2));
}

void List::deletePair(size_t id1, size_t id2) 

{  
   Tuple T(id1, id2);

   std::vector<Tuple>::iterator it = find (id_list.begin(), id_list.end(), T);

   if (it == id_list.end()) {

      throw std::runtime_error("deletePair: tuple not found");

   }

   id_list.erase(it);

}

void List::foreach(ParticlePairComputer& pairComputer) {

   vector<Tuple>::iterator it;

   ParticleStorage::PropertyTraits<Real3D>::Reference 
       pos = storage.getProperty<Real3D>(coordinates);

   for (it = id_list.begin(); it != id_list.end(); it++) {

       ParticleStorage::reference pref1 = storage.getParticleByID(it->first);
       ParticleStorage::reference pref2 = storage.getParticleByID(it->second);

       Real3D dist = bc.getDist(pos[pref1], pos[pref2]);

       pairComputer(dist, pref1, pref2);

   }
}


/*--------------------------------------------------------------------------
List::foreach(ConstParticlePairComputer)
-------------------------------------------------------------------------- */

void List::foreach(ConstParticlePairComputer& pairComputer) const {

   vector<Tuple>::const_iterator it;

   ParticleStorage::PropertyTraits<Real3D>::Reference 
       pos = storage.getProperty<Real3D>(coordinates);

   for (it = id_list.begin(); it != id_list.end(); it++) {

       ParticleStorage::const_reference pref1 = storage.getParticleByID(it->first);
       ParticleStorage::const_reference pref2 = storage.getParticleByID(it->second);

       Real3D dist = bc.getDist(pos[pref1], pos[pref2]);

       pairComputer(dist, pref1, pref2);
   }
}

