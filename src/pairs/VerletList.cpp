#include "VerletList.hpp"
#include "python.hpp"
#include "types.hpp"
#include "estypes.hpp"

#include <stdexcept>
#include <algorithm>

// remove next line
#include <iostream>


using namespace espresso;
using namespace espresso::pairs;
using namespace espresso::storage;
using namespace std;

/* Destructor */
VerletList::~VerletList() {}

/* Constructor */
VerletList::VerletList(bc::BC::SelfPtr _bc,
		       Storage::SelfPtr _storage,
                       real _radius,
		       real _skin)
  : bc(_bc),
    storage1(_storage),
    storage2(_storage),
    radius(_radius),
    skin(_skin)
    {}

/* Constructor */
VerletList::VerletList(bc::BC::SelfPtr _bc, 
		       Storage::SelfPtr _storage1, 
		       Storage::SelfPtr _storage2,
                       real _radius, 
		       real _skin) 
  : bc(_bc),
    storage1(_storage1),
    storage2(_storage2),
    radius(_radius),
    skin(_skin)
    {}

namespace {
  typedef std::pair< storage::ParticleHandle, storage::ParticleHandle > phTuple;

  struct DistanceComputer: public pairs::Computer {
    DistanceComputer(storage::Storage::SelfPtr _stor,
                     bc::BC::SelfPtr _bc,
                     espresso::real _radius,
                     espresso::real _skin,
                     std::vector< phTuple > &_ph_list):
                     stor(_stor), bc(_bc), radius(_radius),
                     skin(_skin), ph_list(_ph_list)
                     {}

    void prepare(Storage::SelfPtr s1, Storage::SelfPtr s2) {
      ph_list.clear();
    }

    bool apply(const Real3D &d, const ParticleHandle p1, const ParticleHandle p2) {
      
      storage::PropertyHandle<Real3D> pos1 = stor->getPositionPropertyHandle();
      storage::PropertyHandle<Real3D> pos2 = stor->getPositionPropertyHandle();

      // replace with radiusPlusSkinSqr
      //Real3D d = bc->getDist(pos1[p1], pos2[p2]);
      if(d.sqr() <= pow(radius + skin, 2)){
        ph_list.push_back(phTuple(p1, p2));
      }
      return true;
    }
   
    void finalize() {
    }

    // what if have two storages?
    storage::Storage::SelfPtr stor;
    bc::BC::SelfPtr bc;
    espresso::real radius;
    espresso::real skin; // probably can delete
    std::vector< phTuple > &ph_list;
  };
}

void VerletList::update() {
  // how to handle both storages?
  DistanceComputer distanceComputer(storage1, bc, radius, skin, ph_list);
  storage1->enclForeachPairWithin(distanceComputer);
}

bool VerletList::foreachPairApply(Computer &computer) {
   vector<phTuple>::const_iterator it;
   storage::PropertyHandle<Real3D> pos1 = storage1->getPositionPropertyHandle();
   storage::PropertyHandle<Real3D> pos2 = storage2->getPositionPropertyHandle();

   update();
   std::cout << "VERLET LIST SIZE = " << ph_list.size() << std::endl;

   for(it = ph_list.begin(); it != ph_list.end(); it++) {
     ParticleHandle p1 = it->first;
     ParticleHandle p2 = it->second;
     if (!p1 || !p2) {
       throw runtime_error("VerletList::foreachPairApply: pair in Verlet list does not exist in storages.");
     }
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
     espresso::real,
     espresso::real >())
    .def(init<bc::BC::SelfPtr,
	 storage::Storage::SelfPtr,
         storage::Storage::SelfPtr,
         espresso::real,
	 real >())
     .add_property("radius", &VerletList::getRadius, &VerletList::setRadius)
     .add_property("skin", &VerletList::getSkin, &VerletList::setSkin)
    ;
}
