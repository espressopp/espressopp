#include "All.hpp"
#include "esutil/virtual_functional.hpp"

#include <cstdio>

using namespace espresso;
using namespace espresso::pairs;
using namespace espresso::particlestorage;

// some abbreviations

typedef ParticleStorage::PropertyTraits<size_t>::ConstReference SizeRef;

typedef ParticleStorage::PropertyTraits<Real3D>::ConstReference RealArrayRef;


// Helper class 1

template<class Reference>
class Traverser1 : public esutil::VirtualUnaryFunction<Reference, void>  {

   typedef ParticlePairComputerBase<Reference> PairComputer;

   private:

    class Traverser2 : public esutil::VirtualUnaryFunction<Reference, void> {
   
      public:

      const espresso::bc::BC& bc;

      const SizeRef id;
      const RealArrayRef pos;

      const Reference pref1;
      const Real3D pos1;
      const size_t id1;

      PairComputer& pairComputer;

      Traverser2(const All* all,
                 const Reference pref,
                 PairComputer& _pairComputer
	  ) :

        bc(all->getBC()), 
        id(all->getSet().getStorage()->getIDProperty()),
        pos(all->getSet().getStorage()->template getProperty<Real3D>(all->getCoordinateProperty())),
        pref1(pref),
        pos1(pos[pref1]),
        id1(id[pref1]),
        pairComputer(_pairComputer)

        {  } 

      virtual void operator()(const Reference pref2) {
   
        if (id1 < id[pref2]) {

          Real3D dist = bc.getDist(pos1, pos[pref2]);
   
          pairComputer(dist, pref1, pref2);
        }
      }

     };

  public:

     PairComputer& pairComputer;

     const All* all;

     Traverser1(const All* _all, PairComputer& _pairComputer) :

       pairComputer(_pairComputer),
       all(_all)

      {  }

      virtual void operator()(const Reference pref) {

         // printf ("Traverser1: will call traverser2\n");

         Traverser2 traverser2(all, pref, pairComputer);

         all->getSet().foreach(traverser2);
      }
};

/*--------------------------------------------------------------------------
All::~All()
-------------------------------------------------------------------------- */

All::~All() {}

/*--------------------------------------------------------------------------
All::All(boundary_conditions, particle_set)
-------------------------------------------------------------------------- */

All::All(espresso::bc::BC& _bc, espresso::particleset::ParticleSet& _set, size_t _coordinates):

  set(_set),
  bc(_bc),
  coordinates(_coordinates)

{ }

/*--------------------------------------------------------------------------
All::foreach(ParticlePairComputer)
-------------------------------------------------------------------------- */

void All::foreach(ParticlePairComputer& pairComputer) {

   // printf ("ParticlePairComputer non-const\n");

   Traverser1<ParticleStorage::reference> traverser1(this, pairComputer);;

   set.foreach(traverser1);

}
       
/*--------------------------------------------------------------------------
All::foreach(ConstParticlePairComputer)
-------------------------------------------------------------------------- */

void All::foreach(ConstParticlePairComputer& pairComputer) const {

   // printf ("ParticlePairComputer const\n");

   Traverser1<ParticleStorage::const_reference> traverser1(this, pairComputer);;

   set.foreach(traverser1);
}

