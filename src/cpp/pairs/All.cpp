
#include "All.hpp"
#include "virtual_functional.hpp"

#include <cstdio>

using namespace espresso;
using namespace espresso::pairs;
using namespace espresso::particlestorage;

// some abbreviations

typedef ParticleStorage::PropertyReference<size_t> SizeRef;

typedef ParticleStorage::ArrayPropertyReference<real> RealArrayRef;


// Helper class 1

template<class Reference>
class Traverser1 : public VirtualUnaryFunction<void, Reference>  {

   typedef ParticlePairComputerBase<Reference> PairComputer;

   private:

      class Traverser2 : public VirtualUnaryFunction<void, Reference> {
   
      public:

      espresso::bc::BC& bc;

      SizeRef id;
      RealArrayRef pos;

      const Reference pref1;
   
      PairComputer& pairComputer;

      Traverser2(const All* all,
                 const Reference pref,
                 PairComputer& _pairComputer) :

        bc(all->bc), 
        id(all->set.getStorage()->getIDProperty()),
        pos(all->set.getStorage()->getPosProperty()),
        pref1(pref),
        pairComputer(_pairComputer)

        {  } 

      virtual void operator()(const Reference pref2) {
   
        Real3D pos1(pos[pref1][0], pos[pref1][1], pos[pref1][2]);
        Real3D pos2(pos[pref2][0], pos[pref2][1], pos[pref2][2]);
   
        Real3D dist = bc.getDist(pos1, pos2);
   
        if (id[pref1] < id[pref2]) {
          pairComputer(dist, pref1, pref2);
        }
      }

     };

  public:

     PairComputer& pairComputer;

     const All* all;

     Traverser1(const All* _all, PairComputer& _pairComputer) :

       all(_all),
       pairComputer(_pairComputer)

      {  }

      virtual void operator()(const Reference pref) {

         // printf ("Traverser1: will call traverser2\n");

         Traverser2 traverser2(all, pref, pairComputer);

         all->set.foreach(traverser2);
      }
};

/*--------------------------------------------------------------------------
All::~All()
-------------------------------------------------------------------------- */

All::~All() {}

/*--------------------------------------------------------------------------
All::All(boundary_conditions, particle_set)
-------------------------------------------------------------------------- */

All::All(espresso::bc::BC& _bc, espresso::particleset::ParticleSet& _set):

  set(_set),
  bc(_bc)

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

