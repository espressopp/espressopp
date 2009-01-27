#ifndef _PAIRS_ALL_HPP
#define _PAIRS_ALL_HPP

#include "ParticlePairs.hpp"
#include "particlestorage/ParticleStorage.hpp"
#include "bc/BC.hpp"

#warning "pairs/All.hpp currently includes cstdio, remove as fast as possible"
#include <cstdio>

namespace espresso {

  namespace pairs {

     class All : public ParticlePairs {
 
       typedef espresso::particlestorage::ParticleStorage ParticleStorage;

       typedef espresso::particlestorage::ParticleStorage::PropertyReference<size_t> SizeRef;
       typedef espresso::particlestorage::ParticleStorage::ArrayPropertyReference<real> RealArrayRef;

     private:

       espresso::particleset::ParticleSet& set;
       espresso::bc::BC& bc;

       class ParticleSetTraverserComputer1 : public espresso::particlestorage::ParticleComputer {

          public:

          ParticlePairComputer& pairComputer;

          All* all;

          ParticleSetTraverserComputer1(All* _all, 
                                        ParticlePairComputer& _pairComputer) :

            all(_all),
            pairComputer(_pairComputer)

            {  }

          virtual void operator()(const ParticleStorage::reference pref) {

             printf ("Traverser1: will call traverser2\n");

             ParticleSetTraverserComputer2 traverser2(all, pref, pairComputer);

             all->set.foreach(traverser2);
          }

          virtual void operator()(const ParticleStorage::const_reference pref) const {

             printf ("const traversing1 particle\n");

          }
       };

       class ParticleSetTraverserComputer2 : public espresso::particlestorage::ParticleComputer {

          public:

          All* all;

          SizeRef id;
          RealArrayRef pos;

          const ParticleStorage::reference pref1;
          
          ParticlePairComputer& pairComputer;

          ParticleSetTraverserComputer2(All* _all,
                                        const ParticleStorage::reference pref,
                                         ParticlePairComputer& _pairComputer) :

            all(_all), 
            id(all->set.getStorage()->getIDProperty()),
            pos(all->set.getStorage()->getPosProperty()),
            pref1(pref),
            pairComputer(_pairComputer)

            {  } 

          virtual void operator()(const ParticleStorage::reference pref2) {

            Real3D pos1(pos[pref1][0], pos[pref1][1], pos[pref1][2]);
            Real3D pos2(pos[pref2][0], pos[pref2][1], pos[pref2][2]);

            Real3D dist = all->bc.getDist(pos1, pos2);

            /*
            printf ("traversing pair ids = (%d, %d), dist = (%f, %f, %f)\n",
                     id[pref1], id[pref2], 
                     dist.getX(), dist.getY(), dist.getZ());
            */

            if (id[pref1] < id[pref2]) {
              pairComputer(dist, pref1, pref2);
            }
          }

          virtual void operator()(const ParticleStorage::const_reference pref) const {

             printf ("const traversing2 particle\n");

          }
       };

     public:

       ~All() {}

       All (espresso::bc::BC& _bc, espresso::particleset::ParticleSet& _set):

          set(_set),
          bc(_bc)

       { }

       virtual void foreach(ParticlePairComputer& pairComputer) {

           ParticleSetTraverserComputer1 traverser1(this, pairComputer);;

           set.foreach(traverser1);
       }
       
       virtual void foreach(ConstParticlePairComputer& pairComputer) const {

           printf ("ParticlePairComputer const");
       }
     };
  }
}

#endif
