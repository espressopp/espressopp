#include "types.hpp"

#include "particles/Storage.hpp"
#include "pairs/Computer.hpp"

namespace espresso {
   namespace pairs {

     /** This class is used to print all particle pairs.

        \code
         PairWriteComputer pairWriteComputer(particlestorage);
         ParticlePairs pairs ...
         pairs->foreach(pairWriteComputer)
        \endcode

        \sa espresso::pairs::ParticlePairComputer
     */
     class PairWriteComputer: public ConstComputer {
       typedef particles::Storage Storage;
       typedef Storage::PropertyTraits<Storage::ParticleId>::ConstReference ConstParticleIdRef;
       typedef Storage::PropertyTraits<Real3D>::ConstReference ConstRealArrayRef;

     public:

       ConstRealArrayRef pos; //<! const property position
       ConstParticleIdRef id; //<! const property id(entification)

       /** Constructor of the pair writer class.
         
           \param particleStorage is needed to access the properties that will be printed

       */
       PairWriteComputer(const espresso::particles::Storage* particleStorage,
	                 espresso::particles::Storage::PropertyId position) :
         pos(particleStorage->getProperty<Real3D>(position)),
         id(particleStorage->getIDProperty())
       {}

       /** Implementation of the pure routine that is applied to each particle pair.

          \sa espresso::ParticlePairComputer::operator()
       */
       virtual void operator()(const Real3D &dist,
                               const espresso::particles::Set::const_reference p1,
                               const espresso::particles::Set::const_reference p2)

       {
          printf("Pair: id = (%ld,%ld) , pos1 = (%f,%f,%f), pos2 = (%f,%f,%f), dist = (%f, %f, %f)\n",

                 size_t(id[p1]), size_t(id[p2]),
                 pos[p1][0], pos[p1][1], pos[p1][2],
                 pos[p2][0], pos[p2][1], pos[p2][2],
                 dist[0], dist[1], dist[2]);
       }
     };
  }
}
