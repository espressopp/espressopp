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
     public:
       particles::ConstPropertyReference<Real3D>                pos; //<! const property position
       particles::ConstPropertyReference<particles::ParticleId> id; //<! const property id(entification)

       /** Constructor of the pair writer class.
         
           \param particleStorage is needed to access the properties that will be printed

       */
       PairWriteComputer(const particles::Storage* particleStorage,
	                 particles::PropertyId position) :
         pos(particleStorage->getPropertyReference<Real3D>(position)),
         id(particleStorage->getIDProperty())
       {}

       /** Implementation of the pure routine that is applied to each particle pair.

          \sa espresso::ParticlePairComputer::operator()
       */
       virtual void operator()(const Real3D &dist,
                               const particles::ConstParticleReference p1,
                               const particles::ConstParticleReference p2)

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
