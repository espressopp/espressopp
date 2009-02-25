#include "types.hpp"

#include "particles/Storage.hpp"
#include "pairs/Computer.hpp"

namespace espresso {
   namespace pairs {

     typedef particles::Storage::PropertyTraits<size_t>::ConstReference ConstSizeRef;
     typedef particles::Storage::PropertyTraits<Real3D>::ConstReference ConstRealArrayRef;

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

       ConstRealArrayRef pos;  //<! const property position
       ConstSizeRef id;        //<! const property id(entification)

       /** Constructor of the pair writer class.
         
           \param particleStorage is needed to access the properties that will be printed

       */
       PairWriteComputer(const espresso::particles::Storage* particleStorage,
	                 size_t position) :
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

                 id[p1], id[p2],
                 pos[p1].getX(), pos[p1].getY(), pos[p1].getZ(),
                 pos[p2].getX(), pos[p2].getY(), pos[p2].getZ(),
                 dist.getX(), dist.getY(), dist.getZ());
       }
     };
  }
}
