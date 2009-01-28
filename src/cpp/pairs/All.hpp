#ifndef _PAIRS_ALL_HPP
#define _PAIRS_ALL_HPP

#include "ParticlePairs.hpp"
#include "ParticlePairComputer.hpp"
#include "particlestorage/ParticleStorage.hpp"
#include "bc/BC.hpp"

namespace espresso {

  namespace pairs {

     /** Class that applies a ParticlePairComputer to all particle pairs of
         a given particle set.
     */

     class All : public ParticlePairs {
 
     public:

       espresso::particleset::ParticleSet& set; 
       espresso::bc::BC& bc;

       /** Destructor. */

       ~All();

       /** Constructor for this class 

         \param bc are the boundary conditions that are needed for distance calculation.
         \param set specifies the set of particles for which pairs will be considered.

       */

       All (espresso::bc::BC& bc, espresso::particleset::ParticleSet& set);

       /** This routine will apply a function operator to all pairs.

         \param pairComputer is the object that provides the function to be applied to all pairs.

       */

       virtual void foreach(ParticlePairComputer& pairComputer);

       virtual void foreach(ConstParticlePairComputer& pairComputer) const;

     };
  }
}

#endif
