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
 
     private:

       espresso::particleset::ParticleSet& set; 

       espresso::bc::BC& bc;

       size_t coordinates;

     public:

       /** Destructor. */

       ~All();

       /** Constructor for this class 

         \param bc are the boundary conditions that are needed for distance calculation.
         \param set specifies the set of particles for which pairs will be considered.
	 \param coordinates the identifier of the coordinates property to use

       */

       All (espresso::bc::BC& bc, espresso::particleset::ParticleSet& set, size_t coordinates);

       /** Getter routine for the boundary conditions. */

       espresso::bc::BC& getBC() const { return bc; }

       /** Getter routine for the set of particles. */

       espresso::particleset::ParticleSet& getSet() const { return set; }

       /** Getter routine for the ID of the coordinate */

       size_t getCoordinateProperty() const {return coordinates; }

       /** This routine will apply a function operator to all pairs.

         \param pairComputer is the object that provides the function to be applied to all pairs.

       */

       virtual void foreach(ParticlePairComputer& pairComputer);

       /** This routine will apply a function operator for read-only particles to all pairs.

         \param pairComputer is the object that provides the read-only function to be applied to all pairs.

       */

       virtual void foreach(ConstParticlePairComputer& pairComputer) const;

     };
  }
}

#endif
