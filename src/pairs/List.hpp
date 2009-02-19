#ifndef _PAIRS_LIST_HPP
#define _PAIRS_LIST_HPP

#include "ParticlePairs.hpp"
#include "ParticlePairComputer.hpp"
#include "particlestorage/ParticleStorage.hpp"
#include "bc/BC.hpp"

namespace espresso {

  namespace pairs {

     /** Class that applies a ParticlePairComputer to a list of particle pairs 
     */

     class List : public ParticlePairs {
 
     private:

       espresso::particlestorage::ParticleStorage& storage; 

       espresso::bc::BC& bc;

       size_t coordinates;

       typedef std::pair<size_t,size_t> Tuple;

       std::vector<Tuple> id_list;

     public:

       /** Destructor. */

       ~List();

       /** Constructor for this class 

         \param bc are the boundary conditions that are needed for distance calculation.
         \param storage specifies the particle storage to which particles belong
	 \param coordinates the identifier of the coordinates property to use

       */

       List (espresso::bc::BC& bc, 
             espresso::particlestorage::ParticleStorage& storage, 
             size_t coordinates);

       int size();

       /** Ask if a particle pair tuple (id1, id2) is in the pair list

           \param id1 is the identificaiton of particle 1
           \param id2 is the identificaiton of particle 2
           \return true if tuple (id1, id2) is in the list
 
       */

       bool findPair(size_t id1, size_t id2);

       /** Adding a particle pair tuple (id1, id2) to the pair list

           \param id1 is the identificaiton of particle 1
           \param id2 is the identificaiton of particle 2
 
           Note: a tuple (id1, id2) can be added several times.
       */

       void addPair(size_t id1, size_t id2);

       /** Deleting a particle pair tuple (id1, id2) from the pair list

           \param id1 is the identificaiton of particle 1
           \param id2 is the identificaiton of particle 2
 
           Particle (id1, i2) must be in the pair list otherwise exception.
       */

       void deletePair(size_t id1, size_t id2);

       /** Getter routine for the boundary conditions. */

       espresso::bc::BC& getBC() const { return bc; }

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
