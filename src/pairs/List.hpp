#ifndef _PAIRS_LIST_HPP
#define _PAIRS_LIST_HPP

#include "Set.hpp"
#include "Computer.hpp"
#include "Particle.hpp"
#include "Property.hpp"
#include "bc/BC.hpp"

namespace espresso {
  namespace pairs {

    /** Class that applies a Computer to a list of pairs */

    class List : public Set {
    public:
      typedef shared_ptr< List > SelfPtr;

      static void registerPython();

      /* 2-parameter constructor */
      List(bc::BC::SelfPtr _bc, 
           storage::Storage::SelfPtr _storage1, 
           storage::Storage::SelfPtr _storage2, 
           Property< Real3D >::SelfPtr _posProperty1,
           Property< Real3D >::SelfPtr _posProperty2);

      /** Constructor
	  \param bc are the boundary conditions that are needed for distance calculation.
	  \param storage specifies the particle storage to which particles belong
	  \param posProperty the identifier of the posProperty property to use
      */
      List(bc::BC::SelfPtr _bc, 
           storage::Storage::SelfPtr _storage,
	   Property< Real3D >::SelfPtr _posProperty);

      /** Destructor. */
      virtual ~List();     

      /** return size of list */
      size_t size() const;

      /** Ask if a particle pair tuple (id1, id2) is in the pair list
	  \param id1 is the identificaiton of particle 1
	  \param id2 is the identificaiton of particle 2
	  \return true if tuple (id1, id2) is in the list
      */
      bool findPair(ParticleId id1, ParticleId id2) const;

      /** Adding a particle pair tuple (id1, id2) to the pair list
	  \param id1 is the identificaiton of particle 1
	  \param id2 is the identificaiton of particle 2
	  Note: a tuple (id1, id2) can be added several times.
      */
      void addPair(ParticleId id1, ParticleId id2);

      /** Deleting a particle pair tuple (id1, id2) from the pair list
	  \param id1 is the identificaiton of particle 1
	  \param id2 is the identificaiton of particle 2
	  Particle (id1, i2) must be in the pair list otherwise exception.
      */
      void deletePair(ParticleId id1, ParticleId id2);

      /** Getter routine for the storage */
      storage::Storage::SelfPtr getStorage() const { return storage; }

    protected:

      /** This routine will apply a function operator to all pairs of the list.
	  \param pairComputer is the object that provides the function to be applied.
      */
      virtual void foreachApply(Computer &computer);

    private:

      storage::Storage::SelfPtr storage;
      bc::BC::SelfPtr bc;
      Property< Real3D >::SelfPtr posProperty;

      typedef std::pair< ParticleId, ParticleId > Tuple;
      std::vector< Tuple > id_list;

    };

  }
}

#endif
