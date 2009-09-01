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
      List(bc::BC::SelfPtr bc,
	   storage::Storage::SelfPtr _storage1,
           storage::Storage::SelfPtr _storage2);

      /** Constructor
	  \param bc are the boundary conditions that are needed for distance calculation.
	  \param storage specifies the particle storage to which particles belong
	  \param posProperty the identifier of the posProperty property to use
      */
      List(bc::BC::SelfPtr bc,
	   storage::Storage::SelfPtr _storage);

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
      TODO: should addPair check to see if id1/2 are in range of storage?
      */
      void addPair(ParticleId id1, ParticleId id2);

      /** Deleting a particle pair tuple (id1, id2) from the pair list
	  \param id1 is the identificaiton of particle 1
	  \param id2 is the identificaiton of particle 2
	  Particle (id1, i2) must be in the pair list otherwise exception.
      */
      void deletePair(ParticleId id1, ParticleId id2);

      /** Getter routine for storage1 */
      virtual storage::Storage::SelfPtr getLeftStorage();
      /** Getter routine for storage2 */
      virtual storage::Storage::SelfPtr getRightStorage();

    protected:
      /** This routine will apply a function operator to all pairs of the list.
	  \param pairComputer is the object that provides the function to be applied.
      */
      virtual bool foreachPairApply(Computer &computer);

    private:
      bc::BC::SelfPtr bc;
      storage::Storage::SelfPtr storage1;
      storage::Storage::SelfPtr storage2;

      typedef std::pair< ParticleId, ParticleId > Tuple;
      std::vector< Tuple > id_list;
    };

  }
}

#endif
