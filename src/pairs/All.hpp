#ifndef _PAIRS_ALL_HPP
#define _PAIRS_ALL_HPP

#include "Set.hpp"
#include "Computer.hpp"
#include "Property.hpp"
#include "bc/BC.hpp"

namespace espresso {
  namespace pairs {
    /** Class that applies a Computer to all particle pairs of
	a given particle set.
    */
    class All : public Set {
    public:
      typedef shared_ptr< All > SelfPtr;

      static void registerPython();

      /** Constructor for this class 

	  \param set specifies the set of particles for which pairs will be considered.

      */
      All(particles::Set::SelfPtr _set);

      virtual ~All() {}

      /** Getter routine for storage1 */
      virtual storage::Storage::SelfPtr getLeftStorage();
      /** Getter routine for storage2 */
      virtual storage::Storage::SelfPtr getRightStorage();

    protected:
      /** This routine will apply a function operator to all pairs.
	  \param pairComputer is the object that provides the function to be applied to all pairs.
      */
      virtual bool foreachPairApply(Computer &computer);

    private:
      particles::Set::SelfPtr set;
    };

  }
}

#endif
