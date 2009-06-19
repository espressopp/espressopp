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
 
     private:
       particles::PSet set;
       bc::PBC bc;
       PReal3DProperty coordinates;

     public:
       static void registerPython();

       /** Destructor. */
       virtual ~All();

       /** Constructor for this class 

         \param bc are the boundary conditions that are needed for distance calculation.
         \param set specifies the set of particles for which pairs will be considered.
	 \param coordinates the identifier of the coordinates property to use

       */
       All(bc::PBC _bc, 
           particles::PSet _set, 
           PReal3DProperty _coordinates);

       /** Getter routine for the boundary conditions. */
       bc::PBC getBC() const { return bc; }

       /** Getter routine for the set of particles. */
       particles::PSet getSet() const { return set; }

       /** Getter routine for the coordinate property */
       PReal3DProperty getCoordinateProperty() const { return coordinates; }

       /** Getter routine for the storage */
       particles::PStorage getStorage() const { return set->getStorage(); }

       /** This routine will apply a function operator to all pairs.
         \param pairComputer is the object that provides the function to be applied to all pairs.
       */
       virtual void foreach(Computer& pairComputer);

       /** This routine will apply a function operator for read-only particles to all pairs.
         \param pairComputer is the object that provides the read-only function to be applied to all pairs.
       */
       virtual void foreach(ConstComputer& pairComputer) const;

     };

    typedef boost::shared_ptr< All > PAll;
  }
}

#endif
