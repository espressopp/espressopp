#ifndef _PAIRS_ALL_HPP
#define _PAIRS_ALL_HPP

#include "Set.hpp"
#include "Computer.hpp"
#include "particles/Storage.hpp"
#include "bc/BC.hpp"

namespace espresso {
  namespace pairs {

     /** Class that applies a Computer to all particle pairs of
         a given particle set.
     */

     class All : public Set {
 
     private:
       particles::Set& set;
       bc::BC& bc;
       particles::PropertyId coordinates;

     public:
       static void registerPython();

       /** Destructor. */
       ~All();

       /** Constructor for this class 

         \param bc are the boundary conditions that are needed for distance calculation.
         \param set specifies the set of particles for which pairs will be considered.
	 \param coordinates the identifier of the coordinates property to use

       */
       All(bc::BC& _bc, particles::Set& _set, particles::PropertyId _coordinates);

       /** Getter routine for the boundary conditions. */
       bc::BC& getBC() const { return bc; }

       /** Getter routine for the set of particles. */
       particles::Set& getSet() const { return set; }

       /** Getter routine for the ID of the coordinate */
       particles::PropertyId getCoordinateProperty() const { return coordinates; }

       /** This routine will apply a function operator to all pairs.
         \param pairComputer is the object that provides the function to be applied to all pairs.
       */
       virtual void foreach(Computer& pairComputer);

       /** This routine will apply a function operator for read-only particles to all pairs.
         \param pairComputer is the object that provides the read-only function to be applied to all pairs.
       */
       virtual void foreach(ConstComputer& pairComputer) const;

     };
  }
}

#endif
