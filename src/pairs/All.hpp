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
       particles::Set::SelfPtr set;
       bc::BC::SelfPtr bc;
       Real3DProperty::SelfPtr coordinates;

     public:
       typedef shared_ptr< All > SelfPtr;

       static void registerPython();

       /** Destructor. */
       virtual ~All();

       /** Constructor for this class 

         \param bc are the boundary conditions that are needed for distance calculation.
         \param set specifies the set of particles for which pairs will be considered.
	 \param coordinates the identifier of the coordinates property to use

       */
       All(bc::BC::SelfPtr _bc, 
           particles::Set::SelfPtr _set, 
           Real3DProperty::SelfPtr _coordinates);

       /** Getter routine for the boundary conditions. */
       bc::BC::SelfPtr getBC() const { return bc; }

       /** Getter routine for the set of particles. */
       particles::Set::SelfPtr getSet() const { return set; }

       /** Getter routine for the coordinate property */
       Real3DProperty::SelfPtr getCoordinateProperty() const { return coordinates; }

       /** Getter routine for the storage */
       particles::Storage::SelfPtr getStorage() const { return set->getStorage(); }

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
