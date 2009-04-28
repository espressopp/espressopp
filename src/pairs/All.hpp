#ifndef _PAIRS_ALL_HPP
#define _PAIRS_ALL_HPP

#include "boost/shared_ptr.hpp"
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
       boost::shared_ptr<particles::Set> set;
       boost::shared_ptr<bc::BC> bc;
       boost::shared_ptr< Property<Real3D> > coordinates;

     public:
       static void registerPython();

       /** Destructor. */
       ~All();

       /** Constructor for this class 

         \param bc are the boundary conditions that are needed for distance calculation.
         \param set specifies the set of particles for which pairs will be considered.
	 \param coordinates the identifier of the coordinates property to use

       */
       All(boost::shared_ptr<bc::BC> _bc, 
           boost::shared_ptr<particles::Set> _set, 
           boost::shared_ptr< Property<Real3D> > _coordinates);

       /** Getter routine for the boundary conditions. */
       boost::shared_ptr<bc::BC> getBC() const { return bc; }

       /** Getter routine for the set of particles. */
       boost::shared_ptr<particles::Set> getSet() const { return set; }

       /** Getter routine for the coordinate property */
       boost::shared_ptr< const Property<Real3D> > getCoordinateProperty() const { return coordinates; }

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
