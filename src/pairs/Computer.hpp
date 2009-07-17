#ifndef _PAIRS_COMPUTER_HPP
#define _PAIRS_COMPUTER_HPP

#include "types.hpp"
#include "particles/ParticleHandle.hpp"
#include "particles/Storage.hpp"
#include "boost/function.hpp"

namespace espresso {
  namespace pairs {
    /** Abstract class that defines the operator() applied to particle pairs
     */
    template< class ParticleHandle >
    class ComputerBase {
    public:
      /// @name extended function object interface
      //@{
      typedef Real3D first_argument_type;
      typedef ParticleHandle second_argument_type;
      typedef ParticleHandle  third_argument_type;
      typedef void                       result_type;
      //@}

      virtual void prepare(particles::Storage::SelfPtr storage1, 
			   particles::Storage::SelfPtr storage2) {}

      /** Interface of the routine that is applied to particle pairs
	  \param dist: distance vector between the two particles
	  \param p1, p2: references to the two particles

	  The pair computer gets the distance vector between the
	  particles, as it is hard to get otherwise and furthermore
	  some pair generators compute it anyway.

	  Note: The references are necessary if more property data of
	  the particles is needed than only the distance.
      */
      virtual void apply(const Real3D dist, 
			 const ParticleHandle p1, 
			 const ParticleHandle p2) = 0;

      virtual void finalize() {}
    };

    /** Abstract class that defines a function on pairs of particles */
    class Computer: 
      public ComputerBase< particles::ParticleHandle > 
    {
    public:
      typedef shared_ptr< Computer > SelfPtr;
      static void registerPython();
    };
    
    /** Abstract class that defines a function on pairs of read-only particles */
    class ConstComputer:
      public ComputerBase< particles::ConstParticleHandle > 
    {
    public:
      typedef shared_ptr< ConstComputer > SelfPtr;
    };
  }
}

#endif
