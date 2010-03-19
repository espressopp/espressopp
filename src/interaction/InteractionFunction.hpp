#ifndef _INTERACTION_INTERACTIONFUNCTION_HPP
#define _INTERACTION_INTERACTIONFUNCTION_HPP

#include "types.hpp"

namespace espresso {
  namespace interaction {

    /**
       An interaction function needs to provide the following
       interface, which is called from the different
       InteractionTemplates. 
    */
    class InteractionFunctionBase {
      /** Computes the energy of the interaction for the two given
	  particles. 
	  Returns whether the energy is non-zero.
      */
      bool _getEnergy(Particle &p1, Particle &p2,
		      const Real3DRef dist,
		      const real distSqr, real& energy) const {
	return false;
      }

      /** Computes the force of the interaction for the two given
	  particles.  Returns whether the force is non-zero.
      */
      bool _getForce(Particle &p1, Particle &p2,
		     const Real3DRef dist,
		     const real distSqr, Real3DRef force) const {
	return false;
      }
    };

    /** Interaction function base class. */
    class DistanceDependentInteractionFunction {
    private:
      real cutoff;
      real cutoffSqr;
      
    public:
      InteractionFunction() { setCutoff(0.0); }

      void setCutoff(real _cutoff) { cutoff = _cutoff; cutoffSqr = cutoff * cutoff; }
      real getCutoff() const { return cutoff; }
      real getCutoffSqr() const { return cutoffSqr; }
    };

    template < class Derived > 
    class InteractionFunctionTemplate {
    public:
      real getEnergy(Particle &p1, Particle &p2,
		     const Real3DRef dist,
		     const real distSqr) const;
      void getForce(Particle &p1, Particle &p2,
		    const Real3DRef dist,
		    const real distSqr, Real3DRef force) const;

    protected:
      Derived* derived_this() {
	return static_cast< Derived* >(this);
      }
      
      const Derived* derived_this() const {
	return static_cast< const Derived* >(this);
      }
    };

    template < class Derived > 
    inline real
    InteractionFunctionTemplate< Derived >::
    getEnergy(Particle &p1, Particle &p2,
	      const Real3DRef dist,
	      const real distSqr) const {
      derived_this()->getEnergy();
    }

  }
}

#endif
