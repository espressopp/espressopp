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
      bool _getEnergy(const Particle &p1, const Particle &p2,
		      const Real3D& dist,
		      const real distSqr, real& energy) const {
	return false;
      }

      /** Computes the force of the interaction for the two given
	  particles.  Returns whether the force is non-zero.
      */
      bool _getForce(const Particle &p1, const Particle &p2,
		     const Real3D& dist,
		     const real distSqr, Real3D& force) const {
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
      real getEnergy(const Particle &p1, const Particle &p2,
		     const Real3D& dist,
		     const real distSqr) const;
      void getForce(const Particle &p1, const Particle &p2,
		    const Real3D& dist,
		    const real distSqr, Real3D& force) const;

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
    getEnergy(const Particle &p1, const Particle &p2,
	      const Real3D& dist,
	      const real distSqr) const {
      derived_this()->getEnergy();
    }

  }
}

#endif
