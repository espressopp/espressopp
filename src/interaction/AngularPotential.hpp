// ESPP_CLASS
#ifndef _INTERACTION_ANGULARPOTENTIAL_HPP
#define _INTERACTION_ANGULARPOTENTIAL_HPP

#include "types.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include <cmath>

namespace espresso {
  namespace interaction {
    class AngularPotential {
    public:
      virtual real computeEnergy(const Particle &p1, const Particle &p2, const Particle &p3) const = 0;
      virtual real computeEnergy(const Real3D& dist12, const Real3D& dist32) const = 0;
      virtual real computeEnergy(real theta) const = 0;

      virtual void computeForce(Real3D& force12, Real3D& force32,
                                const Particle &p1, const Particle &p2, const Particle &p3) const = 0;
      virtual void computeForce(Real3D& force12, Real3D& force32,
                                const Real3D& dist12, const Real3D& dist32) const = 0;
      virtual real computeForce(real theta) const = 0; // used for generating tabular file


      virtual void setCutoff(real _cutoff) = 0;
      virtual real getCutoff() const = 0;

      static void registerPython();
    };

    // enum PotentialType {
    //   Default,
    //   NoScalarDistance,
    //   NoDistance
    // };

    //    template < class Derived, enum PotentialType = Default >
    /** Provides a template for the simple implementation of a
    shifted, absolute distance dependent potential with cutoff.
    */
    template < class Derived >
    class AngularPotentialTemplate : public AngularPotential {
    public:
      AngularPotentialTemplate();

      // Implements the Potential virtual interface
      virtual real computeEnergy(const Particle &p1, const Particle &p2, const Particle &p3) const;
      virtual real computeEnergy(const Real3D& dist12, const Real3D& dist32) const;
      virtual real computeEnergy(const real theta) const;

      virtual void computeForce(Real3D& force12, Real3D& force32,
                                const Particle &p1, const Particle &p2, const Particle &p3) const;
      virtual void computeForce(Real3D& force12, Real3D& force32,
                                const Real3D& dist12, const Real3D& dist32) const;
      virtual real computeForce(real theta) const; // used for generating tabular file

      virtual void setCutoff(real _cutoff);
      virtual real getCutoff() const;

      // Implements the non-virtual interface 
      // (used by e.g. the Interaction templates)
      real _computeEnergy(const Particle &p1, const Particle &p2, const Particle &p3) const;
      real _computeEnergy(const Real3D& dist12, const Real3D& dist32) const;
      real _computeEnergy(real theta) const;

      void _computeForce(Real3D& force12, Real3D& force32,
			 const Particle &p1, const Particle &p2, const Particle &p3) const;
      void _computeForce(Real3D& force12, Real3D& force32,
			 const Real3D& dist12, const Real3D& dist32) const;

      // Requires the following non-virtual interface in Derived
      // real _computeEnergySqrRaw(real distSqr) const;
      // bool _computeForceRaw(const Real3D& dist, Real3D& force) const;

      // void _computeForce(const Particle &p1, const Particle &p2, 
      //                    Real3D& force) const {
      // 	Real3D dist = p1.r.p - p2.r.p;
      // 	derived_this()->_computeForce(dist, force);
      // }

    protected:
      real cutoff;
      real cutoffSqr;

      Derived* derived_this() {
        return static_cast< Derived* >(this);
      }

    const Derived* derived_this() const {
        return static_cast< const Derived* >(this);
      }
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < class Derived > 
    inline
    AngularPotentialTemplate< Derived >::
    AngularPotentialTemplate() 
      : cutoff(infinity), cutoffSqr(infinity)
    {}

    // Shift/cutoff handling
    template < class Derived >
    inline void
    AngularPotentialTemplate< Derived >::
    setCutoff(real _cutoff) {
      cutoff = _cutoff;
      cutoffSqr = cutoff*cutoff;
    }

    template < class Derived >
    inline real
    AngularPotentialTemplate< Derived >::
    getCutoff() const
    { return cutoff; }

    // Energy computation
    template < class Derived >
    inline real
    AngularPotentialTemplate< Derived >::
    computeEnergy(const Particle &p1, const Particle &p2, const Particle &p3) const {
      Real3D dist12 = p1.position() - p2.position();
      Real3D dist32 = p3.position() - p2.position();
      return computeEnergy(dist12, dist32);
    }

    template < class Derived >
    inline real
    AngularPotentialTemplate< Derived >::
    computeEnergy(const Real3D& dist12, const Real3D& dist32) const {
      real dist12Sqr = dist12 * dist12;
      real dist32Sqr = dist32 * dist32;
      real cos_theta = dist12 * dist32 / (sqrt(dist12Sqr) * sqrt(dist32Sqr));
      return computeEnergy(acos(cos_theta));
    }

    template < class Derived >
    inline real 
    AngularPotentialTemplate< Derived >::
    computeEnergy(real theta) const {
      return _computeEnergy(theta); // a bug was here (it was: return computeEnergy(theta);)
    }
    
    template < class Derived > 
    inline real 
    AngularPotentialTemplate< Derived >::
    _computeEnergy(const Particle &p1, const Particle &p2, const Particle &p3) const {
      Real3D dist12 = p1.position() - p2.position();
      Real3D dist32 = p3.position() - p2.position();
      return _computeEnergy(dist12, dist32);
    }

    template < class Derived >
    inline real
    AngularPotentialTemplate< Derived >::
    _computeEnergy(const Real3D& dist12, const Real3D& dist32) const {

      real dist12_sqr = dist12 * dist12;
      real dist32_sqr = dist32 * dist32;
      real cos_theta = dist12 * dist32 / (sqrt(dist12_sqr) * sqrt(dist32_sqr));
      return _computeEnergy(acos(cos_theta));
    }

    template < class Derived > 
    inline real
    AngularPotentialTemplate< Derived >::
    _computeEnergy(real theta) const {
      return derived_this()->_computeEnergyRaw(theta);
    }
    
    // Force computation
    template < class Derived >
    inline void
    AngularPotentialTemplate< Derived >::
    computeForce(Real3D& force12,
                 Real3D& force32,
                 const Particle &p1, const Particle &p2, const Particle &p3) const {
      Real3D dist12 = p1.position() - p2.position();
      Real3D dist32 = p3.position() - p2.position();
      _computeForce(force12, force32, dist12, dist32);
    }

    template < class Derived >
    inline void
    AngularPotentialTemplate< Derived >::
    computeForce(Real3D& force12,
                 Real3D& force32,
                 const Real3D& dist12,
                 const Real3D& dist32) const {
      _computeForce(force12, force32, dist12, dist32);
    }

    template < class Derived >
    inline void
    AngularPotentialTemplate< Derived >::
    _computeForce(Real3D& force12,
                  Real3D& force32,
                  const Particle &p1, const Particle &p2, const Particle &p3) const {
      Real3D dist12 = p1.position() - p2.position();
      Real3D dist32 = p3.position() - p2.position();
      _computeForce(force12, force32, dist12, dist32);
    }

    template < class Derived >
    inline void
    AngularPotentialTemplate< Derived >::
    _computeForce(Real3D& force12,
                  Real3D& force32,
                  const Real3D& dist12,
                  const Real3D& dist32) const {
      derived_this()->_computeForceRaw(force12, force32, dist12, dist32);
    }
    
    // used for generating tabular angular potential
    template < class Derived >
    inline real
    AngularPotentialTemplate< Derived >::
    computeForce(real theta) const {
      return derived_this()->_computeForceRaw(theta);
    }
    
  }
}

#endif
