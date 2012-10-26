// ESPP_CLASS
#ifndef _INTERACTION_ANGULARUNIQUEPOTENTIAL_HPP
#define _INTERACTION_ANGULARUNIQUEPOTENTIAL_HPP

#include "types.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include <cmath>
#include "logging.hpp"

namespace espresso {
  namespace interaction {
    class AngularUniquePotential {
    public:
      virtual real computeEnergy(const Particle &p1, const Particle &p2, const Particle &p3, real cos0) const = 0;
      virtual real computeEnergy(const Real3D& dist12, const Real3D& dist32, real cos0) const = 0;
      virtual real computeEnergy(real theta, real cos0) const = 0;

      virtual void computeForce(Real3D& force12, Real3D& force32,
                                const Particle &p1, const Particle &p2, const Particle &p3, real cos0) const = 0;
      virtual void computeForce(Real3D& force12, Real3D& force32,
                                const Real3D& dist12, const Real3D& dist32, real cos0) const = 0;
      virtual real computeForce(real theta, real cos0) const = 0; // used for generating tabular file


      virtual void setCutoff(real _cutoff) = 0;
      virtual real getCutoff() const = 0;

      static void registerPython();

      static LOG4ESPP_DECL_LOGGER(theLogger);
    };

    template < class Derived >
    class AngularUniquePotentialTemplate : public AngularUniquePotential {
    public:
      AngularUniquePotentialTemplate();

      // Implements the Potential virtual interface
      virtual real computeEnergy(const Particle &p1, const Particle &p2, const Particle &p3, real cos0) const;
      virtual real computeEnergy(const Real3D& dist12, const Real3D& dist32, real cos0) const;
      virtual real computeEnergy(const real theta, real cos0) const;

      virtual void computeForce(Real3D& force12, Real3D& force32,
                                const Particle &p1, const Particle &p2, const Particle &p3, real cos0) const;
      virtual void computeForce(Real3D& force12, Real3D& force32,
                                const Real3D& dist12, const Real3D& dist32, real cos0) const;
      virtual real computeForce(real theta, real cos0) const; // used for generating tabular file

      virtual void setCutoff(real _cutoff);
      virtual real getCutoff() const;

      // Implements the non-virtual interface 
      // (used by e.g. the Interaction templates)
      real _computeEnergy(const Particle &p1, const Particle &p2, const Particle &p3, real cos0) const;
      real _computeEnergy(const Real3D& dist12, const Real3D& dist32, real cos0) const;
      real _computeEnergy(real theta, real cos0) const;

      void _computeForce(Real3D& force12, Real3D& force32,
			 const Particle &p1, const Particle &p2, const Particle &p3, real cos0) const;
      
      bool _computeForce(Real3D& force12, Real3D& force32,
			 const Real3D& dist12, const Real3D& dist32, real cos0) const;

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
    AngularUniquePotentialTemplate< Derived >::
    AngularUniquePotentialTemplate() 
      : cutoff(infinity), cutoffSqr(infinity)
    {}

    // Shift/cutoff handling
    template < class Derived >
    inline void
    AngularUniquePotentialTemplate< Derived >::
    setCutoff(real _cutoff) {
      cutoff = _cutoff;
      cutoffSqr = cutoff*cutoff;
    }

    template < class Derived >
    inline real
    AngularUniquePotentialTemplate< Derived >::
    getCutoff() const
    { return cutoff; }

    // Energy computation
    template < class Derived >
    inline real
    AngularUniquePotentialTemplate< Derived >::
    computeEnergy(const Particle &p1, const Particle &p2, const Particle &p3, real cos0) const {
      Real3D dist12 = p1.position() - p2.position();
      Real3D dist32 = p3.position() - p2.position();
      return computeEnergy(dist12, dist32, cos0);
    }

    template < class Derived >
    inline real
    AngularUniquePotentialTemplate< Derived >::
    computeEnergy(const Real3D& dist12, const Real3D& dist32, real cos0) const {
      real dist12Sqr = dist12 * dist12;
      real dist32Sqr = dist32 * dist32;
      real cos_theta = dist12 * dist32 / (sqrt(dist12Sqr) * sqrt(dist32Sqr));
      return computeEnergy(acos(cos_theta), cos0);
    }

    template < class Derived >
    inline real 
    AngularUniquePotentialTemplate< Derived >::
    computeEnergy(real theta, real cos0) const {
      return _computeEnergy(theta, cos0); // a bug was here (it was: return computeEnergy(theta);)
    }
    
    template < class Derived > 
    inline real 
    AngularUniquePotentialTemplate< Derived >::
    _computeEnergy(const Particle &p1, const Particle &p2, const Particle &p3, real cos0) const {
      Real3D dist12 = p1.position() - p2.position();
      Real3D dist32 = p3.position() - p2.position();
      return _computeEnergy(dist12, dist32, cos0);
    }

    template < class Derived >
    inline real
    AngularUniquePotentialTemplate< Derived >::
    _computeEnergy(const Real3D& dist12, const Real3D& dist32, real cos0) const {
      real dist12_sqr = dist12 * dist12;
      real dist32_sqr = dist32 * dist32;
      if (dist12_sqr >= cutoffSqr || dist32_sqr >= cutoffSqr )
        return 0.0;
      else{
        real cos_theta = dist12 * dist32 / (sqrt(dist12_sqr) * sqrt(dist32_sqr));
        return _computeEnergy(acos(cos_theta), cos0);
      }
    }

    template < class Derived > 
    inline real
    AngularUniquePotentialTemplate< Derived >::
    _computeEnergy(real theta, real cos0) const {
      return derived_this()->_computeEnergyRaw(theta, cos0);
    }
    
    // Force computation
    template < class Derived >
    inline void
    AngularUniquePotentialTemplate< Derived >::
    computeForce(Real3D& force12,
                 Real3D& force32,
                 const Particle &p1, const Particle &p2, const Particle &p3, real cos0) const {
      Real3D dist12 = p1.position() - p2.position();
      Real3D dist32 = p3.position() - p2.position();
      _computeForce(force12, force32, dist12, dist32, cos0);
    }

    template < class Derived >
    inline void
    AngularUniquePotentialTemplate< Derived >::
    computeForce(Real3D& force12,
                 Real3D& force32,
                 const Real3D& dist12,
                 const Real3D& dist32, real cos0) const {
      _computeForce(force12, force32, dist12, dist32, cos0);
    }

    template < class Derived >
    inline void
    AngularUniquePotentialTemplate< Derived >::
    _computeForce(Real3D& force12,
                  Real3D& force32,
                  const Particle &p1, const Particle &p2, const Particle &p3, real cos0) const {
      Real3D dist12 = p1.position() - p2.position();
      Real3D dist32 = p3.position() - p2.position();
      _computeForce(force12, force32, dist12, dist32, cos0);
    }

    template < class Derived >
    inline bool
    AngularUniquePotentialTemplate< Derived >::
    _computeForce(Real3D& force12,
                  Real3D& force32,
                  const Real3D& dist12,
                  const Real3D& dist32, real cos0) const {
      
      return derived_this()->_computeForceRaw(force12, force32, dist12, dist32, cos0);
    }
    
    // used for generating tabular angular potential
    template < class Derived >
    inline real
    AngularUniquePotentialTemplate< Derived >::
    computeForce(real theta, real cos0) const {
      return derived_this()->_computeForceRaw(theta, cos0);
    }
    
  }
}

#endif
