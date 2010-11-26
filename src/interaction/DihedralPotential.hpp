// ESPP_CLASS
#ifndef _INTERACTION_DIHEDRALPOTENTIAL_HPP
#define _INTERACTION_DIHEDRALPOTENTIAL_HPP

#include "types.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include <cmath>

namespace espresso {
  namespace interaction {
    class DihedralPotential {
    public:
      virtual real computeEnergy(const Real3D& dist21,
                                 const Real3D& dist32,
                                 const Real3D& dist43) const = 0;

      virtual real computeEnergy(real phi) const = 0;

      virtual void computeForce(Real3D& force1,
                                Real3D& force2,
                                Real3D& force3,
                                Real3D& force4,
                                const Real3D& dist21,
                                const Real3D& dist32,
                                const Real3D& dist43) const = 0;

      virtual void setCutoff(real _cutoff) = 0;
      virtual real getCutoff() const = 0;

      static void registerPython();
    };

    template < class Derived >
    class DihedralPotentialTemplate : public DihedralPotential {
    public:
      DihedralPotentialTemplate();

      // Implements the Potential virtual interface
      virtual real computeEnergy(const Real3D& dist21,
                                 const Real3D& dist32,
                                 const Real3D& dist43) const;

      virtual real computeEnergy(real phi) const;

      virtual void computeForce(Real3D& force1,
                                Real3D& force2,
                                Real3D& force3,
                                Real3D& force4,
                                const Real3D& dist21,
                                const Real3D& dist32,
                                const Real3D& dist43) const;

      virtual void setCutoff(real _cutoff);
      virtual real getCutoff() const;

      // Implements the non-virtual interface 
      // (used by e.g. the Interaction templates)
      real _computeEnergy(const Real3D& dist21,
                          const Real3D& dist32,
                          const Real3D& dist43) const;

      real _computeEnergy(real phi) const;

      void _computeForce(Real3D& force1,
                         Real3D& force2,
                         Real3D& force3,
                         Real3D& force4,
			 const Real3D& dist21,
			 const Real3D& dist32,
                         const Real3D& dist43) const;

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
    DihedralPotentialTemplate< Derived >::
    DihedralPotentialTemplate() 
      : cutoff(infinity), cutoffSqr(infinity)
    {}

    // Shift/cutoff handling
    template < class Derived >
    inline void
    DihedralPotentialTemplate< Derived >::
    setCutoff(real _cutoff) {
      cutoff = _cutoff;
      cutoffSqr = cutoff*cutoff;
    }

    template < class Derived >
    inline real
    DihedralPotentialTemplate< Derived >::
    getCutoff() const
    { return cutoff; }

    // Energy computation
    template < class Derived >
    inline real
    DihedralPotentialTemplate< Derived >::
    computeEnergy(const Real3D& dist21,
                  const Real3D& dist32,
                  const Real3D& dist43) const {
        // compute phi
      return computeEnergy(-1.0);
    }

    template < class Derived >
    inline real
    DihedralPotentialTemplate< Derived >::
    computeEnergy(real phi) const {
      return _computeEnergy(phi);
    }
    
    template < class Derived >
    inline real
    DihedralPotentialTemplate< Derived >::
    _computeEnergy(const Real3D& dist21,
                   const Real3D& dist32,
                   const Real3D& dist43) const {
        // compute phi
      return _computeEnergy(-1.0);
    }

    template < class Derived > 
    inline real
    DihedralPotentialTemplate< Derived >::
    _computeEnergy(real phi) const {
      return derived_this()->_computeEnergyRaw(phi);
    }
    
    // Force computation
    template < class Derived >
    inline void
    DihedralPotentialTemplate< Derived >::
    computeForce(Real3D& force1,
                 Real3D& force2,
                 Real3D& force3,
                 Real3D& force4,
                 const Real3D& dist21,
                 const Real3D& dist32,
                 const Real3D& dist43) const {
      _computeForce(force1, force2, force3, force4, dist21, dist32, dist43);
    }

    template < class Derived >
    inline void
    DihedralPotentialTemplate< Derived >::
    _computeForce(Real3D& force1,
                  Real3D& force2,
                  Real3D& force3,
                  Real3D& force4,
                  const Real3D& dist21,
                  const Real3D& dist32,
                  const Real3D& dist43) const {
      derived_this()->_computeForceRaw(force1, force2, force3, force4, dist21, dist32, dist43);
    }
  }
}

#endif
