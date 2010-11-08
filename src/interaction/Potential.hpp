// ESPP_CLASS
#ifndef _INTERACTION_POTENTIAL_HPP
#define _INTERACTION_POTENTIAL_HPP

#include "types.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"

namespace espresso {
  namespace interaction {
    class Potential {
    public:
      virtual real computeEnergy(Particle &p1, Particle &p2) const = 0;
      virtual real computeEnergy(ConstReal3DRef dist) const = 0;
      virtual real computeEnergy(real dist) const = 0;
      virtual real computeEnergySqr(real distSqr) const = 0;
      
      virtual Real3D computeForce(Particle &p1, Particle &p2) const = 0;
      virtual Real3D computeForce(ConstReal3DRef dist) const = 0;

      virtual void setCutoff(real _cutoff) = 0;
      virtual real getCutoff() const = 0;

      virtual void setShift(real _shift) = 0;
      virtual real getShift() const = 0;
      virtual real setAutoShift() = 0;

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
    class PotentialTemplate : public Potential {
    public:
      PotentialTemplate();

      // Implements the Potential virtual interface
      virtual real computeEnergy(Particle &p1, Particle &p2) const;
      virtual real computeEnergy(ConstReal3DRef dist) const;
      virtual real computeEnergy(real dist) const;
      virtual real computeEnergySqr(real distsqr) const;
      virtual Real3D computeForce(Particle &p1, Particle &p2) const;
      virtual Real3D computeForce(ConstReal3DRef dist) const;
      virtual void setCutoff(real _cutoff);
      virtual real getCutoff() const;
      virtual void setShift(real _shift);
      virtual real getShift() const;
      virtual real setAutoShift();
      void updateAutoShift();

      // Implements the non-virtual interface 
      // (used by e.g. the Interaction templates)
      real _computeEnergy(Particle &p1, Particle &p2) const;
      real _computeEnergy(ConstReal3DRef dist) const;
      real _computeEnergy(real dist) const;
      real _computeEnergySqr(real distSqr) const;

      bool _computeForce(real force[3], 
			 Particle &p1, Particle &p2) const;
      bool _computeForce(Real3DRef force, 
			 ConstReal3DRef dist) const;

      // Requires the following non-virtual interface in Derived
      // real _computeEnergySqrRaw(real distSqr) const;
      // bool _computeForceRaw(ConstReal3DRef dist, Real3DRef force) const;


      // void _computeForce(Particle &p1, Particle &p2, 
      // 			 Real3DRef force) const {
      // 	Real3D dist = p1.r.p - p2.r.p;
      // 	derived_this()->_computeForce(dist, force);
      // }

    protected:
      real cutoff;
      real cutoffSqr;
      real shift;
      bool autoShift;

      Derived* derived_this() {
	return static_cast< Derived* >(this);
      }

      const Derived* derived_this() const {
	return static_cast< const Derived* >(this);
      }
    };

    // template < class Derived >
    // class PotentialTemplate< Derived, NoScalarDistance > 
    //   : public PotentialTemplate< Derived, Default > {
    // public:
    //   typedef PotentialTemplate< Derived, NoScalarDistance > Super;
    //   using PotentialTemplate< Derived >::_computeEnergy;
    //   using PotentialTemplate< Derived >::computeEnergy;

    //   real _computeEnergySqr(real distsqr) const {
    // 	throw Potential::BadCall();
    //   }

    //   real computeEnergy(ConstReal3DRef dist) {
    // 	return _computeEnergy(dist);
    //   }
    // };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < class Derived > 
    inline
    PotentialTemplate< Derived >::
    PotentialTemplate() 
      : cutoff(infinity), cutoffSqr(infinity),
	shift(0.0), autoShift(false)
    {}

    // Shift/cutoff handling
    template < class Derived > 
    inline void 
    PotentialTemplate< Derived >::
    setCutoff(real _cutoff) {
      cutoff = _cutoff;
      cutoffSqr = cutoff*cutoff;
      updateAutoShift();
    }

    template < class Derived > 
    inline real 
    PotentialTemplate< Derived >::
    getCutoff() const 
    { return cutoff; }

    template < class Derived > 
    inline void 
    PotentialTemplate< Derived >::
    setShift(real _shift) { 
      autoShift = false; 
      shift = _shift; 
    }

    template < class Derived > 
    inline real 
    PotentialTemplate< Derived >::
    getShift() const 
    { return shift; }

    template < class Derived > 
    inline real 
    PotentialTemplate< Derived >::
    setAutoShift() {
      autoShift = true;
      if (cutoffSqr == infinity) 
	shift = 0.0;
      else 
	shift = derived_this()->_computeEnergySqrRaw(cutoffSqr);
      return shift;
    }

    template < class Derived > 
    inline void 
    PotentialTemplate< Derived >::
    updateAutoShift() {
      if (autoShift) setAutoShift();
    }

    // Energy computation
    template < class Derived > 
    inline real 
    PotentialTemplate< Derived >::
    computeEnergy(Particle &p1, Particle &p2) const {
      Real3D dist = Real3DRef(p1.r.p) - Real3DRef(p2.r.p);
      return computeEnergy(dist);
    }

    template < class Derived > 
    inline real 
    PotentialTemplate< Derived >::
    computeEnergy(ConstReal3DRef dist) const {
	return computeEnergySqr(dist.sqr());
      }

    template < class Derived > 
    inline real 
    PotentialTemplate< Derived >::
    computeEnergy(real dist) const {
      return computeEnergySqr(dist*dist);
    }
    
    template < class Derived > 
    inline real 
    PotentialTemplate< Derived >::
    computeEnergySqr(real distsqr) const {
      return _computeEnergySqr(distsqr);
    }

    template < class Derived > 
    inline real 
    PotentialTemplate< Derived >::
    _computeEnergy(Particle &p1, Particle &p2) const {
      Real3D dist = Real3DRef(p1.r.p) - Real3DRef(p2.r.p);
      return _computeEnergy(dist);
    }

    template < class Derived > 
    inline real 
    PotentialTemplate< Derived >::
    _computeEnergy(ConstReal3DRef dist) const {
      return _computeEnergySqr(dist.sqr());
    }

    template < class Derived > 
    inline real 
    PotentialTemplate< Derived >::
    _computeEnergy(real dist) const {
      return _computeEnergySqr(dist*dist);
    }

    template < class Derived > 
    inline real
    PotentialTemplate< Derived >::
    _computeEnergySqr(real distSqr) const {
      if (distSqr > cutoffSqr) 
	return 0.0;
      else
	return derived_this()->_computeEnergySqrRaw(distSqr) - shift;
    }
    
    // Force computation
    template < class Derived > 
    inline Real3D 
    PotentialTemplate< Derived >::
    computeForce(Particle &p1, Particle &p2) const {
      Real3D dist = Real3DRef(p1.r.p) - Real3DRef(p2.r.p);
      return computeForce(dist);
    }
    
    template < class Derived > 
    inline Real3D 
    PotentialTemplate< Derived >::
    computeForce(ConstReal3DRef dist) const {
      Real3D force;
      if (!_computeForce(force, dist))
	force = 0.0;
      return force;
    }

    template < class Derived > 
    inline bool
    PotentialTemplate< Derived >::
    _computeForce(real force[3], Particle &p1, Particle &p2) const {
      real dist[3];
      dist[0] = p1.r.p[0] - p2.r.p[0];
      dist[1] = p1.r.p[1] - p2.r.p[1];
      dist[2] = p1.r.p[2] - p2.r.p[2];
      real distSqr = dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2];
      if (distSqr > cutoffSqr) return false;
      derived_this()->_computeForceRaw(force, dist, distSqr);
      return true;
    }
    
    template < class Derived > 
    inline bool 
    PotentialTemplate< Derived >::
    _computeForce(Real3DRef force, ConstReal3DRef dist) const {
      real distSqr = dist.sqr();
      if (distSqr > cutoffSqr) return false;
      derived_this()->_computeForceRaw(force.get(), dist.get(), distSqr);
      return true;
    }
  }
}

#endif
