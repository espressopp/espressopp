/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

// ESPP_CLASS
#ifndef _INTERACTION_PotentialUniqueDistUNIQUEDIST_HPP
#define _INTERACTION_PotentialUniqueDistUNIQUEDIST_HPP

#include "types.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "logging.hpp"

namespace espressopp {
  namespace interaction {
    class PotentialUniqueDist {
    public:
      virtual ~PotentialUniqueDist() {};
      virtual real computeEnergy(const Particle &p1, const Particle &p2, const real curDist) const = 0;
      virtual real computeEnergy(const Real3D& r21, const real curDist) const = 0;
      virtual real computeEnergy(real dist, real curDist) const = 0;
      virtual real computeEnergySqr(real distSqr, real curDist) const = 0;
      
      virtual Real3D computeForce(const Particle &p1, const Particle &p2, const real curDist) const = 0;
      virtual Real3D computeForce(const Real3D& dist, const real curDist) const = 0;

      virtual void setCutoff(real _cutoff) = 0;
      virtual real getCutoff() const = 0;

      virtual void setShift(real _shift) = 0;
      virtual real getShift() const = 0;
      virtual real setAutoShift() = 0;

      static void registerPython();

      static LOG4ESPP_DECL_LOGGER(theLogger);
    };

    // enum PotentialUniqueDistType {
    //   Default,
    //   NoScalarDistance,
    //   NoDistance
    // };

    //    template < class Derived, enum PotentialUniqueDistType = Default >
    /** Provides a template for the simple implementation of a
    shifted, absolute distance dependent PotentialUniqueDist with cutoff.
    */
    template < class Derived >
    class PotentialUniqueDistTemplate : public PotentialUniqueDist {
    public:
      PotentialUniqueDistTemplate();
      virtual ~PotentialUniqueDistTemplate() {};

      // Implements the PotentialUniqueDist virtual interface
      virtual real computeEnergy(const Particle &p1, const Particle &p2, const real curDist) const;
      virtual real computeEnergy(const Real3D& r21, const real curDist) const;
      virtual real computeEnergy(real dist, real curDist) const;
      virtual real computeEnergySqr(real distsqr, real curDist) const;
      virtual Real3D computeForce(const Particle &p1, const Particle &p2, const real curDist) const;
      virtual Real3D computeForce(const Real3D& r21, const real curDist) const;
      virtual void setCutoff(real _cutoff);
      virtual real getCutoff() const;
      virtual void setShift(real _shift);
      virtual real getShift() const;
      virtual real setAutoShift();
      void updateAutoShift();

      // Implements the non-virtual interface 
      // (used by e.g. the Interaction templates)
      real _computeEnergy(const Particle &p1, const Particle &p2, const real curDist) const;
      real _computeEnergy(const Real3D& r21, const real curDist) const;
      real _computeEnergy(real dist, real curDist) const;
      real _computeEnergySqr(real distSqr, real curDist) const;

      bool _computeForce(Real3D& force, 
			 const Particle &p1, const Particle &p2, const real curDist) const;
      bool _computeForce(Real3D& force, 
			 const Real3D& dist, const real curDist) const;
      
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

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < class Derived > 
    inline
    PotentialUniqueDistTemplate< Derived >::PotentialUniqueDistTemplate() : cutoff(infinity), cutoffSqr(infinity), shift(0.0), autoShift(false){
    }

    // Shift/cutoff handling
    template < class Derived > 
    inline void 
    PotentialUniqueDistTemplate< Derived >::setCutoff(real _cutoff) {
      cutoff = _cutoff;
      cutoffSqr = cutoff*cutoff;
      updateAutoShift();
    }

    template < class Derived > 
    inline real 
    PotentialUniqueDistTemplate< Derived >::getCutoff() const {
        return cutoff;
    }

    template < class Derived > 
    inline void 
    PotentialUniqueDistTemplate< Derived >::setShift(real _shift) { 
      autoShift = false; 
      shift = _shift; 
    }

    template < class Derived > 
    inline real 
    PotentialUniqueDistTemplate< Derived >::
    getShift() const 
    { return shift; }

    template < class Derived > 
    inline real 
    PotentialUniqueDistTemplate< Derived >::
    setAutoShift() {
      std::cout<<"Warning! Auto Shift does not work for this kind of potential"<<std::endl;
      
      /*
      autoShift = true;
      if (cutoffSqr == infinity) 
	    shift = 0.0;
      else 
	    shift = derived_this()->_computeEnergySqrRaw(cutoffSqr);
       */
      return 0.0;
    }

    template < class Derived > 
    inline void 
    PotentialUniqueDistTemplate< Derived >::
    updateAutoShift() {
      if (autoShift) setAutoShift();
    }

    // Energy computation
    template < class Derived > 
    inline real 
    PotentialUniqueDistTemplate< Derived >::
    computeEnergy(const Particle& p1, const Particle& p2, const real curDist) const {
      Real3D dist = p1.position() - p2.position();
      return computeEnergy(dist, curDist);
    }

    template < class Derived > 
    inline real 
    PotentialUniqueDistTemplate< Derived >::
    computeEnergy(const Real3D& r21, const real curDist) const {
	return computeEnergySqr(r21.sqr(), curDist);
      }

    template < class Derived > 
    inline real 
    PotentialUniqueDistTemplate< Derived >::
    computeEnergy(real dist, real curDist) const {
      return computeEnergySqr(dist*dist, curDist);
    }
    
    template < class Derived > 
    inline real 
    PotentialUniqueDistTemplate< Derived >::
    computeEnergySqr(real distsqr, real curDist) const {
      return _computeEnergySqr(distsqr, curDist);
    }

    template < class Derived > 
    inline real 
    PotentialUniqueDistTemplate< Derived >::
    _computeEnergy(const Particle& p1, const Particle& p2, const real curDist) const {
      Real3D dist = p1.position() - p2.position();
      return _computeEnergy(dist, curDist);
    }

    template < class Derived > 
    inline real 
    PotentialUniqueDistTemplate< Derived >::
    _computeEnergy(const Real3D& dist, const real curDist) const {
      return _computeEnergySqr(dist.sqr(), curDist);
    }

    template < class Derived > 
    inline real 
    PotentialUniqueDistTemplate< Derived >::
    _computeEnergy(real dist, real curDist) const {
      return _computeEnergySqr(dist*dist, curDist);
    }

    template < class Derived > 
    inline real
    PotentialUniqueDistTemplate< Derived >::
    _computeEnergySqr(real distSqr, real curDist) const {
      if (distSqr > cutoffSqr) 
        return 0.0;
      else {
        real e = derived_this()->_computeEnergySqrRaw(distSqr, curDist) - shift;
        LOG4ESPP_TRACE(theLogger, "Epot(r*r=" << distSqr << ") = " << e);
        return e;
      }
    }
    
    // Force computation
    template < class Derived > 
    inline Real3D 
    PotentialUniqueDistTemplate< Derived >::
    computeForce(const Particle& p1, const Particle& p2, const real curDist) const {
      Real3D dist = p1.position() - p2.position();
      return computeForce(dist, curDist);
    }
    
    template < class Derived > 
    inline Real3D 
    PotentialUniqueDistTemplate< Derived >::
    computeForce(const Real3D& r21, const real curDist) const {
      Real3D force;
      if(!_computeForce(force, r21, curDist))
        force = 0.0;
      return force;
    }

    template < class Derived > 
    inline bool
    PotentialUniqueDistTemplate< Derived >::
    _computeForce(Real3D& force, const Particle &p1, const Particle &p2, const real curDist) const {
      Real3D dist = p1.position() - p2.position();
      return _computeForce(force, dist, curDist);
    }

    template < class Derived > 
    inline bool
    PotentialUniqueDistTemplate< Derived >::
    _computeForce(Real3D& force, const Real3D& dist, const real curDist) const {
      real distSqr = dist.sqr();
      if (distSqr > cutoffSqr)
        return false;
      else {
        return derived_this()->_computeForceRaw(force, dist, distSqr, curDist);
      }
    }
    
  }
}

#endif
