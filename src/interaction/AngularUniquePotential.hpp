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
#ifndef _INTERACTION_ANGULARUNIQUEPOTENTIAL_HPP
#define _INTERACTION_ANGULARUNIQUEPOTENTIAL_HPP

#include "types.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include <cmath>
#include "logging.hpp"

namespace espressopp {
  namespace interaction {
    class AngularUniquePotential {
    public:
      virtual real computeEnergy(const Particle &p1, const Particle &p2,
                const Particle &p3, real theta0) const = 0;
      virtual real computeEnergy(const Real3D& r12, const Real3D& r32,
                real theta0) const = 0;
      virtual real computeEnergy(real theta, real theta0) const = 0;

      virtual void computeForce(Real3D& force12, Real3D& force32, const Particle &p1,
                const Particle &p2, const Particle &p3, real theta0) const = 0;
      virtual void computeForce(Real3D& force12, Real3D& force32, const Real3D& r12,
                const Real3D& r32, real theta0) const = 0;
      virtual real computeForce(real theta, real theta0) const = 0; // used for generating tabular file


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
      virtual real computeEnergy(const Particle &p1, const Particle &p2,
                const Particle &p3, real theta0) const;
      virtual real computeEnergy(const Real3D& r12,
                const Real3D& r32, real theta0) const;
      virtual real computeEnergy(const real theta, real theta0) const;

      virtual void computeForce(Real3D& force12, Real3D& force32, const Particle &p1,
                const Particle &p2, const Particle &p3, real theta0) const;
      virtual void computeForce(Real3D& force12, Real3D& force32, const Real3D& r12,
                const Real3D& r32, real cos0) const;
      virtual real computeForce(real theta, real cos0) const; // used for generating tabular file

      virtual void setCutoff(real _cutoff);
      virtual real getCutoff() const;

      // Implements the non-virtual interface 
      // (used by e.g. the Interaction templates)
      real _computeEnergy(const Particle &p1, const Particle &p2,
                const Particle &p3, real theta0) const;
      real _computeEnergy(const Real3D& r12, const Real3D& r32, real theta0) const;
      real _computeEnergy(real theta, real theta0) const;

      void _computeForce(Real3D& force12, Real3D& force32, const Particle &p1,
                const Particle &p2, const Particle &p3, real theta0) const;
      
      bool _computeForce(Real3D& force12, Real3D& force32,
			 const Real3D& r12, const Real3D& r32, real theta0) const;

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
    computeEnergy(const Particle &p1, const Particle &p2, const Particle &p3, real theta0) const {
      Real3D r12 = p1.position() - p2.position();
      Real3D r32 = p3.position() - p2.position();
      return computeEnergy(r12, r32, theta0);
    }

    template < class Derived >
    inline real
    AngularUniquePotentialTemplate< Derived >::
    computeEnergy(const Real3D& r12, const Real3D& r32, real theta0) const {
      real dist12Sqr = r12.sqr();
      real dist32Sqr = r32.sqr();
      real cos_theta = r12 * r32 / (sqrt(dist12Sqr) * sqrt(dist32Sqr));
      return computeEnergy(acos(cos_theta), theta0);
    }

    template < class Derived >
    inline real 
    AngularUniquePotentialTemplate< Derived >::
    computeEnergy(real theta, real theta0) const {
      return _computeEnergy(theta, theta0); // a bug was here (it was: return computeEnergy(theta);)
    }
    
    template < class Derived > 
    inline real 
    AngularUniquePotentialTemplate< Derived >::
    _computeEnergy(const Particle &p1, const Particle &p2, const Particle &p3, real theta0) const {
      Real3D r12 = p1.position() - p2.position();
      Real3D r32 = p3.position() - p2.position();
      return _computeEnergy(r12, r32, theta0);
    }

    template < class Derived >
    inline real
    AngularUniquePotentialTemplate< Derived >::
    _computeEnergy(const Real3D& r12, const Real3D& r32, real theta0) const {
      real dist12_sqr = r12.sqr();
      real dist32_sqr = r32.sqr();
      if (dist12_sqr >= cutoffSqr || dist32_sqr >= cutoffSqr )
        return 0.0;
      else{
        real cos_theta = r12 * r32 / (sqrt(dist12_sqr) * sqrt(dist32_sqr));
        return _computeEnergy(acos(cos_theta), theta0);
      }
    }

    template < class Derived > 
    inline real
    AngularUniquePotentialTemplate< Derived >::
    _computeEnergy(real theta, real theta0) const {
      return derived_this()->_computeEnergyRaw(theta, theta0);
    }
    
    // Force computation
    template < class Derived >
    inline void
    AngularUniquePotentialTemplate< Derived >::
    computeForce(Real3D& force12,
                 Real3D& force32,
                 const Particle &p1, const Particle &p2,
                 const Particle &p3, real theta0) const {
      Real3D r12 = p1.position() - p2.position();
      Real3D r32 = p3.position() - p2.position();
      _computeForce(force12, force32, r12, r32, theta0);
    }

    template < class Derived >
    inline void
    AngularUniquePotentialTemplate< Derived >::
    computeForce(Real3D& force12,
                 Real3D& force32,
                 const Real3D& r12,
                 const Real3D& r32, real theta0) const {
      _computeForce(force12, force32, r12, r32, theta0);
    }

    template < class Derived >
    inline void
    AngularUniquePotentialTemplate< Derived >::
    _computeForce(Real3D& force12,
                  Real3D& force32,
                  const Particle &p1, const Particle &p2,
                  const Particle &p3, real theta0) const {
      Real3D r12 = p1.position() - p2.position();
      Real3D r32 = p3.position() - p2.position();
      _computeForce(force12, force32, r12, r32, theta0);
    }

    template < class Derived >
    inline bool
    AngularUniquePotentialTemplate< Derived >::
    _computeForce(Real3D& force12,
                  Real3D& force32,
                  const Real3D& r12,
                  const Real3D& r32, real theta0) const {
      
      return derived_this()->_computeForceRaw(force12, force32, r12, r32, theta0);
    }
    
    // used for generating tabular angular potential
    template < class Derived >
    inline real
    AngularUniquePotentialTemplate< Derived >::
    computeForce(real theta, real theta0) const {
      return derived_this()->_computeForceRaw(theta, theta0);
    }
    
  }
}

#endif
