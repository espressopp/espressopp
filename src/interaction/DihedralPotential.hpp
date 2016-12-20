/*
  Copyright (C) 2012,2013,2016
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
#ifndef _INTERACTION_DIHEDRALPOTENTIAL_HPP
#define _INTERACTION_DIHEDRALPOTENTIAL_HPP

#include "types.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include <cmath>
#include "logging.hpp"

namespace espressopp {
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

      virtual real computeForce(real phi) const = 0;

      virtual void setCutoff(real _cutoff) = 0;
      virtual real getCutoff() const = 0;

      static void registerPython();
      static LOG4ESPP_DECL_LOGGER(theLogger);
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

      virtual real computeForce(real phi) const; // used for generating tabular file

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
        
      return _computeEnergy(dist21, dist32, dist43);
      
      
      
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
    _computeEnergy(const Real3D& r21,
                   const Real3D& r32,
                   const Real3D& r43) const {
                       
                       
/*                       
        // compute phi
        real dist21_sqr = dist21 * dist21;
        real dist32_sqr = dist32 * dist32;
        real dist43_sqr = dist43 * dist43;
        real dist21_magn = sqrt(dist21_sqr);
        real dist32_magn = sqrt(dist32_sqr);
        real dist43_magn = sqrt(dist43_sqr);
        
        // cos0
        real sb1 = 1.0 / dist21_sqr;
        real sb2 = 1.0 / dist32_sqr;
        real sb3 = 1.0 / dist43_sqr;
        real rb1 = sqrt(sb1);
        real rb3 = sqrt(sb3);
        real c0 = dist21 * dist43 * rb1 * rb3;
        
        
        // 1st and 2nd angle
        real ctmp = dist21 * dist32;
        real r12c1 = 1.0 / (dist21_magn * dist32_magn);
        real c1mag = ctmp * r12c1;
        
        ctmp = (-1.0 * dist32) * dist43;
        real r12c2 = 1.0 / (dist32_magn * dist43_magn);
        real c2mag = ctmp * r12c2;
        
        
        //cos and sin of 2 angles and final cos
        real sin2 = 1.0 - c1mag * c1mag;
        if (sin2 < 0) sin2 = 0.0;
        real sc1 = sqrt(sin2);
        sc1 = 1.0 / sc1;
        
        sin2 = 1.0 - c2mag * c2mag;
        if (sin2 < 0) sin2 = 0.0;
        real sc2 = sqrt(sin2);
        sc2 = 1.0 / sc2;
        
        real s1 = sc1 * sc1;
        real s2 = sc2 * sc2;
        real s12 = sc1 * sc2;
        real c = (c0 + c1mag * c2mag) * s12;
        
        Real3D cc = dist21.cross(dist32);
        real cmag = sqrt(cc * cc);
        real dx = cc * dist43 / cmag / dist43_magn;
        
        if (c > 1.0) c = 1.0;
        else if (c < -1.0) c = -1.0;
        
        // phi
        real phi = acos(c);
        if (dx < 0.0) phi *= -1.0;
 */
      
        Real3D rijjk = r21.cross(r32); // [r21 x r32]
        Real3D rjkkn = r32.cross(r43); // [r32 x r43]
        
        real rijjk_sqr = rijjk.sqr();
        real rjkkn_sqr = rjkkn.sqr();
        
        real rijjk_abs = sqrt(rijjk_sqr);
        real rjkkn_abs = sqrt(rjkkn_sqr);
        
        real inv_rijjk = 1.0 / rijjk_abs;
        real inv_rjkkn = 1.0 / rjkkn_abs;
        
        // cosine between planes
        real cos_phi = (rijjk * rjkkn) * (inv_rijjk * inv_rjkkn);
        if (cos_phi > 1.0) cos_phi = 1.0;
        else if (cos_phi < -1.0) cos_phi = -1.0;
        
        real phi = acos(cos_phi);
        //get sign of phi
        //positive if (rij x rjk) x (rjk x rkn) is in the same direction as rjk, negative otherwise (see DLPOLY manual)
        Real3D rcross = rijjk.cross(rjkkn); //(rij x rjk) x (rjk x rkn)
        real signcheck = rcross * r32;
        if (signcheck < 0.0) phi *= -1.0;
        
        return _computeEnergy(phi);
      
      
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
    
    // used for generating tabular file
    template < class Derived >
    inline real
    DihedralPotentialTemplate< Derived >::
    computeForce(real phi) const {
      return derived_this()->_computeForceRaw(phi);
    }
  }
}

#endif
