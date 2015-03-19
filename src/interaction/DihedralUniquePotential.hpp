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
#ifndef _INTERACTION_DIHEDRALUNIQUEPOTENTIAL_HPP
#define _INTERACTION_DIHEDRALUNIQUEPOTENTIAL_HPP

#include "types.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include <cmath>
#include "logging.hpp"

namespace espressopp {
  namespace interaction {
    class DihedralUniquePotential {
    public:
      virtual real computeEnergy(const Real3D& r21,
                                 const Real3D& r32,
                                 const Real3D& r43, const real phi0) const = 0;

      virtual real computeEnergy(real phi, real phi0) const = 0;

      virtual void computeForce(Real3D& force1,
                                Real3D& force2,
                                Real3D& force3,
                                Real3D& force4,
                                const Real3D& r21,
                                const Real3D& r32,
                                const Real3D& r43, const real phi0) const = 0;

      virtual real computeForce(real phi, real phi0) const = 0;

      virtual void setCutoff(real _cutoff) = 0;
      virtual real getCutoff() const = 0;

      static void registerPython();
      static LOG4ESPP_DECL_LOGGER(theLogger);
    };

    template < class Derived >
    class DihedralUniquePotentialTemplate : public DihedralUniquePotential {
    public:
      DihedralUniquePotentialTemplate();

      // Implements the Potential virtual interface
      virtual real computeEnergy(const Real3D& r21,
                                 const Real3D& r32,
                                 const Real3D& r43, const real phi0) const;

      virtual real computeEnergy(real phi, real phi0) const;

      virtual void computeForce(Real3D& force1,
                                Real3D& force2,
                                Real3D& force3,
                                Real3D& force4,
                                const Real3D& r21,
                                const Real3D& r32,
                                const Real3D& r43, const real phi0) const;

      virtual real computeForce(real phi, real phi0) const; // used for generating tabular file

      virtual void setCutoff(real _cutoff);
      virtual real getCutoff() const;

      // Implements the non-virtual interface 
      // (used by e.g. the Interaction templates)
      real _computeEnergy(const Real3D& r21,
                          const Real3D& r32,
                          const Real3D& r43, const real phi0) const;


      real _computeEnergy(real phi, real phi0) const;

      void _computeForce(Real3D& force1,
                         Real3D& force2,
                         Real3D& force3,
                         Real3D& force4,
                         const Real3D& r21,
                         const Real3D& r32,
                         const Real3D& r43, const real phi0) const;


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
    DihedralUniquePotentialTemplate< Derived >::
    DihedralUniquePotentialTemplate() 
      : cutoff(infinity), cutoffSqr(infinity)
    {}

    // Shift/cutoff handling
    template < class Derived >
    inline void
    DihedralUniquePotentialTemplate< Derived >::
    setCutoff(real _cutoff) {
      cutoff = _cutoff;
      cutoffSqr = cutoff*cutoff;
    }

    template < class Derived >
    inline real
    DihedralUniquePotentialTemplate< Derived >::
    getCutoff() const
    { return cutoff; }

    // Energy computation
    template < class Derived >
    inline real
    DihedralUniquePotentialTemplate< Derived >::
    computeEnergy(const Real3D& r21,
                  const Real3D& r32,
                  const Real3D& r43, const real phi0) const {
        
      return _computeEnergy(r21, r32, r43, phi0);
      
      
      
    }

    template < class Derived >
    inline real
    DihedralUniquePotentialTemplate< Derived >::
    computeEnergy(real phi, real phi0) const {
      return _computeEnergy(phi, phi0);
    }
    
    template < class Derived >
    inline real
    DihedralUniquePotentialTemplate< Derived >::
    _computeEnergy(const Real3D& r21,
                   const Real3D& r32,
                   const Real3D& r43, const real phi0) const {
      
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
        
        return _computeEnergy(phi, phi0);
    }

    template < class Derived > 
    inline real
    DihedralUniquePotentialTemplate< Derived >::
    _computeEnergy(real phi, real phi0) const {
      return derived_this()->_computeEnergyRaw(phi, phi0);
    }
    
    // Force computation
    template < class Derived >
    inline void
    DihedralUniquePotentialTemplate< Derived >::
    computeForce(Real3D& force1,
                 Real3D& force2,
                 Real3D& force3,
                 Real3D& force4,
                 const Real3D& r21,
                 const Real3D& r32,
                 const Real3D& r43, const real phi0) const {
      
        _computeForce(force1, force2, force3, force4, r21, r32, r43, phi0);
    }

    template < class Derived >
    inline void
    DihedralUniquePotentialTemplate< Derived >::
    _computeForce(Real3D& force1,
                  Real3D& force2,
                  Real3D& force3,
                  Real3D& force4,
                  const Real3D& r21,
                  const Real3D& r32,
                  const Real3D& r43, const real phi0) const {
      derived_this()->_computeForceRaw(force1, force2, force3, force4, r21, r32, r43, phi0);
    }
    
    // used for generating tabular file
    template < class Derived >
    inline real
    DihedralUniquePotentialTemplate< Derived >::
    computeForce(real phi, real phi0) const {
      return derived_this()->_computeForceRaw(phi, phi0);
    }
  }
}

#endif
