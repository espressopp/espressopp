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
#ifndef _INTERACTION_POTENTIALVSPHEREPAIR_HPP
#define _INTERACTION_POTENTIALVSPHEREPAIR_HPP

#include "types.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "logging.hpp"

namespace espressopp {
  namespace interaction {
    class PotentialVSpherePair {
    public:
      virtual ~PotentialVSpherePair() {};
      virtual real computeEnergy(const Particle &p1, const Particle &p2) const = 0;
      virtual real computeEnergy(const Real3D& dist, real& sigmaij) const = 0;
      virtual real computeEnergy(real dist, real sigmaij) const = 0;
      virtual real computeEnergySqr(real distSqr, real sigmaij) const = 0;
      
      virtual python::list computeForce(const Particle &p1, const Particle &p2) const = 0;
      virtual python::list computeForce(const Real3D& dist, const real& sigmai, const real& sigmaj) const = 0;

      virtual void setCutoff(real _cutoff) = 0;
      virtual real getCutoff() const = 0;

      virtual void setShift(real _shift) = 0;
      virtual real getShift() const = 0;
      virtual real setAutoShift() = 0;

      static void registerPython();

      static LOG4ESPP_DECL_LOGGER(theLogger);
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
    class PotentialVSpherePairTemplate : public PotentialVSpherePair {
    public:
      PotentialVSpherePairTemplate();
      virtual ~PotentialVSpherePairTemplate() {};

      // Implements the Potential virtual interface
      virtual real computeEnergy(const Particle &p1, const Particle &p2) const;
      virtual real computeEnergy(const Real3D& dist, real &sigmaij) const;
      virtual real computeEnergy(real dist, real sigmaij) const;
      virtual real computeEnergySqr(real distsqr, real sigmaij) const;
      virtual python::list computeForce(const Particle &p1, const Particle &p2) const;
      virtual python::list computeForce(const Real3D& dist, const real& sigmai, const real& sigmaj) const;
      virtual void setCutoff(real _cutoff);
      virtual real getCutoff() const;
      virtual void setShift(real _shift);
      virtual real getShift() const;
      virtual real setAutoShift();
      void updateAutoShift();

      // Implements the non-virtual interface 
      // (used by e.g. the Interaction templates)
      real _computeEnergy(const Particle &p1, const Particle &p2) const;
      real _computeEnergy(const Real3D& dist, real& sigmaij) const;
      real _computeEnergy(real dist, real sigmaij) const;
      real _computeEnergySqr(real distSqr, real sigmaij) const;

      bool _computeForce(Real3D& force, real& sigiforce, real& sigjforce,
			 const Particle &p1, const Particle &p2) const;
      bool _computeForce(Real3D& force, real& sigiforce, real& sigjforce,
	        const Real3D& dist, const real& sigmai, const real& sigmaj) const;
      
      //bool _computeForce(CellList realcells) const;
      
      // Requires the following non-virtual interface in Derived
      // real _computeEnergySqrRaw(real distSqr) const;
      // bool _computeForceRaw(const Real3D& dist, const Real3D& force) const;


      // void _computeForce(const Particle &p1, const Particle &p2, 
      //                    Real3D& force) const {
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
    // class PotentialVSpherePairTemplate< Derived, NoScalarDistance >
    //   : public PotentialVSpherePairTemplate< Derived, Default > {
    // public:
    //   typedef PotentialVSpherePairTemplate< Derived, NoScalarDistance > Super;
    //   using PotentialVSpherePairTemplate< Derived >::_computeEnergy;
    //   using PotentialVSpherePairTemplate< Derived >::computeEnergy;

    //   real _computeEnergySqr(real distsqr) const {
    // 	throw Potential::BadCall();
    //   }

    //   real computeEnergy(const Real3D& dist) {
    // 	return _computeEnergy(dist);
    //   }
    // };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < class Derived > 
    inline
    PotentialVSpherePairTemplate< Derived >::PotentialVSpherePairTemplate() : cutoff(infinity), cutoffSqr(infinity), shift(0.0), autoShift(false){
    }

    // Shift/cutoff handling
    template < class Derived > 
    inline void 
    PotentialVSpherePairTemplate< Derived >::setCutoff(real _cutoff) {
      cutoff = _cutoff;
      cutoffSqr = cutoff*cutoff;
      updateAutoShift();
    }

    template < class Derived > 
    inline real 
    PotentialVSpherePairTemplate< Derived >::getCutoff() const {
        return cutoff;
    }

    template < class Derived > 
    inline void 
    PotentialVSpherePairTemplate< Derived >::setShift(real _shift) {
      autoShift = false; 
      shift = _shift; 
    }

    template < class Derived > 
    inline real 
    PotentialVSpherePairTemplate< Derived >::
    getShift() const 
    { return shift; }

    template < class Derived > 
    inline real 
    PotentialVSpherePairTemplate< Derived >::
    setAutoShift() {
      autoShift = true;
      if (cutoffSqr == infinity) {
	    shift = 0.0;
      } else {
        // WARNING setAutoShift may not make sense in PotentialVShere, using dummy for sigmaij
    	real dummy = 0.0;
	    shift = derived_this()->_computeEnergySqrRaw(cutoffSqr, dummy);
      }
      return shift;
    }

    template < class Derived > 
    inline void 
    PotentialVSpherePairTemplate< Derived >::
    updateAutoShift() {
      if (autoShift) setAutoShift();
    }

    // Energy computation
    template < class Derived > 
    inline real 
    PotentialVSpherePairTemplate< Derived >::
    computeEnergy(const Particle& p1, const Particle& p2) const {
      Real3D dist = p1.position() - p2.position();
      real r1 = p1.radius();
      real r2 = p2.radius();
      real sigmaij = r1*r1 + r2*r2;
      return computeEnergy(dist, sigmaij);
    }

    template < class Derived > 
    inline real 
    PotentialVSpherePairTemplate< Derived >::
    computeEnergy(const Real3D& dist, real& sigmaij) const {
	return computeEnergySqr(dist.sqr(), sigmaij);
      }

    template < class Derived > 
    inline real 
    PotentialVSpherePairTemplate< Derived >::
    computeEnergy(real dist, real sigmaij) const {
      return computeEnergySqr(dist*dist, sigmaij);
    }
    
    template < class Derived > 
    inline real 
    PotentialVSpherePairTemplate< Derived >::
    computeEnergySqr(real distsqr, real sigmaij) const {
      return _computeEnergySqr(distsqr, sigmaij);
    }

    template < class Derived > 
    inline real 
    PotentialVSpherePairTemplate< Derived >::
    _computeEnergy(const Particle& p1, const Particle& p2) const {
      Real3D dist = p1.position() - p2.position();
      real r1 = p1.radius();
      real r2 = p2.radius();
      real sigmaij = r1*r1 + r2*r2;
      return _computeEnergy(dist, sigmaij);
    }

    template < class Derived > 
    inline real 
    PotentialVSpherePairTemplate< Derived >::
    _computeEnergy(const Real3D& dist, real& sigmaij) const {
      return _computeEnergySqr(dist.sqr(), sigmaij);
    }

    template < class Derived > 
    inline real 
    PotentialVSpherePairTemplate< Derived >::
    _computeEnergy(real dist, real sigmaij) const {
      return _computeEnergySqr(dist*dist, sigmaij);
    }

    template < class Derived > 
    inline real
    PotentialVSpherePairTemplate< Derived >::
    _computeEnergySqr(real distSqr, real sigmaij) const {
      if (distSqr > cutoffSqr) 
        return 0.0;
      else {
        real e = derived_this()->_computeEnergySqrRaw(distSqr, sigmaij) - shift;
        LOG4ESPP_TRACE(theLogger, "Epot(r*r=" << distSqr << ") = " << e);
        return e;
      }
    }
    
    // Force computation
    template < class Derived > 
    inline python::list
    PotentialVSpherePairTemplate< Derived >::
    computeForce(const Particle& p1, const Particle& p2) const {
      Real3D fr;
      real fsi, fsj;
      python::list pl;
      bool ret;
      ret = _computeForce(fr, fsi, fsj, p1, p2);
      if (ret) {
        pl.append(fr[0]);
        pl.append(fr[1]);
        pl.append(fr[2]);
        pl.append(fsi);
        pl.append(fsj);
      }
      return pl;
    }

    template < class Derived > 
    inline python::list
    PotentialVSpherePairTemplate< Derived >::
    computeForce(const Real3D& dist, const real& sigmai, const real& sigmaj) const {
      Real3D fr;
      real fsi, fsj;
      python::list pl;
      bool ret;
      ret = _computeForce(fr, fsi, fsj, dist, sigmai, sigmaj);
      if (ret) {
        pl.append(fr[0]);
        pl.append(fr[1]);
        pl.append(fr[2]);
        pl.append(fsi);
        pl.append(fsj);
      }
      return pl;
    }

    template < class Derived > 
    inline bool
    PotentialVSpherePairTemplate< Derived >::
    _computeForce(Real3D& force, real& fsi, real& fsj, const Particle &p1, const Particle &p2) const {
      Real3D dist = p1.position() - p2.position();
      real sigmai = p1.radius();
      real sigmaj = p2.radius();
      return _computeForce(force, fsi, fsj, dist, sigmai, sigmaj);
    }


    template < class Derived > 
    inline bool
    PotentialVSpherePairTemplate< Derived >::
    _computeForce(Real3D& force, real& fsi, real& fsj, const Real3D& dist, const real& sigmai, const real& sigmaj) const {
      real distSqr = dist.sqr();
      if (distSqr > cutoffSqr)
        return false;
      else {
        return derived_this()->_computeForceRaw(force, fsi, fsj, dist, distSqr, sigmai, sigmaj);
      }
    }
    

  }
}

#endif
