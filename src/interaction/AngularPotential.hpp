/*
  Copyright (C) 2012,2013,2018
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
#ifndef _INTERACTION_ANGULARPOTENTIAL_HPP
#define _INTERACTION_ANGULARPOTENTIAL_HPP

#include "types.hpp"
#include "Real3D.hpp"
#include "RealND.hpp"
#include "Particle.hpp"
#include <cmath>
#include "logging.hpp"
#include "FixedPairList.hpp"
#include "FixedTripleList.hpp"
#include "FixedQuadrupleList.hpp"

namespace espressopp {
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

      virtual void setColVarBondList(const shared_ptr<FixedPairList>& fpl) = 0;
      virtual shared_ptr<FixedPairList> getColVarBondList() const = 0;
      virtual void setColVarAngleList(const shared_ptr<FixedTripleList>& fpl) = 0;
      virtual shared_ptr<FixedTripleList> getColVarAngleList() const = 0;
      virtual void setColVarDihedList(const shared_ptr<FixedQuadrupleList>& fpl) = 0;
      virtual shared_ptr<FixedQuadrupleList> getColVarDihedList() const = 0;

      virtual void setColVar(const RealND& cv) = 0;
      virtual void setColVar(const Real3D& dist12,
          const Real3D& dist32, const bc::BC& bc) = 0;
      virtual RealND getColVar() const = 0;

      virtual void computeColVarWeights(const Real3D& dist12,
          const Real3D& dist32, const bc::BC& bc) = 0;

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

      virtual void setColVarBondList(const shared_ptr<FixedPairList>& fpl);
      virtual shared_ptr<FixedPairList> getColVarBondList() const;
      virtual void setColVarAngleList(const shared_ptr<FixedTripleList>& fpl);
      virtual shared_ptr<FixedTripleList> getColVarAngleList() const;
      virtual void setColVarDihedList(const shared_ptr<FixedQuadrupleList>& fpl);
      virtual shared_ptr<FixedQuadrupleList> getColVarDihedList() const;

      virtual void setColVar(const RealND& cv);
      virtual void setColVar(const Real3D& dist12,
          const Real3D& dist32, const bc::BC& bc);
      virtual RealND getColVar() const;

      virtual void computeColVarWeights(const Real3D& dist12,
          const Real3D& dist32, const bc::BC& bc);

      // Implements the non-virtual interface
      // (used by e.g. the Interaction templates)
      real _computeEnergy(const Particle &p1, const Particle &p2, const Particle &p3) const;
      real _computeEnergy(const Real3D& dist12, const Real3D& dist32) const;
      real _computeEnergy(real theta) const;

      void _computeForce(Real3D& force12, Real3D& force32,
			 const Particle &p1, const Particle &p2, const Particle &p3) const;

      bool _computeForce(Real3D& force12, Real3D& force32,
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
      // List of bonds that correlate with the angle potential
      shared_ptr<FixedPairList> colVarBondList;
      // List of angles that correlate with the bond potential
      shared_ptr<FixedTripleList> colVarAngleList;
      // List of dihedrals that correlate with the dihedral potential
      shared_ptr<FixedQuadrupleList> colVarDihedList;
      // Collective variables: first itself, then bonds and dihedrals
      RealND colVar;

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
      : cutoff(infinity), cutoffSqr(infinity), colVar(1, 0.) {   }

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

    // Pair list for collective variables
    template < class Derived >
    inline void
    AngularPotentialTemplate< Derived >::
    setColVarBondList(const shared_ptr < FixedPairList >& _fpl) {
      colVarBondList = _fpl;
    }

    template < class Derived >
    inline shared_ptr < FixedPairList >
    AngularPotentialTemplate< Derived >::
    getColVarBondList() const
    { return colVarBondList; }

    // Triple list for collective variables
    template < class Derived >
    inline void
    AngularPotentialTemplate< Derived >::
    setColVarAngleList(const shared_ptr < FixedTripleList >& _fpl) {
      colVarAngleList = _fpl;
    }

    template < class Derived >
    inline shared_ptr < FixedTripleList >
    AngularPotentialTemplate< Derived >::
    getColVarAngleList() const
    { return colVarAngleList; }

    // Quadruple list for collective variables
    template < class Derived >
    inline void
    AngularPotentialTemplate< Derived >::
    setColVarDihedList(const shared_ptr < FixedQuadrupleList >& _fpl) {
      colVarDihedList = _fpl;
    }

    template < class Derived >
    inline shared_ptr < FixedQuadrupleList >
    AngularPotentialTemplate< Derived >::
    getColVarDihedList() const
    { return colVarDihedList; }

    // Collective variables
    template < class Derived >
    inline void
    AngularPotentialTemplate< Derived >::
    setColVar(const RealND& cv)
    { colVar = cv; }

    // Collective variables
    template < class Derived >
    inline void
    AngularPotentialTemplate< Derived >::
    setColVar(const Real3D& dist12,
        const Real3D& dist32, const bc::BC& bc)
    {
      // In general, do nothing
    }

    template < class Derived >
    inline RealND
    AngularPotentialTemplate< Derived >::
    getColVar() const
    { return colVar; }

    template < class Derived >
    inline void
    AngularPotentialTemplate< Derived >::
    computeColVarWeights(const Real3D& dist12,
        const Real3D& dist32, const bc::BC& bc) {
        // in general do nothing
    }

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
      if (dist12_sqr >= cutoffSqr || dist32_sqr >= cutoffSqr )
        return 0.0;
      else{
        real cos_theta = dist12 * dist32 / (sqrt(dist12_sqr) * sqrt(dist32_sqr));
        return _computeEnergy(acos(cos_theta));
      }
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
    inline bool
    AngularPotentialTemplate< Derived >::
    _computeForce(Real3D& force12,
                  Real3D& force32,
                  const Real3D& dist12,
                  const Real3D& dist32) const {

      return derived_this()->_computeForceRaw(force12, force32, dist12, dist32);
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
