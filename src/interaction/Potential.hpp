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
#ifndef _INTERACTION_POTENTIAL_HPP
#define _INTERACTION_POTENTIAL_HPP

#include "types.hpp"
#include "Real3D.hpp"
#include "RealND.hpp"
#include "Particle.hpp"
#include "logging.hpp"
#include "FixedPairList.hpp"
#include "FixedTripleList.hpp"
#include "FixedQuadrupleList.hpp"

namespace espressopp {
  namespace interaction {
    class Potential {
    public:
      virtual ~Potential() {};
      virtual real computeEnergy(const Particle &p1, const Particle &p2) const = 0;
      virtual real computeEnergy(const Real3D& dist) const = 0;
      virtual real computeEnergy(real dist) const = 0;
      virtual real computeEnergySqr(real distSqr) const = 0;

      virtual Real3D computeForce(const Particle &p1, const Particle &p2) const = 0;
      virtual Real3D computeForce(const Real3D& dist) const = 0;

      virtual void setCutoff(real _cutoff) = 0;
      virtual real getCutoff() const = 0;

      virtual void setColVarBondList(const shared_ptr<FixedPairList>& fpl) = 0;
      virtual shared_ptr<FixedPairList> getColVarBondList() const = 0;
      virtual void setColVarAngleList(const shared_ptr<FixedTripleList>& fpl) = 0;
      virtual shared_ptr<FixedTripleList> getColVarAngleList() const = 0;
      virtual void setColVarDihedList(const shared_ptr<FixedQuadrupleList>& fpl) = 0;
      virtual shared_ptr<FixedQuadrupleList> getColVarDihedList() const = 0;

      virtual void setColVar(const RealND& cv) = 0;
      virtual void setColVar(const Real3D& dist, const bc::BC& bc) = 0;
      virtual RealND getColVar() const = 0;

      virtual void computeColVarWeights(const Real3D& dist, const bc::BC& bc) = 0;

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
    class PotentialTemplate : public Potential {
    public:
      PotentialTemplate();
      virtual ~PotentialTemplate() {};

      // Implements the Potential virtual interface
      virtual real computeEnergy(const Particle &p1, const Particle &p2) const;
      virtual real computeEnergy(const Real3D& dist) const;
      virtual real computeEnergy(real dist) const;
      virtual real computeEnergySqr(real distsqr) const;
      virtual Real3D computeForce(const Particle &p1, const Particle &p2) const;
      virtual Real3D computeForce(const Real3D& dist) const;
      virtual void setCutoff(real _cutoff);
      virtual real getCutoff() const;
      virtual void setShift(real _shift);
      virtual real getShift() const;
      virtual real setAutoShift();
      void updateAutoShift();

      virtual void setColVarBondList(const shared_ptr<FixedPairList>& fpl);
      virtual shared_ptr<FixedPairList> getColVarBondList() const;
      virtual void setColVarAngleList(const shared_ptr<FixedTripleList>& fpl);
      virtual shared_ptr<FixedTripleList> getColVarAngleList() const;
      virtual void setColVarDihedList(const shared_ptr<FixedQuadrupleList>& fpl);
      virtual shared_ptr<FixedQuadrupleList> getColVarDihedList() const;

      virtual void setColVar(const RealND& cv);
      virtual void setColVar(const Real3D& dist, const bc::BC& bc);
      virtual RealND getColVar() const;

      virtual void computeColVarWeights(const Real3D& dist, const bc::BC& bc);

      // Implements the non-virtual interface
      // (used by e.g. the Interaction templates)
      real _computeEnergy(const Particle &p1, const Particle &p2) const;
      real _computeEnergy(const Real3D& dist) const;
      real _computeEnergy(real dist) const;
      real _computeEnergy(const Particle &p1, const Particle &p2, const Real3D& dist) const;
      real _computeEnergySqr(real distSqr) const;

      real _computeEnergyDeriv(const Particle &p1, const Particle &p2) const;

      bool _computeForce(Real3D& force,
			 const Particle &p1, const Particle &p2) const;
      bool _computeForce(Real3D& force,
			 const Real3D& dist) const;
      bool _computeForce(Real3D& force,
                         const Particle &p1, const Particle &p2, const Real3D& dist) const;

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
      // List of bonds that correlate with the bond potential
      shared_ptr<FixedPairList> colVarBondList;
      // List of angles that correlate with the bond potential
      shared_ptr<FixedTripleList> colVarAngleList;
      // List of dihedrals that correlate with the angle potential
      shared_ptr<FixedQuadrupleList> colVarDihedList;
      // Collective variables: first itself, then angles
      RealND colVar;

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

    //   real computeEnergy(const Real3D& dist) {
    // 	return _computeEnergy(dist);
    //   }
    // };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < class Derived >
    inline
    PotentialTemplate< Derived >::PotentialTemplate() : cutoff(infinity),
      cutoffSqr(infinity), shift(0.0), autoShift(false), colVar(1, 0.) {
    }

    // Shift/cutoff handling
    template < class Derived >
    inline void
    PotentialTemplate< Derived >::setCutoff(real _cutoff) {
      cutoff = _cutoff;
      cutoffSqr = cutoff*cutoff;
      LOG4ESPP_INFO(Derived::theLogger, " cutoff=" << cutoff);
      updateAutoShift();
    }

    template < class Derived >
    inline real
    PotentialTemplate< Derived >::getCutoff() const {
        return cutoff;
    }

    template < class Derived >
    inline void
    PotentialTemplate< Derived >::setShift(real _shift) {
      autoShift = false;
      shift = _shift;
      LOG4ESPP_INFO(Derived::theLogger, " (manual) shift=" << shift);
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
      LOG4ESPP_INFO(Derived::theLogger, " (auto) shift=" << shift);
      return shift;
    }

    template < class Derived >
    inline void
    PotentialTemplate< Derived >::
    updateAutoShift() {
      if (autoShift) setAutoShift();
    }

    // Pair list for collective variables
    template < class Derived >
    inline void
    PotentialTemplate< Derived >::
    setColVarBondList(const shared_ptr < FixedPairList >& _fpl) {
      colVarBondList = _fpl;
    }

    template < class Derived >
    inline shared_ptr < FixedPairList >
    PotentialTemplate< Derived >::
    getColVarBondList() const
    { return colVarBondList; }

    // Triple list for collective variables
    template < class Derived >
    inline void
    PotentialTemplate< Derived >::
    setColVarAngleList(const shared_ptr < FixedTripleList >& _fpl) {
      colVarAngleList = _fpl;
    }

    template < class Derived >
    inline shared_ptr < FixedTripleList >
    PotentialTemplate< Derived >::
    getColVarAngleList() const
    { return colVarAngleList; }

    // Quadruple list for collective variables
    template < class Derived >
    inline void
    PotentialTemplate< Derived >::
    setColVarDihedList(const shared_ptr < FixedQuadrupleList >& _fpl) {
      colVarDihedList = _fpl;
    }

    template < class Derived >
    inline shared_ptr < FixedQuadrupleList >
    PotentialTemplate< Derived >::
    getColVarDihedList() const
    { return colVarDihedList; }

    // Collective variables
    template < class Derived >
    inline void
    PotentialTemplate< Derived >::
    setColVar(const RealND& cv)
    { colVar = cv; }

    // Collective variables
    template < class Derived >
    inline void
    PotentialTemplate< Derived >::
    setColVar(const Real3D& dist, const bc::BC& bc)
    {
      // In general, do nothing
    }

    template < class Derived >
    inline RealND
    PotentialTemplate< Derived >::
    getColVar() const
    { return colVar; }

    template < class Derived >
    inline void
    PotentialTemplate< Derived >::
    computeColVarWeights(const Real3D& dist, const bc::BC& bc) {
        // in general do nothing
    }

    // Energy computation
    template < class Derived >
    inline real
    PotentialTemplate< Derived >::
    computeEnergy(const Particle& p1, const Particle& p2) const {
      Real3D dist = p1.position() - p2.position();
      return computeEnergy(dist);
    }

    template < class Derived >
    inline real
    PotentialTemplate< Derived >::
    computeEnergy(const Real3D& dist) const {
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
    _computeEnergy(const Particle& p1, const Particle& p2) const {
      Real3D dist = p1.position() - p2.position();
      return _computeEnergy(dist);
    }

    template < class Derived >
    inline real
    PotentialTemplate< Derived >::
    _computeEnergy(const Particle& p1, const Particle& p2, const Real3D& dist) const {
      return _computeEnergy(dist);
    }

    template < class Derived >
    inline real
    PotentialTemplate< Derived >::
    _computeEnergy(const Real3D& dist) const {
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
      else {
        real e = derived_this()->_computeEnergySqrRaw(distSqr) - shift;
        LOG4ESPP_TRACE(theLogger, "Epot(r*r=" << distSqr << ") = " << e);
        return e;
      }
    }

    //computation of derivative wrt TI lambda
    template < class Derived >
    inline real
    PotentialTemplate< Derived >::
    _computeEnergyDeriv(const Particle& p1, const Particle& p2) const {
      std::cout<<"Calculation of energy derivative wrt lambda not implemented for all potentials!"<<std::endl;
      return 0.0;
    }


    // Force computation
    template < class Derived >
    inline Real3D
    PotentialTemplate< Derived >::
    computeForce(const Particle& p1, const Particle& p2) const {
      Real3D dist = p1.position() - p2.position();
      return computeForce(dist);
    }

    template < class Derived >
    inline Real3D
    PotentialTemplate< Derived >::
    computeForce(const Real3D& dist) const {
      Real3D force;
      if(!_computeForce(force, dist))
        force = 0.0;
      return force;
    }

    template < class Derived >
    inline bool
    PotentialTemplate< Derived >::
    _computeForce(Real3D& force, const Particle &p1, const Particle &p2) const {
      Real3D dist = p1.position() - p2.position();
      return _computeForce(force, dist);
    }

    template < class Derived >
    inline bool
    PotentialTemplate< Derived >::
    _computeForce(Real3D& force, const Particle& p1, const Particle& p2, const Real3D& dist) const {
      return _computeForce(force, dist);
    }

    template < class Derived >
    inline bool
    PotentialTemplate< Derived >::
    _computeForce(Real3D& force, const Real3D& dist) const {
      real distSqr = dist.sqr();
      if (distSqr > cutoffSqr)
        return false;
      else {
        return derived_this()->_computeForceRaw(force, dist, distSqr);
      }
    }

  }
}

#endif
