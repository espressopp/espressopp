/*
  Copyright (C) 2012,2013,2016,2018
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
#include "RealND.hpp"
#include "Particle.hpp"
#include <cmath>
#include "logging.hpp"
#include "FixedPairList.hpp"
#include "FixedTripleList.hpp"
#include "FixedQuadrupleList.hpp"

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

      virtual void setColVarBondList(const shared_ptr<FixedPairList>& fpl) = 0;
      virtual shared_ptr<FixedPairList> getColVarBondList() const = 0;
      virtual void setColVarAngleList(const shared_ptr<FixedTripleList>& ftl) = 0;
      virtual shared_ptr<FixedTripleList> getColVarAngleList() const = 0;
      virtual void setColVarDihedList(const shared_ptr<FixedQuadrupleList>& fql) = 0;
      virtual shared_ptr<FixedQuadrupleList> getColVarDihedList() const = 0;

      virtual void setColVar(const RealND& cv) = 0;
      virtual void setColVar(const Real3D& dist21,
          const Real3D& dist32, const Real3D& dist43, const bc::BC& bc) = 0;
      virtual RealND getColVar() const = 0;

      virtual void computeColVarWeights(const Real3D& dist21,
          const Real3D& dist32, const Real3D& dist43,
          const bc::BC& bc) = 0;

      static real computePhi(const Real3D& dist21,
                             const Real3D& dist32,
                             const Real3D& dist43);
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

      virtual void setColVarBondList(const shared_ptr<FixedPairList>& fpl);
      virtual shared_ptr<FixedPairList> getColVarBondList() const;
      virtual void setColVarAngleList(const shared_ptr<FixedTripleList>& ftl);
      virtual shared_ptr<FixedTripleList> getColVarAngleList() const;
      virtual void setColVarDihedList(const shared_ptr<FixedQuadrupleList>& fql);
      virtual shared_ptr<FixedQuadrupleList> getColVarDihedList() const;

      virtual void setColVar(const RealND& cv);
      virtual void setColVar(const Real3D& dist21,
          const Real3D& dist32, const Real3D& dist43, const bc::BC& bc);
      virtual RealND getColVar() const;

      virtual void computeColVarWeights(const Real3D& dist21,
          const Real3D& dist32, const Real3D& dist43,
          const bc::BC& bc);

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
      // List of bonds that correlate with the bond potential
      shared_ptr<FixedPairList> colVarBondList;
      // List of angles that correlate with the bond potential
      shared_ptr<FixedTripleList> colVarAngleList;
      // List of dihedrals that correlate with the dihedral potential
      shared_ptr<FixedQuadrupleList> colVarDihedList;
      // Collective variables: first itself, then bonds, then angles
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
    DihedralPotentialTemplate< Derived >::
    DihedralPotentialTemplate()
      : cutoff(infinity), cutoffSqr(infinity), colVar(1, 0.)
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

    // Pair list for collective variables
    template < class Derived >
    inline void
    DihedralPotentialTemplate< Derived >::
    setColVarBondList(const shared_ptr < FixedPairList >& _fpl) {
      colVarBondList = _fpl;
    }

    template < class Derived >
    inline shared_ptr < FixedPairList >
    DihedralPotentialTemplate< Derived >::
    getColVarBondList() const
    { return colVarBondList; }

    // Pair list for collective variables
    template < class Derived >
    inline void
    DihedralPotentialTemplate< Derived >::
    setColVarAngleList(const shared_ptr < FixedTripleList >& _fpl) {
      colVarAngleList = _fpl;
    }

    template < class Derived >
    inline shared_ptr < FixedTripleList >
    DihedralPotentialTemplate< Derived >::
    getColVarAngleList() const
    { return colVarAngleList; }

    // Quadruple list for collective variables
    template < class Derived >
    inline void
    DihedralPotentialTemplate< Derived >::
    setColVarDihedList(const shared_ptr < FixedQuadrupleList >& _fql) {
      colVarDihedList = _fql;
    }

    template < class Derived >
    inline shared_ptr < FixedQuadrupleList >
    DihedralPotentialTemplate< Derived >::
    getColVarDihedList() const
    { return colVarDihedList; }

    // Collective variables
    template < class Derived >
    inline void
    DihedralPotentialTemplate< Derived >::
    setColVar(const RealND& cv)
    { colVar = cv; }

    // Collective variables
    template < class Derived >
    inline void
    DihedralPotentialTemplate< Derived >::
    setColVar(const Real3D& dist21, const Real3D& dist32,
        const Real3D& dist43, const bc::BC& bc)
    {
      // In general, do nothing
    }

    template < class Derived >
    inline RealND
    DihedralPotentialTemplate< Derived >::
    getColVar() const
    { return colVar; }

    template < class Derived >
    inline void
    DihedralPotentialTemplate< Derived >::
    computeColVarWeights(const Real3D& dist21, const Real3D& dist32,
        const Real3D& dist43, const bc::BC& bc) {
        // in general do nothing
    }

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
        real phi = computePhi(r21, r32, r43);
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
