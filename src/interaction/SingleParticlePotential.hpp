/*
  Copyright (C) 2014
      Pierre de Buyl
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
#ifndef _INTERACTION_SINGLEPARTICLEPOTENTIAL_HPP
#define _INTERACTION_SINGLEPARTICLEPOTENTIAL_HPP

#include "types.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "logging.hpp"

namespace espresso {
  namespace interaction {

    /** This class is used to define single-particle interactions, typically
        used for external forces on the system.
     */
    class SingleParticlePotential {
    public:
      virtual real computeEnergy(const Real3D& position, const bc::BC& bc) const = 0;
      virtual Real3D computeForce(const Real3D& position, const bc::BC& bc) const = 0;
      virtual real getMaxCutoff() = 0;

      static void registerPython();
      static LOG4ESPP_DECL_LOGGER(theLogger);
    };

    //    template < class Derived, enum PotentialType = Default >
    /** Provides a template for the simple implementation of a
        position-dependent potential.
    */
    template < class Derived >
    class SingleParticlePotentialTemplate : public SingleParticlePotential {
    public:
      SingleParticlePotentialTemplate();
      virtual ~SingleParticlePotentialTemplate() {};

      // Implements the SingleParticlePotential virtual interface
      virtual real computeEnergy(const Real3D& position, const bc::BC& bc) const;
      virtual Real3D computeForce(const Real3D& position, const bc::BC& bc) const;

      // Implements the non-virtual interface
      // (used by e.g. the Interaction templates)
      real _computeEnergy(const Real3D& position, const bc::BC& bc) const;
      bool _computeForce(Real3D& force, const Real3D& position, const bc::BC& bc) const;

    protected:
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
    SingleParticlePotentialTemplate< Derived >::SingleParticlePotentialTemplate() {}

    // Energy computation
    template < class Derived >
    inline real
    SingleParticlePotentialTemplate< Derived >::
    computeEnergy(const Real3D& position, const bc::BC& bc) const {
        return _computeEnergy(position, bc);
      }

    template < class Derived >
    inline real
    SingleParticlePotentialTemplate< Derived >::
    _computeEnergy(const Real3D& position, const bc::BC& bc) const {
      return derived_this()->_computeEnergyRaw(position, bc);
    }

    // Force computation
    template < class Derived >
    inline Real3D
    SingleParticlePotentialTemplate< Derived >::
    computeForce(const Real3D& position, const bc::BC& bc) const {
      Real3D force;
      if(!_computeForce(force, position, bc)) {
        force = 0.0;
      }
      return force;
    }

    template < class Derived >
    inline bool
    SingleParticlePotentialTemplate< Derived >::
    _computeForce(Real3D& force, const Real3D &position, const bc::BC& bc) const {
      return derived_this()->_computeForceRaw(force, position, bc);
    }

  }
}

#endif
