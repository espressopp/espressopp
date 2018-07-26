/*
  Copyright (C) 2012-2018
      Max Planck Institute for Polymer Research
  Copyright (C) 2008-2011
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
#ifndef _INTERACTION_INTERACTION_HPP
#define _INTERACTION_INTERACTION_HPP

#include "types.hpp"
#include "logging.hpp"
#include "esutil/ESPPIterator.hpp"

namespace espressopp {
  namespace interaction {

    enum bondTypes {unused, Nonbonded, Single, Pair, Angular, Dihedral, NonbondedSlow};

    /** Interaction base class. */

    class Interaction {

    public:
      virtual ~Interaction() {};
      virtual void addForces() = 0;
      virtual real computeEnergy() = 0;
      virtual real computeEnergyDeriv() = 0;
      virtual real computeEnergyAA() = 0;
      virtual real computeEnergyCG() = 0;
      virtual real computeEnergyAA(int atomtype) = 0;
      virtual real computeEnergyCG(int atomtype) = 0;
      virtual real computeVirial() = 0;
      virtual void computeVirialTensor(Tensor& w) = 0;
      virtual void computeVirialX(std::vector<real> &p_xx_total, int bins) = 0;
      // this should compute the virial locally around a surface which crosses the box at
      // z (according to the method of Irving and Kirkwood)
      virtual void computeVirialTensor(Tensor& w, real z) = 0;
      // the same Irving - Kirkwood method, but Z direction is divided by n planes
      virtual void computeVirialTensor(Tensor *w, int n) = 0;

      /** This method returns the maximal cutoff defined for one type pair. */
      virtual real getMaxCutoff() = 0;
      virtual int bondType() = 0;

      static void registerPython();

    protected:
      /** Logger */
      static LOG4ESPP_DECL_LOGGER(theLogger);
    };

    struct InteractionList
      : public std::vector< shared_ptr< Interaction > > {
      typedef esutil::ESPPIterator< std::vector< Interaction > > Iterator;
    };


  }
}

#endif
