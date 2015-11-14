/*
  Copyright (C) 2012,2013,2015
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
#ifndef _INTERACTION_COULOMBTRUNCATEDUNIQUECHARGE_HPP
#define _INTERACTION_COULOMBTRUNCATEDUNIQUECHARGE_HPP

#include "Potential.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    /* This class provides methods to compute forces and energies of
	the truncated Coulomb potential.
    */
    class CoulombTruncatedUniqueCharge : public PotentialTemplate< CoulombTruncatedUniqueCharge > {
    private:
      real qq;

    public:
      static void registerPython();

      CoulombTruncatedUniqueCharge()
	: qq(0.0) {
	setShift(0.0);
	setCutoff(infinity);
      }

      CoulombTruncatedUniqueCharge(real _qq,
		   real _cutoff, real _shift)
	: qq(_qq) {
	setShift(_shift);
	setCutoff(_cutoff);
      }

      CoulombTruncatedUniqueCharge(real _qq,
		   real _cutoff)
	: qq(_qq)
      {
	autoShift = false;
	setCutoff(_cutoff);
	setAutoShift();
      }

      // Setter and getter
      void setQQ(real _qq) {
	qq = _qq;
	updateAutoShift();
      }
      real getQQ() const { return qq; }

      real _computeEnergySqrRaw(real distSqr) const {
	real energy = qq / sqrt(distSqr);
	return energy;
      }

      bool _computeForceRaw(Real3D& force,
                            const Real3D& dist,
                            real distSqr) const {

        real ffactor;
	ffactor = qq / pow(sqrt(distSqr), 3);
        force = dist * ffactor;
        return true;
      }

    };
    // provide pickle support
    struct CoulombTruncatedUniqueCharge_pickle : boost::python::pickle_suite
    {
      static
      boost::python::tuple
      getinitargs(CoulombTruncatedUniqueCharge const& pot)
      {
    	  real q2;
          real rc;
          real sh;
          q2 =pot.getQQ();
          rc =pot.getCutoff();
          sh =pot.getShift();
          return boost::python::make_tuple(q2, rc, sh);
      }
    };
  }
}

#endif
