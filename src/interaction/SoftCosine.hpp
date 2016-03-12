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
#ifndef _INTERACTION_SOFTCOSINE_HPP
#define _INTERACTION_SOFTCOSINE_HPP

#include "Potential.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
	the SoftCosine potential.

	\f[
  V(r) = A \left[ 1.0 + cos \left( \frac{\pi r}{r_c} \right) \right]
	\f]
    */
    class SoftCosine : public PotentialTemplate< SoftCosine > {
    private:
      real A;

    public:
      static void registerPython();

      SoftCosine() : A(0.0) {
	setShift(0.0);
	setCutoff(infinity);
        preset();
      }

      SoftCosine(real _A, real _cutoff, real _shift) : A(_A) {
	setShift(_shift);
	setCutoff(_cutoff);
        preset();
      }

      SoftCosine(real _A, real _cutoff) : A(_A) {	
	autoShift = false;
	setCutoff(_cutoff);
	setAutoShift(); 
        preset();
      }

      void preset() { }

      // Setter and getter
      void setA(real _A) { 
	A = _A; 
	updateAutoShift();
        preset();
      }

      real getA() const { return A; }

      real _computeEnergySqrRaw(real distSqr) const {
        real r = sqrt(distSqr);
	real energy = A * (1.0 + cos(M_PI * r / getCutoff()));
	return energy;
      }

      bool _computeForceRaw(Real3D& force,
                            const Real3D& dist,
                            real distSqr) const {
        real r = sqrt(distSqr);
        real rc = getCutoff();
        real ffactor = (A * M_PI) * sin(M_PI * r / rc) / (rc * r);
        force = dist * ffactor;
        return true;
      }
    };
    // provide pickle support
    struct SoftCosine_pickle : boost::python::pickle_suite
    {
      static
      boost::python::tuple
      getinitargs(SoftCosine const& pot)
      {
          real a;
          real rc;
          real sh;
          a=pot.getA();
          rc =pot.getCutoff();
          sh =pot.getShift();
          return boost::python::make_tuple(a, rc, sh);
      }
    };

  }
}

#endif
