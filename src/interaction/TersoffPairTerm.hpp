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
#ifndef _INTERACTION_TERSOFFPAIRTERM_HPP
#define _INTERACTION_TERSOFFPAIRTERM_HPP

#include "Potential.hpp"

#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

#include <cmath>

#ifndef M_PIl
#define M_PIl 3.1415926535897932384626433832795029L
#endif

namespace espressopp {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
	2 body term of Tersoff potential.
    */
    class TersoffPairTerm : public PotentialTemplate< TersoffPairTerm > {
    private:
      // 2 body
      real A, lambda1, R, D;
      real Pi_2D;
    public:
      static void registerPython();

      TersoffPairTerm(): A(0.0), lambda1(0.0), R(0.0), D(0.0) {
        setShift(0.0);
        setCutoff(infinity);
        preset();
      }

      TersoffPairTerm(real _A, real _lambda1, real _R, real _D, real _cutoff):
                A(_A), lambda1(_lambda1), R(_R), D(_D) {
        setShift(0.0);
        setCutoff(_cutoff);
        preset();
      }

      virtual ~TersoffPairTerm() {};

      void setA(real _A) { A = _A; }
      real getA() const { return A; }

      void setLambda1(real _lambda1) { lambda1 = _lambda1; }
      real getLambda1() const { return lambda1; }
      
      void setR(real _R) { R = _R; }
      real getR() const { return R; }
      
      void setD(real _D) {
        D = _D;
        preset();
      }
      real getD() const { return D; }
      
      // Setter and getter
      
      void preset() {
        Pi_2D = 0.5*(M_PIl/D);
      }

      real _computeEnergySqrRaw(real distSqr) const {

        real d12 = sqrt(distSqr);

        if (d12 > R + D)
          return 0.0;
        else {
          real fR = A*exp(-lambda1*d12);
          if (d12 < R - D)
            return fR;
          else {
             real arg_sin = 0.5*Pi_2D*(d12-R);
             real fC = 0.5*(1.0 - sin(arg_sin));
             real energy = fC * fR;
             return energy;
          }
        }
      }

      bool _computeForceRaw(Real3D& force, const Real3D& r12, real distSqr) const {

        real d12 = sqrt(distSqr);
        real ffactor = 0.0;

        if (d12 > R + D) {
//          force = r12 * ffactor;
          force = 0.0;
        }
        else {
          real inv_d12 = 1.0 / d12;
          real fR = A * exp(-lambda1*d12);
          if (d12 < R - D)
            ffactor = lambda1*fR;
          else {
            real arg_sin = 0.5 * Pi_2D * (d12-R);
            real fC = 0.5 * (1.0 - sin(arg_sin));
            real fC_r = -0.5 * Pi_2D * cos(arg_sin);
            ffactor = -fC_r * fR + lambda1 * fR * fC;
          }      
        force = r12 * inv_d12 * ffactor;
        }
        return true;
      }
    };
  }
}

#endif
