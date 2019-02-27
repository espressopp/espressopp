/*
  Copyright (C) 2017,2018
  Max Planck Institute for Polymer Research

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
#ifndef _INTERACTION_SQUAREWELL_HPP
#define _INTERACTION_SQUAREWELL_HPP

#include "Potential.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"
#include "FixedPairListTypesInteractionTemplate.hpp"
#include <cmath>

namespace espressopp {

  namespace interaction {
    /** This class provides methods to compute forces and energies of
         a squarewell potential based on Leitold and Dellago J. Chem. Phys. 141. 134901 (2014)
    */
    class SmoothSquareWell: public PotentialTemplate < SmoothSquareWell > {
    private:
      real lambda;
      real epsilon;
      real halfepsilon;
      real a;
      real sigma;
      // rb=lambda*sigma. Rightside bounary of the squarewell
      real rb;
    public:
      static void registerPython();
      SmoothSquareWell(): epsilon(0.0), sigma(0.0), a(0.0), lambda(0.0), halfepsilon(0.0) {
        setShift(0.0);
        setCutoff(infinity);
      }

      SmoothSquareWell(real _epsilon, real _sigma, real _cutoff, real _shift): epsilon(_epsilon), sigma(_sigma) {
        halfepsilon = _epsilon/2;
        setShift(_shift);
        setCutoff(_cutoff);
      }

      SmoothSquareWell(real _epsilon, real _sigma, real _cutoff): epsilon(_epsilon), sigma(_sigma) {
        halfepsilon = _epsilon/2;
        autoShift = false;
        setCutoff(_cutoff);
        setAutoShift();
      }

      virtual ~SmoothSquareWell() {};

      void setEpsilon(real _epsilon) {
        epsilon = _epsilon;
        halfepsilon = _epsilon/2;
        updateAutoShift();
      }

      real getEpsilon() const {return epsilon;}

      void setLambda(real _lambda) {
        lambda = _lambda;
        rb = lambda * sigma;
        updateAutoShift();
      }

      real getLambda() const {return lambda;}

      void setSigma(real _sigma) {
        sigma = _sigma;
        updateAutoShift();
      }

      real getSigma() const {return sigma;}

      void setA(real _a) {
        a = _a * sigma;
        updateAutoShift();
      }

      real getA() const {return a;}

      real _computeEnergySqrRaw(real distSqr) const {
        real r = sqrt(distSqr);
        real energy = halfepsilon*(exp((sigma-r)/a)+tanh((r-rb)/a)-1);
        return energy;
      }

      bool _computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) const {
        real r = sqrt(distSqr);
        real ffactor = -halfepsilon/r/a*(-exp((sigma-r)/a)+pow(cosh((r-rb)/a), -2));
        force = dist * ffactor;
        return true;
      }
    };

    // provide pickle support
    struct SmoothSquareWell_pickle : boost::python::pickle_suite
    {
      static
      boost::python::tuple
      getinitargs(SmoothSquareWell const& pot)
      {
        real eps;
        real sig;
        real rc;
        real sh;
        eps=pot.getEpsilon();
        sig=pot.getSigma();
        rc=pot.getCutoff();
        sh=pot.getShift();
        return boost::python::make_tuple(eps, sig, rc, sh);
      }
    };

  }
}

#endif
