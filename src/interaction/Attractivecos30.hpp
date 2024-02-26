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
#ifndef _INTERACTION_ATTRACTIVECOS30_HPP
#define _INTERACTION_ATTRACTIVECOS30_HPP

#include "FixedPairListInteractionTemplate.hpp"
#include "FixedPairListTypesInteractionTemplate.hpp"
#include "Potential.hpp"

namespace espressopp {
  namespace interaction {
/** This class provides methods to compute forces and energies of
        the short-range nonbonded attractive potential between r_cut and r^a_cut.

        \f[ V(r) = \alpha [\cos (\mu (r/r_{cut})^2)]
        \f]

*/

    class Attractivecos30 : public PotentialTemplate< Attractivecos30 > {
    private:
      real alpha;
      real mu;
      real ff, ef;
      real rcut,rc2,racut;

    public:
      static void registerPython();

      Attractivecos30() : alpha(0.0), mu(0.0)
      {
         setShift(0.0);
         setCutoff(infinity);
         preset();
      }

      Attractivecos30(real _alpha, real _mu, real _cutoff, real _shift) : alpha(_alpha), mu(_mu)
      {
         setShift(_shift);
         setCutoff(_cutoff);
         preset();
      }
 
      Attractivecos30(real _alpha, real _mu, real _cutoff) : alpha(_alpha), mu(_mu)
      {
         autoShift = false;
         setCutoff(_cutoff);
         preset();
         setAutoShift();
      }
      virtual ~Attractivecos30() {};

      void preset()
      {
         rcut = pow (2,(1.0/6.0));
         rc2=rcut*rcut;
         ff = 2.0*alpha*mu/rc2;
         ef = alpha;
         racut = cutoff;
      }

    // Setter and getter
       void setAlpha(real _alpha) {
         alpha = _alpha;
         LOG4ESPP_INFO(theLogger, "alpha=" << alpha);
         updateAutoShift();
         preset();
       }

      real getAlpha() const { return alpha; }

      void setMu(real _mu) {
        mu = _mu;
        LOG4ESPP_INFO(theLogger, "mu=" << mu);
        updateAutoShift();
        preset();
      }
      real getMu() const { return mu; }

      real _computeEnergySqrRaw(real distSqr) const
      {
        real energy;
/* shift the energy to zero at the cutoff (auto) */
        if (distSqr>=rc2) energy= ef*(cos(mu*distSqr/rc2));
//          if (mu>0) {std::cout<< "r2=" << distSqr << " mu=" << mu << " rc2=" << rc2 << " ef=" << ef << " energy = " << energy << "\n";}
        else energy=0;

        return energy;

    }

    bool _computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) const
    {

        if (distSqr>=rc2) force = ff * sin(mu*distSqr/rc2) *dist;
//             if (ff>0) {std::cout<< "force=" << force << "\n";}}

        return true;

      }
      static LOG4ESPP_DECL_LOGGER(theLogger);
    };


// provide pickle support
struct Attractivecos30_pickle : boost::python::pickle_suite
{
    static boost::python::tuple getinitargs(Attractivecos30 const& pot)
    {
        real alp;
        real muu;
        real src;

        alp=pot.getAlpha();
        muu=pot.getMu();
        src =pot.getCutoff();
        return boost::python::make_tuple(alp, muu, src);
    }
};

}  // namespace interaction
}  // namespace espressopp

#endif

