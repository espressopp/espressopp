/*
  Copyright (C) 2021
      Max Planck Institute for Polymer Research & JGU Mainz
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
#ifndef VEC_INTERACTION_LENNARDJONESCAPPED_HPP
#define VEC_INTERACTION_LENNARDJONESCAPPED_HPP

#include "interaction/Potential.hpp"

namespace espressopp
{
namespace vec
{
namespace interaction
{
using espressopp::interaction::PotentialTemplate;

/** This class provides methods to compute forces and energies of
    the Lennard Jones potential with capped forces.

    \f[ V(r) = 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
    \left( \frac{\sigma}{r} \right)^{6} \right]
    \f]

*/
class LennardJonesCapped : public PotentialTemplate<LennardJonesCapped>
{
private:
    real epsilon;
    real sigma;
    real caprad;
    real ff1, ff2;
    real ef1, ef2;

public:
    inline const real& getCutoffSqr() const { return cutoffSqr; }
    inline const real& getff1() const { return ff1; }
    inline const real& getff2() const { return ff2; }

    static void registerPython();

    LennardJonesCapped() : epsilon(0.0), sigma(0.0)
    {
        setShift(0.0);
        setCutoff(infinity);
        preset();
    }

    LennardJonesCapped(real _epsilon, real _sigma, real _cutoff, real _caprad, real _shift)
        : epsilon(_epsilon), sigma(_sigma), caprad(_caprad)
    {
        setShift(_shift);
        setCutoff(_cutoff);
        preset();
    }

    LennardJonesCapped(real _epsilon, real _sigma, real _cutoff, real _caprad)
        : epsilon(_epsilon), sigma(_sigma), caprad(_caprad)
    {
        autoShift = false;
        setCutoff(_cutoff);
        preset();
        setAutoShift();
    }

    virtual ~LennardJonesCapped(){};

    void preset()
    {
        real sig2 = sigma * sigma;
        real sig6 = sig2 * sig2 * sig2;
        ff1 = 48.0 * epsilon * sig6 * sig6;
        ff2 = 24.0 * epsilon * sig6;
        ef1 = 4.0 * epsilon * sig6 * sig6;
        ef2 = 4.0 * epsilon * sig6;
    }

    // Setter and getter
    void setEpsilon(real _epsilon)
    {
        epsilon = _epsilon;
        updateAutoShift();
        preset();
    }

    real getEpsilon() const { return epsilon; }

    void setCaprad(real _caprad)
    {
        caprad = _caprad;
        updateAutoShift();
        preset();
    }
    real getCaprad() const { return caprad; }

    void setSigma(real _sigma)
    {
        sigma = _sigma;
        updateAutoShift();
        preset();
    }
    real getSigma() const { return sigma; }

    real _computeEnergySqrRaw(real distSqr) const
    {
        real capradSqr = caprad * caprad;

        if (distSqr > capradSqr)
        {
            real frac2 = sigma * sigma / distSqr;
            real frac6 = frac2 * frac2 * frac2;
            real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
            return energy;
        }
        else
        {  // capped
            real frac2 = sigma * sigma / capradSqr;
            real frac6 = frac2 * frac2 * frac2;
            real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
            return energy;
        }
    }

    bool _computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) const
    {
        real capradSqr = caprad * caprad;

        if (distSqr > capradSqr)
        {
            real frac2 = 1.0 / distSqr;
            real frac6 = frac2 * frac2 * frac2;
            real ffactor = frac6 * (ff1 * frac6 - ff2) * frac2;
            force = dist * ffactor;
            return true;
        }
        else
        {  // capped part
            real frac2 = (sigma / caprad) * (sigma / caprad);
            real frac6 = frac2 * frac2 * frac2;
            real ffactor = 48.0 * epsilon * frac6 * (frac6 - 0.5) / (caprad * sqrt(distSqr));
            force = dist * ffactor;
            return true;

            /*real frac2 = 1.0 / sqrt(distSqr)*caprad;
            real frac6 = frac2 * frac2 * frac2;
            real ffactor = frac6 * (ff1 * frac6 - ff2) * frac2;
            force = dist * ffactor;
            return true;*/
        }
    }
};
// provide pickle support
struct LennardJonesCapped_pickle : boost::python::pickle_suite
{
    static boost::python::tuple getinitargs(LennardJonesCapped const& pot)
    {
        real eps;
        real sig;
        real rc;
        real sh;
        eps = pot.getEpsilon();
        sig = pot.getSigma();
        rc = pot.getCutoff();
        sh = pot.getShift();
        return boost::python::make_tuple(eps, sig, rc, sh);
    }
};

}  // namespace interaction
}  // namespace vec
}  // namespace espressopp

#endif  // VEC_INTERACTION_LENNARDJONESCAPPED_HPP
