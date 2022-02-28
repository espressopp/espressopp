/*
  Copyright (C) 2019-2021
      Max Planck Institute for Polymer Research & JGU Mainz
  Copyright (C) 2012,2013,2014,2015,2016,2017,2018
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
#ifndef VEC_INTERACTION_VERLETLISTLENNARDJONESCAPPED_HPP
#define VEC_INTERACTION_VERLETLISTLENNARDJONESCAPPED_HPP

#include "types.hpp"
#include "interaction/Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "storage/Storage.hpp"

#include "LennardJonesCapped.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "vec/VerletList.hpp"
#include "vec/Vectorization.hpp"

namespace espressopp
{
namespace vec
{
namespace interaction
{
typedef VerletListInteractionTemplate<LennardJonesCapped> VerletListLennardJonesCappedBase;

class VerletListLennardJonesCapped : public VerletListLennardJonesCappedBase
{
protected:
    typedef LennardJonesCapped _Potential;
    typedef _Potential Potential;
    typedef VerletListLennardJonesCappedBase base;

    struct LJCoefficients
    {
        LJCoefficients(Potential const& p)
            : ff1(p.getff1()),
              ff2(p.getff2()),
              caprad(p.getCaprad()),
              capradSqr(p.getCaprad() * p.getCaprad()),
              sigma(p.getSigma()),
              epsilon(p.getEpsilon()),
              cfrac2((sigma / caprad) * (sigma / caprad)),
              cfrac6(cfrac2 * cfrac2 * cfrac2)
        {
        }

        const real ff1;
        const real ff2;
        const real caprad;
        const real capradSqr;
        const real sigma;
        const real epsilon;
        const real cfrac2;
        const real cfrac6;
    };

public:
    VerletListLennardJonesCapped(std::shared_ptr<VerletList> _verletList)
        : base(_verletList), np_types(0), p_types(0)
    {
    }

    void rebuildPotential()
    {
        np_types = potentialArray.size_n();
        p_types = potentialArray.size_m();
        ffs.clear();
        ffs.reserve(np_types * p_types);
        cutoffSqr.clear();
        cutoffSqr.reserve(np_types * p_types);
        for (auto& p : potentialArray)
        {
            ffs.push_back(LJCoefficients(p));
            cutoffSqr.push_back(p.getCutoffSqr());
        }
        needRebuildPotential = false;
    }
    virtual void addForces();

protected:
    template <bool ONETYPE>
    void addForces_impl(ParticleArray& particleArray, VerletList::NeighborList const& neighborList);

    size_t np_types, p_types;
    AlignedVector<LJCoefficients> ffs;
    AlignedVector<real> cutoffSqr;
    bool needRebuildPotential = true;
};

//////////////////////////////////////////////////
// INLINE IMPLEMENTATION
//////////////////////////////////////////////////
inline void VerletListLennardJonesCapped::addForces()
{
    LOG4ESPP_DEBUG(_Potential::theLogger, "loop over verlet list pairs and add forces");

    // lookup table for LJ variables
    // ideal for low number of types vs number of particle pairs
    // trigger rebuild on setPotential and on size modification from getPotential
    auto& pa = verletList->getVectorization()->particles;
    const auto& nl = verletList->getNeighborList();
    const auto vlmaxtype = nl.max_type;

    Potential max_pot = getPotential(vlmaxtype, vlmaxtype);
    if (needRebuildPotential) rebuildPotential();

    if (np_types == 1 && p_types == 1)
    {
        addForces_impl<true>(pa, nl);
    }
    else
    {
        addForces_impl<false>(pa, nl);
    }
}

template <bool ONETYPE>
inline void VerletListLennardJonesCapped::addForces_impl(
    ParticleArray& particleArray, VerletList::NeighborList const& neighborList)
{
    using namespace vec::storage;

    {
        real ff1_, ff2_, caprad_, capradSqr_, epsilon_, cutoffSqr_, cfrac6_;
        if (ONETYPE)
        {
            ff1_ = ffs[0].ff1;
            ff2_ = ffs[0].ff2;
            caprad_ = ffs[0].caprad;
            capradSqr_ = ffs[0].capradSqr;
            epsilon_ = ffs[0].epsilon;
            cfrac6_ = ffs[0].cfrac6;
            cutoffSqr_ = cutoffSqr[0];
        }

        const size_t* __restrict pa_type = particleArray.type.data();
        const real* __restrict pa_p_x = particleArray.p_x.data();
        const real* __restrict pa_p_y = particleArray.p_y.data();
        const real* __restrict pa_p_z = particleArray.p_z.data();
        real* __restrict pa_f_x = particleArray.f_x.data();
        real* __restrict pa_f_y = particleArray.f_y.data();
        real* __restrict pa_f_z = particleArray.f_z.data();

        const auto* __restrict plist = neighborList.plist.data();
        const auto* __restrict prange = neighborList.prange.data();
        const auto* __restrict nplist = neighborList.nplist.data();
        const int ip_max = neighborList.plist.size();

        for (int ip = 0; ip < ip_max; ip++)
        {
            int p = plist[ip];
            int p_lookup;
            if (!ONETYPE)
            {
                p_lookup = pa_type[p] * np_types;
            }
            const real p_x = pa_p_x[p];
            const real p_y = pa_p_y[p];
            const real p_z = pa_p_z[p];

            real f_x = 0.0;
            real f_y = 0.0;
            real f_z = 0.0;

            const int in_min = prange[ip].first;
            const int in_max = prange[ip].second;

            ESPP_VEC_PRAGMAS
            for (int in = in_min; in < in_max; in++)
            {
                auto np_ii = nplist[in];
                {
                    int np_lookup;
                    real dist_x, dist_y, dist_z;

                    {
                        dist_x = p_x - pa_p_x[np_ii];
                        dist_y = p_y - pa_p_y[np_ii];
                        dist_z = p_z - pa_p_z[np_ii];
                        if (!ONETYPE) np_lookup = pa_type[np_ii] + p_lookup;
                    }

                    real distSqr = dist_x * dist_x + dist_y * dist_y + dist_z * dist_z;
                    if (!ONETYPE)
                    {
                        cutoffSqr_ = cutoffSqr[np_lookup];
                    }

                    if (distSqr <= cutoffSqr_)
                    {
                        bool NO_CAP;
                        if (ONETYPE)
                        {
                            NO_CAP = distSqr > capradSqr_;
                        }
                        else
                        {
                            NO_CAP = distSqr > ffs[np_lookup].capradSqr;
                        }

                        real ffactor;
                        if (NO_CAP)
                        {
                            const real frac2 = 1.0 / distSqr;
                            const real frac6 = frac2 * frac2 * frac2;
                            if (ONETYPE)
                                ffactor = ff1_ * frac6 - ff2_;
                            else
                                ffactor = ffs[np_lookup].ff1 * frac6 - ffs[np_lookup].ff2;
                            ffactor = frac6 * ffactor * frac2;
                        }
                        else
                        {
                            if (ONETYPE)
                                ffactor = 48.0 * epsilon_ * cfrac6_ * (cfrac6_ - 0.5) /
                                          (caprad_ * sqrt(distSqr));
                            else
                                ffactor = 48.0 * ffs[np_lookup].epsilon * ffs[np_lookup].cfrac6 *
                                          (ffs[np_lookup].cfrac6 - 0.5) /
                                          (ffs[np_lookup].caprad * sqrt(distSqr));
                        }

                        f_x += dist_x * ffactor;
                        f_y += dist_y * ffactor;
                        f_z += dist_z * ffactor;

                        {
                            pa_f_x[np_ii] -= dist_x * ffactor;
                            pa_f_y[np_ii] -= dist_y * ffactor;
                            pa_f_z[np_ii] -= dist_z * ffactor;
                        }
                    }
                }
            }

            {
                pa_f_x[p] += f_x;
                pa_f_y[p] += f_y;
                pa_f_z[p] += f_z;
            }
        }
    }
}
}  // namespace interaction
}  // namespace vec
}  // namespace espressopp

#endif  // VEC_INTERACTION_VERLETLISTLENNARDJONESCAPPED_HPP
