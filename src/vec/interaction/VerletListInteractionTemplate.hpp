/*
  Copyright (C) 2021
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
#ifndef VEC_INTERACTION_VERLETLISTINTERACTIONTEMPLATE_HPP
#define VEC_INTERACTION_VERLETLISTINTERACTIONTEMPLATE_HPP

#include "types.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "storage/Storage.hpp"
#include "interaction/Interaction.hpp"
#include "vec/VerletList.hpp"
#include "vec/Vectorization.hpp"

namespace espressopp
{
namespace vec
{
namespace interaction
{
using espressopp::interaction::Interaction;
using espressopp::interaction::Nonbonded;
using espressopp::vec::VerletList;

template <typename _Potential>
class VerletListInteractionTemplate : public Interaction
{
protected:
    typedef _Potential Potential;

public:
    VerletListInteractionTemplate(std::shared_ptr<VerletList> _verletList)
        : verletList(_verletList), ntypes(0), potentialArray(0, 0, Potential())
    {
    }

    virtual ~VerletListInteractionTemplate(){};

    void setVerletList(std::shared_ptr<VerletList> _verletList) { verletList = _verletList; }

    std::shared_ptr<VerletList> getVerletList() { return verletList; }

    void setPotential(int type1, int type2, const Potential& potential)
    {
        // typeX+1 because i<ntypes
        ntypes = std::max(ntypes, std::max(type1 + 1, type2 + 1));
        potentialArray.at(type1, type2) = potential;
        LOG4ESPP_INFO(_Potential::theLogger,
                      "added potential for type1=" << type1 << " type2=" << type2);
        if (type1 != type2)
        {  // add potential in the other direction
            potentialArray.at(type2, type1) = potential;
            LOG4ESPP_INFO(_Potential::theLogger, "automatically added the same potential for type1="
                                                     << type2 << " type2=" << type1);
        }
    }

    // this is used in the innermost force-loop
    Potential& getPotential(int type1, int type2) { return potentialArray.at(type1, type2); }

    // this is mainly used to access the potential from Python (e.g. to change parameters of the
    // potential)
    std::shared_ptr<Potential> getPotentialPtr(int type1, int type2)
    {
        return make_shared<Potential>(potentialArray.at(type1, type2));
    }

    virtual void addForces();
    virtual real computeEnergy();
    virtual real computeEnergyDeriv();
    virtual real computeEnergyAA();
    virtual real computeEnergyCG();
    virtual real computeEnergyAA(int atomtype);
    virtual real computeEnergyCG(int atomtype);
    virtual void computeVirialX(std::vector<real>& p_xx_total, int bins);
    virtual real computeVirial();
    virtual void computeVirialTensor(Tensor& w);
    virtual void computeVirialTensor(Tensor& w, real z);
    virtual void computeVirialTensor(Tensor* w, int n);
    virtual real getMaxCutoff();
    virtual int bondType() { return Nonbonded; }

protected:
    std::shared_ptr<VerletList> verletList;
    int ntypes;
    esutil::Array2D<Potential, esutil::enlarge> potentialArray;
};

//////////////////////////////////////////////////
// INLINE IMPLEMENTATION
//////////////////////////////////////////////////
template <typename _Potential>
inline void VerletListInteractionTemplate<_Potential>::addForces()
{
    auto& particles = verletList->getVectorization()->particles;
    const auto& neighborList = verletList->getNeighborList();
    const auto* __restrict plist = neighborList.plist.data();
    const auto* __restrict prange = neighborList.prange.data();
    const auto* __restrict nplist = neighborList.nplist.data();

    const auto vlmaxtype = neighborList.max_type;
    Potential max_pot = potentialArray.at(vlmaxtype, vlmaxtype);  // force a resize

    const int ip_max = neighborList.plist.size();
    for (int ip = 0; ip < ip_max; ip++)
    {
        const int p1 = plist[ip];
        const auto type1 = particles.getType(p1);
        const auto pos1 = particles.getPosition(p1);

        const int in_min = prange[ip].first;
        const int in_max = prange[ip].second;
        for (int in = in_min; in < in_max; in++)
        {
            const int p2 = nplist[in];
            const auto type2 = particles.getType(p2);
            const auto pos2 = particles.getPosition(p2);

            const Potential& potential = potentialArray(type1, type2);

            Real3D force(0.0);
            const Real3D r21 = pos1 - pos2;
            if (potential._computeForce(force, r21))
            {
                particles.addForce(p1, force);
                particles.subForce(p2, force);
            }
        }
    }
}

template <typename _Potential>
inline real VerletListInteractionTemplate<_Potential>::computeEnergy()
{
    real e = 0.0;
    real es = 0.0;

    const auto& particles = verletList->getVectorization()->particles;
    const auto& neighborList = verletList->getNeighborList();
    const auto* __restrict plist = neighborList.plist.data();
    const auto* __restrict prange = neighborList.prange.data();
    const auto* __restrict nplist = neighborList.nplist.data();

    const int ip_max = neighborList.plist.size();
    for (int ip = 0; ip < ip_max; ip++)
    {
        const int p1 = plist[ip];
        const auto type1 = particles.getType(p1);
        const auto pos1 = particles.getPosition(p1);

        const int in_min = prange[ip].first;
        const int in_max = prange[ip].second;
        for (int in = in_min; in < in_max; in++)
        {
            const int p2 = nplist[in];
            const auto type2 = particles.getType(p2);
            const auto pos2 = particles.getPosition(p2);

            const Potential& potential = getPotential(type1, type2);
            e = potential._computeEnergy(pos1 - pos2);
            es += e;
        }
    }

    // reduce over all CPUs
    real esum;
    boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, es, esum, std::plus<real>());
    return esum;
    return 0.0;
}

template <typename _Potential>
inline real VerletListInteractionTemplate<_Potential>::computeEnergyDeriv()
{
    LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeEnergyDeriv() is not yet implemented.");
    return 0.0;
}

template <typename _Potential>
inline real VerletListInteractionTemplate<_Potential>::computeEnergyAA()
{
    LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeEnergyAA() is not yet implemented.");
    return 0.0;
}

template <typename _Potential>
inline real VerletListInteractionTemplate<_Potential>::computeEnergyAA(int atomtype)
{
    LOG4ESPP_WARN(_Potential::theLogger,
                  "Warning! computeEnergyAA(int atomtype) is not yet implemented.");
    return 0.0;
}

template <typename _Potential>
inline real VerletListInteractionTemplate<_Potential>::computeEnergyCG()
{
    LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeEnergyCG() is not yet implemented.");
    return 0.0;
}

template <typename _Potential>
inline real VerletListInteractionTemplate<_Potential>::computeEnergyCG(int atomtype)
{
    LOG4ESPP_WARN(_Potential::theLogger,
                  "Warning! computeEnergyCG(int atomtype) is not yet implemented.");
    return 0.0;
}

template <typename _Potential>
inline void VerletListInteractionTemplate<_Potential>::computeVirialX(std::vector<real>& p_xx_total,
                                                                      int bins)
{
    LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeVirialX() is not yet implemented.");
}

template <typename _Potential>
inline real VerletListInteractionTemplate<_Potential>::computeVirial()
{
    real w = 0.0;

    const auto& particles = verletList->getVectorization()->particles;
    const auto& neighborList = verletList->getNeighborList();
    const auto* __restrict plist = neighborList.plist.data();
    const auto* __restrict prange = neighborList.prange.data();
    const auto* __restrict nplist = neighborList.nplist.data();

    const int ip_max = neighborList.plist.size();
    for (int ip = 0; ip < ip_max; ip++)
    {
        const int p1 = plist[ip];
        const auto type1 = particles.getType(p1);
        const auto pos1 = particles.getPosition(p1);

        const int in_min = prange[ip].first;
        const int in_max = prange[ip].second;
        for (int in = in_min; in < in_max; in++)
        {
            const int p2 = nplist[in];
            const auto type2 = particles.getType(p2);
            const auto pos2 = particles.getPosition(p2);

            const Potential& potential = getPotential(type1, type2);

            Real3D force(0.0, 0.0, 0.0);
            const Real3D r21 = pos1 - pos2;
            if (potential._computeForce(force, r21))
            {
                w = w + r21 * force;
            }
        }
    }

    // reduce over all CPUs
    real wsum;
    boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
    return wsum;
}

template <typename _Potential>
inline void VerletListInteractionTemplate<_Potential>::computeVirialTensor(Tensor& w)
{
    Tensor wlocal(0.0);

    const auto& particles = verletList->getVectorization()->particles;
    const auto& neighborList = verletList->getNeighborList();
    const auto* __restrict plist = neighborList.plist.data();
    const auto* __restrict prange = neighborList.prange.data();
    const auto* __restrict nplist = neighborList.nplist.data();

    const int ip_max = neighborList.plist.size();
    for (int ip = 0; ip < ip_max; ip++)
    {
        const int p1 = plist[ip];
        const auto type1 = particles.getType(p1);
        const auto pos1 = particles.getPosition(p1);

        const int in_min = prange[ip].first;
        const int in_max = prange[ip].second;
        for (int in = in_min; in < in_max; in++)
        {
            const int p2 = nplist[in];
            const auto type2 = particles.getType(p2);
            const auto pos2 = particles.getPosition(p2);

            const Potential& potential = getPotential(type1, type2);

            Real3D force(0.0, 0.0, 0.0);
            const Real3D r21 = pos1 - pos2;
            if (potential._computeForce(force, r21))
            {
                wlocal += Tensor(r21, force);
            }
        }
    }

    // reduce over all CPUs
    Tensor wsum(0.0);
    boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
    w += wsum;
}

// local pressure tensor for layer, plane is defined by z coordinate
template <typename _Potential>
inline void VerletListInteractionTemplate<_Potential>::computeVirialTensor(Tensor& w, real z)
{
    System& system = verletList->getSystemRef();
    Real3D Li = system.bc->getBoxL();

    real rc_cutoff = verletList->getVerletCutoff();

    // boundaries should be taken into account
    bool ghost_layer = false;
    real zghost = -100.0;
    if (z < rc_cutoff)
    {
        zghost = z + Li[2];
        ghost_layer = true;
    }
    else if (z >= Li[2] - rc_cutoff)
    {
        zghost = z - Li[2];
        ghost_layer = true;
    }

    Tensor wlocal(0.0);

    const auto& particles = verletList->getVectorization()->particles;
    const auto& neighborList = verletList->getNeighborList();
    const auto* __restrict plist = neighborList.plist.data();
    const auto* __restrict prange = neighborList.prange.data();
    const auto* __restrict nplist = neighborList.nplist.data();

    const int ip_max = neighborList.plist.size();
    for (int ip = 0; ip < ip_max; ip++)
    {
        const int p1 = plist[ip];
        const auto type1 = particles.getType(p1);
        const auto pos1 = particles.getPosition(p1);

        const int in_min = prange[ip].first;
        const int in_max = prange[ip].second;
        for (int in = in_min; in < in_max; in++)
        {
            const int p2 = nplist[in];
            const auto type2 = particles.getType(p2);
            const auto pos2 = particles.getPosition(p2);

            if ((pos1[2] > z && pos2[2] < z) || (pos1[2] < z && pos2[2] > z) ||
                (ghost_layer && ((pos1[2] > zghost && pos2[2] < zghost) ||
                                 (pos1[2] < zghost && pos2[2] > zghost))))
            {
                const Potential& potential = getPotential(type1, type2);

                Real3D force(0.0, 0.0, 0.0);
                Real3D r21 = pos1 - pos2;
                if (potential._computeForce(force, r21))
                {
                    wlocal += Tensor(r21, force) / fabs(r21[2]);
                }
            }
        }
    }
    // reduce over all CPUs
    Tensor wsum(0.0);
    boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
    w += wsum;
}

// it will calculate the pressure in 'n' layers along Z axis
// the first layer has coordinate 0.0 the last - (Lz - Lz/n)
template <typename _Potential>
inline void VerletListInteractionTemplate<_Potential>::computeVirialTensor(Tensor* w, int n)
{
    LOG4ESPP_DEBUG(
        _Potential::theLogger,
        "loop over verlet list pairs and sum up virial tensor in bins along z-direction");

    System& system = verletList->getSystemRef();
    Real3D Li = system.bc->getBoxL();

    real z_dist = Li[2] / float(n);  // distance between two layers
    Tensor* wlocal = new Tensor[n];
    for (int i = 0; i < n; i++) wlocal[i] = Tensor(0.0);

    const auto& particles = verletList->getVectorization()->particles;
    const auto& neighborList = verletList->getNeighborList();
    const auto* __restrict plist = neighborList.plist.data();
    const auto* __restrict prange = neighborList.prange.data();
    const auto* __restrict nplist = neighborList.nplist.data();

    const int ip_max = neighborList.plist.size();
    for (int ip = 0; ip < ip_max; ip++)
    {
        const int p1 = plist[ip];
        const auto type1 = particles.getType(p1);
        const auto pos1 = particles.getPosition(p1);

        const int in_min = prange[ip].first;
        const int in_max = prange[ip].second;
        for (int in = in_min; in < in_max; in++)
        {
            const int p2 = nplist[in];
            const auto type2 = particles.getType(p2);
            const auto pos2 = particles.getPosition(p2);

            const Potential& potential = getPotential(type1, type2);

            Real3D force(0.0, 0.0, 0.0);
            Real3D r21 = pos1 - pos2;
            Tensor ww;
            if (potential._computeForce(force, r21))
            {
                ww = Tensor(r21, force) / fabs(r21[2]);

                int position1 = (int)(pos1[2] / z_dist);
                int position2 = (int)(pos2[2] / z_dist);

                int maxpos = std::max(position1, position2);
                int minpos = std::min(position1, position2);

                // boundaries should be taken into account
                bool boundaries1 = false;
                bool boundaries2 = false;
                if (minpos < 0)
                {
                    minpos += n;
                    boundaries1 = true;
                }
                if (maxpos >= n)
                {
                    maxpos -= n;
                    boundaries2 = true;
                }

                if (boundaries1 || boundaries2)
                {
                    for (int i = 0; i <= maxpos; i++)
                    {
                        wlocal[i] += ww;
                    }
                    for (int i = minpos + 1; i < n; i++)
                    {
                        wlocal[i] += ww;
                    }
                }
                else
                {
                    for (int i = minpos + 1; i <= maxpos; i++)
                    {
                        wlocal[i] += ww;
                    }
                }
            }
        }
    }

    // reduce over all CPUs
    Tensor* wsum = new Tensor[n];
    boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, n, (double*)&wsum, std::plus<double>());

    for (int j = 0; j < n; j++)
    {
        w[j] += wsum[j];
    }

    delete[] wsum;
    delete[] wlocal;

    LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeVirialTensor() is not yet implemented.");
}

template <typename _Potential>
inline real VerletListInteractionTemplate<_Potential>::getMaxCutoff()
{
    real cutoff = 0.0;
    for (int i = 0; i < ntypes; i++)
    {
        for (int j = 0; j < ntypes; j++)
        {
            cutoff = std::max(cutoff, getPotential(i, j).getCutoff());
        }
    }
    return cutoff;
}
}  // namespace interaction
}  // namespace vec
}  // namespace espressopp

#endif  // VEC_INTERACTION_VERLETLISTINTERACTIONTEMPLATE_HPP
