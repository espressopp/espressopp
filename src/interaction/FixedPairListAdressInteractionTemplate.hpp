/*
  Copyright (c) 2015,
      Jakub Krajniak (jkrajniak (at) gmail.com)

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
#ifndef _INTERACTION_FIXEDPAIRLISTADRESSINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDPAIRLISTADRESSINTERACTIONTEMPLATE_HPP

#include <algorithm>
#include <functional>
#include <vector>
#include "mpi.hpp"
#include "integrator/Adress.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "FixedPairList.hpp"
#include "FixedPairListAdress.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "Interaction.hpp"
#include "types.hpp"

namespace espressopp
{
namespace interaction
{
template <typename _Potential>
class FixedPairListAdressInteractionTemplate : public Interaction, SystemAccess
{
protected:
    typedef _Potential Potential;

public:
    FixedPairListAdressInteractionTemplate(std::shared_ptr<System> system,
                                           std::shared_ptr<FixedPairList> _fixedpairList,
                                           std::shared_ptr<Potential> _potential,
                                           bool _cgpotential)
        : SystemAccess(system),
          fixedpairList(_fixedpairList),
          potential(_potential),
          cgPotential(_cgpotential)
    {
        if (!potential)
        {
            LOG4ESPP_ERROR(theLogger, "NULL potential");
        }
        scaleFactor_ = 1.0;
    }

    virtual ~FixedPairListAdressInteractionTemplate() {}

    void setFixedPairList(std::shared_ptr<FixedPairList> _fixedpairList)
    {
        fixedpairList = _fixedpairList;
    }

    std::shared_ptr<FixedPairList> getFixedPairList() { return fixedpairList; }

    void setPotential(std::shared_ptr<Potential> _potential)
    {
        if (_potential)
        {
            potential = _potential;
        }
        else
        {
            LOG4ESPP_ERROR(theLogger, "NULL potential");
        }
    }

    std::shared_ptr<Potential> getPotential() { return potential; }

    void setScaleFactor(real s) { scaleFactor_ = s; }
    real scaleFactor() { return scaleFactor_; }

    virtual void addForces();
    virtual real computeEnergy();
    virtual real computeEnergyDeriv();
    virtual real computeEnergyAA();
    virtual real computeEnergyCG();
    virtual real computeEnergyAA(int atomtype);
    virtual real computeEnergyCG(int atomtype);
    virtual void computeVirialX(std::vector<real> &p_xx_total, int bins);
    virtual real computeVirial();
    virtual void computeVirialTensor(Tensor &w);
    virtual void computeVirialTensor(Tensor &w, real z);
    virtual void computeVirialTensor(Tensor *w, int n);
    virtual real getMaxCutoff();
    virtual int bondType() { return Pair; }

protected:
    int ntypes;
    std::shared_ptr<FixedPairList> fixedpairList;
    std::shared_ptr<Potential> potential;
    bool cgPotential;
    real scaleFactor_;
};

//////////////////////////////////////////////////
// INLINE IMPLEMENTATION
//////////////////////////////////////////////////
template <typename _Potential>
inline void FixedPairListAdressInteractionTemplate<_Potential>::addForces()
{
    LOG4ESPP_INFO(_Potential::theLogger, "Adding forces of FixedPairListAdress");
    const bc::BC &bc = *getSystemRef().bc;  // boundary conditions
    real ltMaxBondSqr = fixedpairList->getLongtimeMaxBondSqr();
    for (FixedPairList::PairList::Iterator it(*fixedpairList); it.isValid(); ++it)
    {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;

        Real3D dist;
        bc.getMinimumImageVectorBox(dist, p1.position(), p2.position());
        Real3D force;
        real d = dist.sqr();
        if (d > ltMaxBondSqr)
        {
            fixedpairList->setLongtimeMaxBondSqr(d);
            ltMaxBondSqr = d;
        }

        // Scaling for the Adress extension. If this potential works between CG beads
        // then the Force is scalded by (1-w12) otherwise it works across the beads,
        // with scaling w12.
        real p1lambda = p1.lambda();
        real p2lambda = p2.lambda();
        if (!cgPotential && (p1lambda < 0.0 || p2lambda < 0.0)) continue;

        real w12 = p1lambda * p2lambda;
        if (w12 < 0.0) w12 = 0.0;
        real forcescale12 = w12;
        if (cgPotential)
        {
            forcescale12 = (1 - w12);
        }

        forcescale12 *= scaleFactor_;

        if (!isAlmostZero(forcescale12))
        {
            if (potential->_computeForce(force, dist))
            {
                p1.force() += forcescale12 * force;
                p2.force() -= forcescale12 * force;
            }
        }
    }
}

template <typename _Potential>
inline real FixedPairListAdressInteractionTemplate<_Potential>::computeEnergy()
{
    LOG4ESPP_INFO(theLogger, "compute energy of the FixedPairListAdress pairs");

    real e_local = 0.0;
    const bc::BC &bc = *getSystemRef().bc;  // boundary conditions
    for (FixedPairList::PairList::Iterator it(*fixedpairList); it.isValid(); ++it)
    {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;

        real p1lambda = p1.lambda();
        real p2lambda = p2.lambda();
        if (!cgPotential && (p1lambda < 0.0 || p2lambda < 0.0)) continue;

        real w12 = p1lambda * p2lambda;
        if (w12 < 0.0) w12 = 0.0;
        real energyscale12 = w12;
        if (cgPotential) energyscale12 = 1.0 - w12;

        energyscale12 *= scaleFactor_;

        if (!isAlmostZero(energyscale12))
        {
            Real3D r21;
            bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
            e_local += energyscale12 * potential->_computeEnergy(r21);
        }
    }
    real esum = 0.0;
    boost::mpi::all_reduce(*mpiWorld, e_local, esum, std::plus<real>());
    return esum;
}

template <typename _Potential>
inline real FixedPairListAdressInteractionTemplate<_Potential>::computeEnergyDeriv()
{
    std::cout << "Warning! At the moment computeEnergyDeriv() in "
                 "FixedPairListAdressInteractionTemplate does not work."
              << std::endl;
    return 0.0;
}

template <typename _Potential>
inline real FixedPairListAdressInteractionTemplate<_Potential>::computeEnergyAA()
{
    if (!cgPotential) return computeEnergy();
    return 0.0;
}

template <typename _Potential>
inline real FixedPairListAdressInteractionTemplate<_Potential>::computeEnergyCG()
{
    if (cgPotential) return computeEnergy();
    return 0.0;
}

template <typename _Potential>
inline real FixedPairListAdressInteractionTemplate<_Potential>::computeEnergyAA(int atomtype)
{
    std::cout << "Warning! At the moment computeEnergyAA() in "
                 "FixedPairListAdressInteractionTemplate does not work."
              << std::endl;
    return 0.0;
}

template <typename _Potential>
inline real FixedPairListAdressInteractionTemplate<_Potential>::computeEnergyCG(int atomtype)
{
    std::cout << "Warning! At the moment computeEnergyCG() in "
                 "FixedPairListAdressInteractionTemplate does not work."
              << std::endl;
    return 0.0;
}

template <typename _Potential>
inline void FixedPairListAdressInteractionTemplate<_Potential>::computeVirialX(
    std::vector<real> &p_xx_total, int bins)
{
    std::cout << "Warning! At the moment computeVirialX() in "
              << "FixedPairListAdressInteractionTemplate does not work." << std::endl;
}

template <typename _Potential>
inline real FixedPairListAdressInteractionTemplate<_Potential>::computeVirial()
{
    LOG4ESPP_INFO(theLogger, "compute the virial for the FixedPair List");

    real w_virial = 0.0;
    const bc::BC &bc = *getSystemRef().bc;  // boundary conditions
    for (FixedPairList::PairList::Iterator it(*fixedpairList); it.isValid(); ++it)
    {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;

        real p1lambda = p1.lambda();
        real p2lambda = p2.lambda();
        if (!cgPotential && (p1lambda < 0.0 || p2lambda < 0.0)) continue;

        real w12 = p1lambda * p2lambda;
        if (w12 < 0.0) w12 = 0.0;
        real forcescale = w12;
        if (cgPotential)
        {
            forcescale = (1.0 - w12);
        }

        forcescale *= scaleFactor_;

        if (!isAlmostZero(forcescale))
        {
            Real3D r21;
            bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
            Real3D force;
            if (potential->_computeForce(force, r21))
            {
                w_virial += r21 * forcescale * force;
            }
        }
    }

    real wsum;
    boost::mpi::all_reduce(*mpiWorld, w_virial, wsum, std::plus<real>());
    return wsum;
}

template <typename _Potential>
inline void FixedPairListAdressInteractionTemplate<_Potential>::computeVirialTensor(Tensor &w)
{
    LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");

    Tensor w_wlocal = 0.0;
    const bc::BC &bc = *getSystemRef().bc;  // boundary conditions
    for (FixedPairList::PairList::Iterator it(*fixedpairList); it.isValid(); ++it)
    {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;

        real p1lambda = p1.lambda();
        real p2lambda = p2.lambda();
        if (!cgPotential && (p1lambda < 0.0 || p2lambda < 0.0)) continue;

        real w12 = p1lambda * p2lambda;
        if (w12 < 0.0) w12 = 0.0;
        real forcescale = w12;
        if (cgPotential)
        {
            forcescale = (1.0 - w12);
        }

        forcescale *= scaleFactor_;

        if (!isAlmostZero(forcescale))
        {
            Real3D r21;
            bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
            Real3D force;
            if (potential->_computeForce(force, r21))
            {
                w_wlocal += Tensor(r21, forcescale * force);
            }
        }
    }

    // reduce over all CPUs
    Tensor wsum(0.0);
    boost::mpi::all_reduce(*mpiWorld, reinterpret_cast<double *>(&w_wlocal), 6,
                           reinterpret_cast<double *>(&wsum), std::plus<double>());
    w += wsum;
}

template <typename _Potential>
inline void FixedPairListAdressInteractionTemplate<_Potential>::computeVirialTensor(Tensor &w,
                                                                                    real z)
{
    LOG4ESPP_ERROR(theLogger,
                   "compute the virial tensor for the FixedPair List Adress not implemented");
}

template <typename _Potential>
inline void FixedPairListAdressInteractionTemplate<_Potential>::computeVirialTensor(Tensor *w,
                                                                                    int n)
{
    LOG4ESPP_INFO(theLogger,
                  "compute the virial tensor for the FixedPair List Adress not implemented");
}

template <typename _Potential>
inline real FixedPairListAdressInteractionTemplate<_Potential>::getMaxCutoff()
{
    return potential->getCutoff();
}

}  // end namespace interaction
}  // end namespace espressopp
#endif
