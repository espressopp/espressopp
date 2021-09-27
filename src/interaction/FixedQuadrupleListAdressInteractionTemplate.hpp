/*
  Copyright (c) 2015
      Jakub Krajniak (jkrajniak at gmail.com)

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
#ifndef _INTERACTION_FIXEDQUADRUPLELISTADRESSINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDQUADRUPLELISTADRESSINTERACTIONTEMPLATE_HPP

#include <functional>
#include <vector>

#include "mpi.hpp"
#include "integrator/Adress.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "FixedQuadrupleList.hpp"
#include "FixedQuadrupleListAdress.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "types.hpp"

namespace espressopp
{
namespace interaction
{
template <typename _DihedralPotential>
class FixedQuadrupleListAdressInteractionTemplate : public Interaction, SystemAccess
{
protected:
    typedef _DihedralPotential Potential;

public:
    FixedQuadrupleListAdressInteractionTemplate(
        std::shared_ptr<System> _system,
        std::shared_ptr<FixedQuadrupleList> _fixedquadrupleList,
        std::shared_ptr<Potential> _potential,
        bool _cgPotential)
        : SystemAccess(_system),
          fixedquadrupleList(_fixedquadrupleList),
          potential(_potential),
          cgPotential(_cgPotential)
    {
        if (!potential)
        {
            LOG4ESPP_ERROR(theLogger, "FixedQuadrupleListAdressInteraction potential");
        }
    }

    void setFixedQuadrupleList(std::shared_ptr<FixedQuadrupleList> _fixedquadrupleList)
    {
        fixedquadrupleList = _fixedquadrupleList;
    }

    virtual ~FixedQuadrupleListAdressInteractionTemplate() {}

    std::shared_ptr<FixedQuadrupleList> getFixedQuadrupleList() { return fixedquadrupleList; }

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
    virtual int bondType() { return Dihedral; }

protected:
    int ntypes;
    std::shared_ptr<FixedQuadrupleList> fixedquadrupleList;
    std::shared_ptr<Potential> potential;
    bool cgPotential;
};

/** Inline implementation */
template <typename _DihedralPotential>
inline void FixedQuadrupleListAdressInteractionTemplate<_DihedralPotential>::addForces()
{
    LOG4ESPP_INFO(theLogger, "add forces computed by FixedQuadrupleList");

    const bc::BC &bc = *getSystemRef().bc;  // boundary conditions

    for (FixedQuadrupleList::QuadrupleList::Iterator it(*fixedquadrupleList); it.isValid(); ++it)
    {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        Particle &p3 = *it->third;
        Particle &p4 = *it->fourth;

        real p1lambda = p1.lambda();
        real p2lambda = p2.lambda();
        real p3lambda = p3.lambda();
        real p4lambda = p4.lambda();
        if (!cgPotential && (p1lambda < 0.0 || p2lambda < 0.0 || p3lambda < 0.0 || p4lambda < 0.0))
            continue;

        real w1234 = pow(p1lambda * p2lambda * p3lambda * p4lambda, 0.5);
        if (w1234 < 0.0) w1234 = 0.0;
        real forcescale1234 = w1234;
        if (cgPotential)
        {
            forcescale1234 = (1.0 - w1234);
        }

        if (!isAlmostZero(forcescale1234))
        {
            LOG4ESPP_DEBUG(theLogger,
                           "scalling quadruple list force with weight " << forcescale1234);
            Real3D dist21, dist32, dist43;

            bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
            bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
            bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

            Real3D force1, force2, force3, force4;  // result forces

            potential->_computeForce(force1, force2, force3, force4, dist21, dist32, dist43);
            p1.force() += forcescale1234 * force1;
            p2.force() += forcescale1234 * force2;
            p3.force() += forcescale1234 * force3;
            p4.force() += forcescale1234 * force4;
        }
    }
}

template <typename _DihedralPotential>
inline real FixedQuadrupleListAdressInteractionTemplate<_DihedralPotential>::computeEnergy()
{
    LOG4ESPP_INFO(theLogger, "compute energy of the quadruples");

    const bc::BC &bc = *getSystemRef().bc;  // boundary conditions
    real e = 0.0;
    for (FixedQuadrupleList::QuadrupleList::Iterator it(*fixedquadrupleList); it.isValid(); ++it)
    {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        const Particle &p4 = *it->fourth;

        real p1lambda = p1.lambda();
        real p2lambda = p2.lambda();
        real p3lambda = p3.lambda();
        real p4lambda = p4.lambda();
        if (!cgPotential && (p1lambda < 0.0 || p2lambda < 0.0 || p3lambda < 0.0 || p4lambda < 0.0))
            continue;

        real w1234 = pow(p1lambda * p2lambda * p3lambda * p4lambda, 0.5);
        if (w1234 < 0.0) w1234 = 0.0;
        real forcescale1234 = w1234;
        if (cgPotential)
        {
            forcescale1234 = (1.0 - w1234);
        }

        if (!isAlmostZero(forcescale1234))
        {
            Real3D dist21, dist32, dist43;

            bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
            bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
            bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

            e += forcescale1234 * potential->_computeEnergy(dist21, dist32, dist43);
        }
    }
    real esum;
    boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
    return esum;
}

template <typename _Potential>
inline real FixedQuadrupleListAdressInteractionTemplate<_Potential>::computeEnergyDeriv()
{
    std::cout << "Warning! At the moment computeEnergyDeriv() in "
                 "FixedQuadrupleListAdressInteractionTemplate does not work."
              << std::endl;
    return 0.0;
}

template <typename _DihedralPotential>
inline real FixedQuadrupleListAdressInteractionTemplate<_DihedralPotential>::computeEnergyAA()
{
    if (!cgPotential) return computeEnergy();
    return 0.0;
}

template <typename _DihedralPotential>
inline real FixedQuadrupleListAdressInteractionTemplate<_DihedralPotential>::computeEnergyCG()
{
    if (cgPotential) return computeEnergy();
    return 0.0;
}

template <typename _Potential>
inline real FixedQuadrupleListAdressInteractionTemplate<_Potential>::computeEnergyAA(int atomtype)
{
    std::cout << "Warning! At the moment computeEnergyAA() in "
                 "FixedQuadrupleListAdressInteractionTemplate does not work."
              << std::endl;
    return 0.0;
}

template <typename _Potential>
inline real FixedQuadrupleListAdressInteractionTemplate<_Potential>::computeEnergyCG(int atomtype)
{
    std::cout << "Warning! At the moment computeEnergyCG() in "
                 "FixedQuadrupleListAdressInteractionTemplate does not work."
              << std::endl;
    return 0.0;
}

template <typename _DihedralPotential>
inline void FixedQuadrupleListAdressInteractionTemplate<_DihedralPotential>::computeVirialX(
    std::vector<real> &p_xx_total, int bins)
{
    std::cout << "Warning! At the moment computeVirialX in ";
    std::cout << "FixedQuadrupleListAdressInteractionTemplate does not work." << std::endl
              << "Therefore, the corresponding interactions won't be included in calculation."
              << std::endl;
}

template <typename _DihedralPotential>
inline real FixedQuadrupleListAdressInteractionTemplate<_DihedralPotential>::computeVirial()
{
    LOG4ESPP_INFO(theLogger, "compute scalar virial of the quadruples");

    real w = 0.0;
    const bc::BC &bc = *getSystemRef().bc;  // boundary conditions
    for (FixedQuadrupleList::QuadrupleList::Iterator it(*fixedquadrupleList); it.isValid(); ++it)
    {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        const Particle &p4 = *it->fourth;

        real p1lambda = p1.lambda();
        real p2lambda = p2.lambda();
        real p3lambda = p3.lambda();
        real p4lambda = p4.lambda();
        if (!cgPotential && (p1lambda < 0.0 || p2lambda < 0.0 || p3lambda < 0.0 || p4lambda < 0.0))
            continue;

        real w1234 = pow(p1lambda * p2lambda * p3lambda * p4lambda, 0.5);
        if (w1234 < 0.0) w1234 = 0.0;
        real forcescale1234 = w1234;
        if (cgPotential)
        {
            forcescale1234 = (1.0 - w1234);
        }

        if (!isAlmostZero(forcescale1234))
        {
            Real3D dist21, dist32, dist43;

            bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
            bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
            bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

            Real3D force1, force2, force3, force4;

            potential->_computeForce(force1, force2, force3, force4, dist21, dist32, dist43);

            w += dist21 * force1 + dist32 * force2;
        }
    }

    real wsum;
    boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
    return w;
}

template <typename _DihedralPotential>
inline void FixedQuadrupleListAdressInteractionTemplate<_DihedralPotential>::computeVirialTensor(
    Tensor &w)
{
    LOG4ESPP_INFO(theLogger, "compute the virial tensor of the quadruples");

    std::cout << "Warning!!! computeVirialTensor in specified volume doesn't work for "
                 "FixedQuadrupleListAdressInteractionTemplate at the moment"
              << std::endl;
}

template <typename _DihedralPotential>
inline void FixedQuadrupleListAdressInteractionTemplate<_DihedralPotential>::computeVirialTensor(
    Tensor &w, real z)
{
    LOG4ESPP_INFO(theLogger, "compute the virial tensor of the quadruples");

    std::cout << "Warning!!! computeVirialTensor in "
              << "specified volume doesn't work for "
              << "FixedQuadrupleListAdressInteractionTemplate at the moment" << std::endl;
}

template <typename _DihedralPotential>
inline void FixedQuadrupleListAdressInteractionTemplate<_DihedralPotential>::computeVirialTensor(
    Tensor *w, int n)
{
    std::cout << "Warning!!! computeVirialTensor in specified volume doesn't work for "
                 "FixedQuadrupleListAdressInteractionTemplate at the moment"
              << std::endl;
}

template <typename _DihedralPotential>
inline real FixedQuadrupleListAdressInteractionTemplate<_DihedralPotential>::getMaxCutoff()
{
    return potential->getCutoff();
}

}  // end namespace interaction
}  // end namespace espressopp
#endif
