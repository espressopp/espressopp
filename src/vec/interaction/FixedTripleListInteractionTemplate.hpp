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
#ifndef VEC_INTERACTION_FIXEDTRIPLELISTINTERACTIONTEMPLATE_HPP
#define VEC_INTERACTION_FIXEDTRIPLELISTINTERACTIONTEMPLATE_HPP

#include "vec/FixedTripleList.hpp"

#include "mpi.hpp"
#include "interaction/Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
// #include "FixedTripleListAdress.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "types.hpp"

namespace espressopp { namespace vec {
  namespace interaction {

    template < typename _AngularPotential >
    class FixedTripleListInteractionTemplate
      : public espressopp::interaction::Interaction, SystemAccess
    {
    protected:
      typedef _AngularPotential Potential;

    public:
      FixedTripleListInteractionTemplate(
        std::shared_ptr < System > _system,
        std::shared_ptr < vec::FixedTripleList > _fixedtripleList,
        std::shared_ptr < Potential > _potential)
        : SystemAccess(_system),
          vectorization(getSystem()->vectorization),
          fixedtripleList(_fixedtripleList),
          potential(_potential)
      {
        if (!potential) {
          LOG4ESPP_ERROR(theLogger, "NULL potential");
        }
      }

      virtual ~FixedTripleListInteractionTemplate() {};

      void setFixedTripleList(std::shared_ptr <FixedTripleList> _fixedtripleList)
      {
        fixedtripleList = _fixedtripleList;
      }

      std::shared_ptr<FixedTripleList> getFixedTripleList()
      {
        return fixedtripleList;
      }

      void setPotential(std::shared_ptr<Potential> _potential)
      {
        if (_potential) {
          potential = _potential;
        } else {
          LOG4ESPP_ERROR(theLogger, "NULL potential");
        }
      }

      std::shared_ptr<Potential> getPotential()
      {
        return potential;
      }

      virtual void addForces();
      virtual real computeEnergy();
      virtual real computeEnergyDeriv();
      virtual real computeEnergyAA();
      virtual real computeEnergyCG();
      virtual real computeEnergyAA(int atomtype);
      virtual real computeEnergyCG(int atomtype);
      virtual void computeVirialX(std::vector<real> &p_xx_total, int bins);
      virtual real computeVirial();
      virtual void computeVirialTensor(Tensor& w);
      virtual void computeVirialTensor(Tensor& w, real z);
      virtual void computeVirialTensor(Tensor *w, int n);
      virtual real getMaxCutoff();
      virtual int bondType() { return espressopp::interaction::Angular; }

    protected:
      int ntypes;
      std::shared_ptr < vec::Vectorization > vectorization;
      std::shared_ptr < FixedTripleList > fixedtripleList;
      std::shared_ptr < Potential > potential;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _AngularPotential > inline void
    FixedTripleListInteractionTemplate <_AngularPotential>::
    addForces()
    {
      LOG4ESPP_INFO(theLogger, "add forces computed by FixedTripleList");

      auto const& bc              = *getSystemRef().bc;
      auto& particles             = vectorization->particles;
      auto& ftl                   = *fixedtripleList;
      auto& pot                   = *potential;
      const real* __restrict p_x  = particles.p_x.data();
      const real* __restrict p_y  = particles.p_y.data();
      const real* __restrict p_z  = particles.p_z.data();
      real* __restrict f_x        = particles.f_x.data();
      real* __restrict f_y        = particles.f_y.data();
      real* __restrict f_z        = particles.f_z.data();
      const size_t* __restrict id = particles.id.data();

      for(const auto& triple: ftl)
      {
        const auto p1 = std::get<0>(triple);
        const auto p2 = std::get<1>(triple);
        const auto p3 = std::get<2>(triple);
        Real3D dist12, dist32;
        bc.getMinimumImageVectorBox(dist12, {p_x[p1],p_y[p1],p_z[p1]}, {p_x[p2],p_y[p2],p_z[p2]});
        bc.getMinimumImageVectorBox(dist32, {p_x[p3],p_y[p3],p_z[p3]}, {p_x[p2],p_y[p2],p_z[p2]});
        Real3D force12, force32;
        pot.computeColVarWeights(dist12, dist32, bc);
        pot._computeForce(force12, force32, dist12, dist32);
        {
          f_x[p1] += force12.get()[0];
          f_y[p1] += force12.get()[1];
          f_z[p1] += force12.get()[2];

          f_x[p2] -= force12.get()[0] + force32.get()[0];
          f_y[p2] -= force12.get()[1] + force32.get()[1];
          f_z[p2] -= force12.get()[2] + force32.get()[2];

          f_x[p3] += force32.get()[0];
          f_y[p3] += force32.get()[1];
          f_z[p3] += force32.get()[2];
        }
      }
    }

    template < typename _AngularPotential > inline real
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeEnergy()
    {
      LOG4ESPP_INFO(theLogger, "compute energy of the triples");

      auto const& bc  = *getSystemRef().bc;
      auto& pot       = *potential;
      auto& particles = vectorization->particles;

      real e = 0.0;
      for(const auto& triple: *fixedtripleList)
      {
        const auto p1 = std::get<0>(triple);
        const auto p2 = std::get<1>(triple);
        const auto p3 = std::get<2>(triple);

        const Real3D p1_pos = particles.getPosition(p1);
        const Real3D p2_pos = particles.getPosition(p2);
        const Real3D p3_pos = particles.getPosition(p3);

        const Real3D dist12 = bc.getMinimumImageVector(p1_pos, p2_pos);
        const Real3D dist32 = bc.getMinimumImageVector(p3_pos, p2_pos);

        pot.computeColVarWeights(dist12, dist32, bc);
        e += pot._computeEnergy(dist12, dist32);
      }

      real esum;
      boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
      return esum;
    }

    template < typename _AngularPotential > inline real
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeEnergyDeriv() {
      std::cout << "Warning! At the moment computeEnergyDeriv() in FixedTripleListInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _AngularPotential > inline real
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeEnergyAA() {
      std::cout << "Warning! At the moment computeEnergyAA() in FixedTripleListInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _AngularPotential > inline real
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeEnergyAA(int atomtype) {
      std::cout << "Warning! At the moment computeEnergyAA(int atomtype) in FixedTripleListInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _AngularPotential > inline real
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeEnergyCG() {
      std::cout << "Warning! At the moment computeEnergyCG() in FixedTripleListInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _AngularPotential > inline real
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeEnergyCG(int atomtype) {
      std::cout << "Warning! At the moment computeEnergyCG(int atomtype) in FixedTripleListInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _AngularPotential >
    inline void
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeVirialX(std::vector<real> &p_xx_total, int bins) {
        std::cout << "Warning! At the moment computeVirialX in FixedTripleListInteractionTemplate does not work."<< std::endl
          << "Therefore, the corresponding interactions won't be included in calculation." << std::endl;
    }

    template < typename _AngularPotential > inline real
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute scalar virial of the triples");

      auto const& bc  = *getSystemRef().bc;
      auto& pot       = *potential;
      auto& particles = vectorization->particles;

      real w = 0.0;
      for(const auto& triple: *fixedtripleList)
      {
        const auto p1 = std::get<0>(triple);
        const auto p2 = std::get<1>(triple);
        const auto p3 = std::get<2>(triple);

        const Real3D p1_pos = particles.getPosition(p1);
        const Real3D p2_pos = particles.getPosition(p2);
        const Real3D p3_pos = particles.getPosition(p3);

        Real3D dist12, dist32;
        bc.getMinimumImageVectorBox(dist12, p1_pos, p2_pos);
        bc.getMinimumImageVectorBox(dist32, p3_pos, p2_pos);
        Real3D force12, force32;
        pot.computeColVarWeights(dist12, dist32, bc);
        pot._computeForce(force12, force32, dist12, dist32);
        w += dist12 * force12 + dist32 * force32;
      }

      real wsum;
      boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
      return wsum;
    }

    template < typename _AngularPotential > inline void
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeVirialTensor(Tensor& w) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor of the triples");

      auto const& bc  = *getSystemRef().bc;
      auto& pot       = *potential;
      auto& particles = vectorization->particles;

      Tensor wlocal(0.0);
      for(const auto& triple: *fixedtripleList)
      {
        const auto p1 = std::get<0>(triple);
        const auto p2 = std::get<1>(triple);
        const auto p3 = std::get<2>(triple);

        const Real3D p1_pos = particles.getPosition(p1);
        const Real3D p2_pos = particles.getPosition(p2);
        const Real3D p3_pos = particles.getPosition(p3);

        Real3D r12, r32;
        bc.getMinimumImageVectorBox(r12, p1_pos, p2_pos);
        bc.getMinimumImageVectorBox(r32, p3_pos, p2_pos);
        Real3D force12, force32;
        pot.computeColVarWeights(r12, r32, bc);
        pot._computeForce(force12, force32, r12, r32);
        wlocal += Tensor(r12, force12) + Tensor(r32, force32);
      }

      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal,6, (double*)&wsum, std::plus<double>());
      w += wsum;
    }

    template < typename _AngularPotential > inline void
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeVirialTensor(Tensor& w, real z) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor of the triples");

      std::cout << "Warning! At the moment IK computeVirialTensor for fixed triples does'n work"<<std::endl;
    }

    template < typename _AngularPotential > inline void
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeVirialTensor(Tensor *w, int n) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor of the triples");

      std::cout << "Warning! At the moment IK computeVirialTensor for fixed triples does'n work"<<std::endl;
    }

    template < typename _AngularPotential >
    inline real
    FixedTripleListInteractionTemplate< _AngularPotential >::
    getMaxCutoff() {
      return potential->getCutoff();
    }
  }
}}

#endif//VEC_INTERACTION_FIXEDTRIPLELISTINTERACTIONTEMPLATE_HPP
