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
#ifndef VEC_INTERACTION_FIXEDPAIRLISTINTERACTIONTEMPLATE_HPP
#define VEC_INTERACTION_FIXEDPAIRLISTINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "interaction/Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "vec/FixedPairList.hpp"
// #include "FixedPairListAdress.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "interaction/Interaction.hpp"
#include "types.hpp"

namespace espressopp { namespace vec {
  namespace interaction {

    template < typename _Potential >
    class FixedPairListInteractionTemplate
      : public espressopp::interaction::Interaction, SystemAccess
    {

    protected:
      typedef _Potential Potential;

    public:
      FixedPairListInteractionTemplate
      (shared_ptr < vec::Vectorization > vectorization,
       shared_ptr < vec::FixedPairList > _fixedpairList,
       shared_ptr < Potential > _potential)
        : vectorization(vectorization),
          SystemAccess(vectorization->getSystem()),
          fixedpairList(_fixedpairList),
          potential(_potential)
      {
        if (! potential) {
          LOG4ESPP_ERROR(theLogger, "NULL potential");
        }
      }

      virtual ~FixedPairListInteractionTemplate() {};

      void
      setFixedPairList(shared_ptr < FixedPairList > _fixedpairList) {
        fixedpairList = _fixedpairList;
      }

      shared_ptr < FixedPairList > getFixedPairList() {
        return fixedpairList;
      }

      void
      setPotential(shared_ptr < Potential > _potential) {
        if (_potential) {
          potential = _potential;
        } else {
          LOG4ESPP_ERROR(theLogger, "NULL potential");
        }
      }

      shared_ptr < Potential > getPotential() {
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
      virtual int bondType() { return espressopp::interaction::Pair; }

    protected:
      int ntypes;
      shared_ptr < vec::Vectorization > vectorization;
      shared_ptr < vec::FixedPairList > fixedpairList;
      shared_ptr < Potential > potential;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _Potential > inline void
    FixedPairListInteractionTemplate < _Potential >::addForces() {
      LOG4ESPP_INFO(_Potential::theLogger, "adding forces of FixedPairList");
      auto const& bc    = *getSystemRef().bc;
      auto& particles   = vectorization->particles;
      auto& fpl         = *fixedpairList;
      auto& pot         = *potential;
      real ltMaxBondSqr = fpl.getLongtimeMaxBondSqr();

      const real* __restrict p_x  = particles.p_x.data();
      const real* __restrict p_y  = particles.p_y.data();
      const real* __restrict p_z  = particles.p_z.data();
      real* __restrict f_x        = particles.f_x.data();
      real* __restrict f_y        = particles.f_y.data();
      real* __restrict f_z        = particles.f_z.data();
      const size_t* __restrict id = particles.id.data();

      for(const auto& pair: fpl)
      {
        const auto p1 = pair.first;
        const auto p2 = pair.second;

        Real3D dist;
        bc.getMinimumImageVectorBox(dist, {p_x[p1],p_y[p1],p_z[p1]}, {p_x[p2],p_y[p2],p_z[p2]});
        Real3D force;
        const real d = dist.sqr();
        if (d > ltMaxBondSqr) {
          fpl.setLongtimeMaxBondSqr(d);
          ltMaxBondSqr = d;
        }
        pot.computeColVarWeights(dist, bc);
        if(pot._computeForce(force, dist)){
          f_x[p1] += force.get()[0];
          f_y[p1] += force.get()[1];
          f_z[p1] += force.get()[2];
          f_x[p2] -= force.get()[0];
          f_y[p2] -= force.get()[1];
          f_z[p2] -= force.get()[2];

          LOG4ESPP_DEBUG(_Potential::theLogger, "p" << id[p1] << "(" << p_x[p1] << "," << p_y[p1] << "," << p_z[p1] << ") "
                                              << "p" << id[p2] << "(" << p_x[p2] << "," << p_y[p2] << "," << p_z[p2] << ") "
                                              << "dist=" << sqrt(dist*dist) << " "
                                              << "force=(" << force[0] << "," << force[1] << "," << force[2] << ")" );
        }
      }
    }

    template < typename _Potential > inline real
    FixedPairListInteractionTemplate < _Potential >::
    computeEnergy()
    {
      LOG4ESPP_INFO(theLogger, "compute energy of the FixedPairList pairs");

      auto const& bc  = *getSystemRef().bc;
      auto& pot       = *potential;
      auto& particles = vectorization->particles;

      real e = 0.0;
      for(const auto& pair: *fixedpairList)
      {
        const auto p1 = pair.first;
        const auto p2 = pair.second;

        const Real3D p1pos = particles.getPosition(p1);
        const Real3D p2pos = particles.getPosition(p2);

        Real3D r21;
        bc.getMinimumImageVectorBox(r21, p1pos, p2pos);
        pot.computeColVarWeights(r21, bc);
        e += pot._computeEnergy(r21);
      }
      real esum;
      boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
      return esum;
    }

    template < typename _Potential > inline real
    FixedPairListInteractionTemplate < _Potential >::
    computeEnergyDeriv() {
      std::cout << "Warning! At the moment computeEnergyDeriv() in FixedPairListInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline real
    FixedPairListInteractionTemplate < _Potential >::
    computeEnergyAA() {
      std::cout << "Warning! At the moment computeEnergyAA() in FixedPairListInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline real
    FixedPairListInteractionTemplate < _Potential >::
    computeEnergyAA(int atomtype) {
      std::cout << "Warning! At the moment computeEnergyAA(int atomtype) in FixedPairListInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline real
    FixedPairListInteractionTemplate < _Potential >::
    computeEnergyCG() {
      std::cout << "Warning! At the moment computeEnergyCG() in FixedPairListInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline real
    FixedPairListInteractionTemplate < _Potential >::
    computeEnergyCG(int atomtype) {
      std::cout << "Warning! At the moment computeEnergyCG(int atomtype) in FixedPairListInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential >
    inline void
    FixedPairListInteractionTemplate < _Potential >::
    computeVirialX(std::vector<real> &p_xx_total, int bins) {
      LOG4ESPP_INFO(theLogger, "compute virial p_xx of the pressure tensor slabwise");
      std::cout << "Warning! At the moment computeVirialX in FixedPairListInteractionTemplate does not work." << std::endl << "Therefore, the corresponding interactions won't be included in calculation." << std::endl;
    }

    template < typename _Potential > inline real
    FixedPairListInteractionTemplate < _Potential >::
    computeVirial()
    {
      LOG4ESPP_INFO(theLogger, "compute the virial for the FixedPair List");

      auto const& bc  = *getSystemRef().bc;
      auto& pot       = *potential;
      auto& particles = vectorization->particles;

      real w = 0.0;
      for(const auto& pair: *fixedpairList)
      {
        const auto& p1 = pair.first;
        const auto& p2 = pair.second;

        const Real3D p1pos = particles.getPosition(p1);
        const Real3D p2pos = particles.getPosition(p2);

        Real3D r21;
        bc.getMinimumImageVectorBox(r21, p1pos, p2pos);

        Real3D force;
        pot.computeColVarWeights(r21, bc);
        if(pot._computeForce(force, r21)) {
          w += r21 * force;
        }
      }

      real wsum;
      boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
      return wsum;
    }

    template < typename _Potential > inline void
    FixedPairListInteractionTemplate < _Potential >::computeVirialTensor(Tensor& w)
    {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");

      auto const& bc  = *getSystemRef().bc;
      auto& pot       = *potential;
      auto& particles = vectorization->particles;

      Tensor wlocal(0.0);
      for(const auto& pair: *fixedpairList)
      {
        const auto& p1 = pair.first;
        const auto& p2 = pair.second;

        const Real3D p1pos = particles.getPosition(p1);
        const Real3D p2pos = particles.getPosition(p2);

        Real3D r21;
        bc.getMinimumImageVectorBox(r21, p1pos, p2pos);

        Real3D force;
        pot.computeColVarWeights(r21, bc);
        if(pot._computeForce(force, r21)) {
          wlocal += Tensor(r21, force);
        }
      }

      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
      w += wsum;
    }

    template < typename _Potential > inline void
    FixedPairListInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor& w, real z)
    {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");

      auto const& bc  = *getSystemRef().bc;
      auto& pot       = *potential;
      auto& particles = vectorization->particles;

      Tensor wlocal(0.0);
      for(const auto& pair: *fixedpairList)
      {
        const auto& p1 = pair.first;
        const auto& p2 = pair.second;

        const Real3D p1pos = particles.getPosition(p1);
        const Real3D p2pos = particles.getPosition(p2);

        if(  (p1pos[2]>=z && p2pos[2]<=z) ||
             (p1pos[2]<=z && p2pos[2]>=z) ){
          Real3D r21;
          bc.getMinimumImageVectorBox(r21, p1pos, p2pos);
          Real3D force;
          pot.computeColVarWeights(r21, bc);
          if(pot._computeForce(force, r21)) {
            wlocal += Tensor(r21, force);
          }
        }
      }

      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
      w += wsum;
    }

    template < typename _Potential > inline void
    FixedPairListInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor *w, int n)
    {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");

      auto const& bc  = *getSystemRef().bc;
      auto& pot       = *potential;
      auto& particles = vectorization->particles;

      Real3D Li = bc.getBoxL();
      Tensor *wlocal = new Tensor[n];
      for(int i=0; i<n; i++) wlocal[i] = Tensor(0.0);

      for(const auto& pair: *fixedpairList)
      {
        const auto& p1 = pair.first;
        const auto& p2 = pair.second;

        const Real3D p1pos = particles.getPosition(p1);
        const Real3D p2pos = particles.getPosition(p2);

        int position1 = (int)( n * p1pos[2]/Li[2]);
        int position2 = (int)( n * p1pos[2]/Li[2]);

        int maxpos = std::max(position1, position2);
        int minpos = std::min(position1, position2);

        Real3D r21;
        bc.getMinimumImageVectorBox(r21, p1pos, p2pos);
        Real3D force;
        Tensor ww;
        if(pot._computeForce(force, r21)) {
          ww = Tensor(r21, force);
        }

        int i = minpos + 1;
        while(i<=maxpos){
          wlocal[i] += ww;
          i++;
        }
      }

      Tensor *wsum = new Tensor[n];
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, n, (double*)&wsum, std::plus<double>());

      for(int j=0; j<n; j++){
        w[j] += wsum[j];
      }

      delete [] wsum;
      delete [] wlocal;
    }

    template < typename _Potential >
    inline real
    FixedPairListInteractionTemplate< _Potential >::
    getMaxCutoff() {
      return potential->getCutoff();
    }
  }
}}

#endif//VEC_INTERACTION_FIXEDPAIRLISTINTERACTIONTEMPLATE_HPP
