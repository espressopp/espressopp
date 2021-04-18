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
      template<bool VEC_MODE_SOA>
      void addForces_impl();

      template<bool VEC_MODE_SOA>
      real computeEnergy_impl();

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
      auto const modeSOA = vectorization->modeSOA();
      if(modeSOA) {
        addForces_impl<1>();
      } else {
        addForces_impl<0>();
      }
    }

    template < typename _Potential>
    template < bool VEC_MODE_SOA >
    inline void
    FixedPairListInteractionTemplate < _Potential >::addForces_impl< VEC_MODE_SOA >()
    {
      LOG4ESPP_INFO(_Potential::theLogger, "adding forces of FixedPairList");
      auto const& bc    = *getSystemRef().bc;
      auto& particles   = vectorization->particles;
      auto& fpl         = *fixedpairList;
      auto& pot         = *potential;
      real ltMaxBondSqr = fpl.getLongtimeMaxBondSqr();

      const Real3DInt *pa_pos     = particles.position.data();
      Real4D *pa_force            = particles.force.data();
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
        if(VEC_MODE_SOA){
          bc.getMinimumImageVectorBox(dist, {p_x[p1],p_y[p1],p_z[p1]}, {p_x[p2],p_y[p2],p_z[p2]});
        } else {
          bc.getMinimumImageVectorBox(dist, pa_pos[p1].to_Real3D(), pa_pos[p2].to_Real3D());
        }
        Real3D force;
        const real d = dist.sqr();
        if (d > ltMaxBondSqr) {
          fpl.setLongtimeMaxBondSqr(d);
          ltMaxBondSqr = d;
        }
        pot.computeColVarWeights(dist, bc);
        if(pot._computeForce(force, dist)){
          if(VEC_MODE_SOA){
            f_x[p1] += force.get()[0];
            f_y[p1] += force.get()[1];
            f_z[p1] += force.get()[2];
            f_x[p2] -= force.get()[0];
            f_y[p2] -= force.get()[1];
            f_z[p2] -= force.get()[2];
          } else {
            *((Real3D*)(&pa_force[p1])) += force;
            *((Real3D*)(&pa_force[p2])) -= force;
          }

          if(VEC_MODE_SOA){
            LOG4ESPP_DEBUG(_Potential::theLogger, "p" << id[p1] << "(" << p_x[p1] << "," << p_y[p1] << "," << p_z[p1] << ") "
                                               << "p" << id[p2] << "(" << p_x[p2] << "," << p_y[p2] << "," << p_z[p2] << ") "
                                               << "dist=" << sqrt(dist*dist) << " "
                                               << "force=(" << force[0] << "," << force[1] << "," << force[2] << ")" );
          } else {
            LOG4ESPP_DEBUG(_Potential::theLogger, "p" << id[p1] << "(" << pa_pos[p1].x << "," << pa_pos[p1].y << "," << pa_pos[p1].z << ") "
                                               << "p" << id[p2] << "(" << pa_pos[p2].x << "," << pa_pos[p2].y << "," << pa_pos[p2].z << ") "
                                               << "dist=" << sqrt(dist*dist) << " "
                                               << "force=(" << force[0] << "," << force[1] << "," << force[2] << ")" );
          }
        }
      }
    }

    template < typename _Potential > inline real
    FixedPairListInteractionTemplate < _Potential >::
    computeEnergy() {
      auto const modeSOA = vectorization->modeSOA();
      if(modeSOA) {
        return computeEnergy_impl<1>();
      } else {
        return computeEnergy_impl<0>();
      }
    }

    template < typename _Potential >
    template < bool VEC_MODE_SOA >
    inline real
    FixedPairListInteractionTemplate < _Potential >::
    computeEnergy_impl< VEC_MODE_SOA >() {
      LOG4ESPP_INFO(theLogger, "compute energy of the FixedPairList pairs");

      auto const& bc  = *getSystemRef().bc;
      auto& particles = vectorization->particles;
      auto& fpl         = *fixedpairList;
      auto& pot         = *potential;
      real ltMaxBondSqr = fpl.getLongtimeMaxBondSqr();

      const Real3DInt *pa_pos     = particles.position.data();
      const real* __restrict p_x  = particles.p_x.data();
      const real* __restrict p_y  = particles.p_y.data();
      const real* __restrict p_z  = particles.p_z.data();

      real e = 0.0;
      for(const auto& pair: fpl)
      {
        const auto p1 = pair.first;
        const auto p2 = pair.second;

        Real3D r21;
        if(VEC_MODE_SOA){
          bc.getMinimumImageVectorBox(r21, {p_x[p1],p_y[p1],p_z[p1]}, {p_x[p2],p_y[p2],p_z[p2]});
        } else {
          bc.getMinimumImageVectorBox(r21, pa_pos[p1].to_Real3D(), pa_pos[p2].to_Real3D());
        }
        pot.computeColVarWeights(r21, bc);
        e += pot._computeEnergy(r21);
      #if 0
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        Real3D r21;
        bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
        potential->computeColVarWeights(r21, bc);
        e += potential->_computeEnergy(r21);
      #endif
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
              //std::cout << "Warning! At the moment computeVirialX in FixedPairListInteractionTemplate does not work." << std::endl << "Therefore, the corresponding interactions won't be included in calculation." << std::endl;

       /*int i = 0;
       int bin1 = 0;
       int bin2 = 0;

       const bc::BC& bc = *getSystemRef().bc;
       Real3D Li = bc.getBoxL();
       real Delta_x = Li[0] / (real)bins;
       real Volume = Li[1] * Li[2] * Delta_x;

       size_t size = bins;
       std::vector <real> p_xx_local(size);
       for (i = 0; i < bins; ++i)
       {
          p_xx_local.at(i) = 0.0;
       }

       for (FixedPairList::PairList::Iterator it(*fixedpairList);
           it.isValid(); ++it) {
         const Particle &p1 = *it->first;
         const Particle &p2 = *it->second;
         Real3D dist(0.0,0.0,0.0);
         Real3D force(0.0,0.0,0.0);
         if(potential->_computeForce(force, p1, p2)) {
             Real3D dist = p1.position() - p2.position();
             real vir_temp = 0.5 * dist[0] * force[0];

             if (p1.position()[0] > Li[0])
             {
                  real p1_wrap = p1.position()[0] - Li[0];
                  bin1 = floor (p1_wrap / Delta_x);
             }
             else if (p1.position()[0] < 0.0)
             {
                  real p1_wrap = p1.position()[0] + Li[0];
                  bin1 = floor (p1_wrap / Delta_x);
             }
             else
             {
                  bin1 = floor (p1.position()[0] / Delta_x);
             }

             if (p2.position()[0] > Li[0])
             {
                  real p2_wrap = p2.position()[0] - Li[0];
                  bin2 = floor (p2_wrap / Delta_x);
             }
             else if (p2.position()[0] < 0.0)
             {
                  real p2_wrap = p2.position()[0] + Li[0];
                  bin2 = floor (p2_wrap / Delta_x);
             }
             else
             {
                  bin2 = floor (p2.position()[0] / Delta_x);
             }

             if (bin1 >= p_xx_local.size() || bin2 >= p_xx_local.size()){
                  std::cout << "p_xx_local.size() " << p_xx_local.size() << "\n";
                  std::cout << "bin1 " << bin1 << " bin2 " << bin2 << "\n";
                  std::cout << "p1.position()[0] " << p1.position()[0] << " p2.position()[0]" << p2.position()[0] << "\n";
                  std::cout << "FATAL ERROR: computeVirialX error" << "\n";
                  exit(0);
             }
             p_xx_local.at(bin1) += vir_temp;
             p_xx_local.at(bin2) += vir_temp;
         }
       }

       std::vector <real> p_xx_sum(size);
       for (i = 0; i < bins; ++i)
       {
           p_xx_sum.at(i) = 0.0;
           boost::mpi::all_reduce(*mpiWorld, p_xx_local.at(i), p_xx_sum.at(i), std::plus<real>());
       }

       std::transform(p_xx_sum.begin(), p_xx_sum.end(), p_xx_sum.begin(),std::bind2nd(std::divides<real>(),Volume));
       for (i = 0; i < bins; ++i)
       {
          p_xx_total.at(i) += p_xx_sum.at(i);         // TO EXCLUDE THEM
       }*/
    }


    template < typename _Potential > inline real
    FixedPairListInteractionTemplate < _Potential >::
    computeVirial() {
      LOG4ESPP_WARN(_Potential::theLogger, "Warning! "<<__FUNCTION__<<"() is not yet implemented.");

      LOG4ESPP_INFO(theLogger, "compute the virial for the FixedPair List");

      real w = 0.0;
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      auto& particles = vectorization->particles;
      for(const auto& pair: *fixedpairList)
      {
        const auto& p1 = pair.first;
        const auto& p2 = pair.second;
      #if 0
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;

        Real3D r21;
        bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
        Real3D force;
        potential->computeColVarWeights(r21, bc);
        if(potential->_computeForce(force, r21)) {
          w += r21 * force;
        }
      #endif
      }

      real wsum;
      boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
      return wsum;
    }

    template < typename _Potential > inline void
    FixedPairListInteractionTemplate < _Potential >::computeVirialTensor(Tensor& w){
      LOG4ESPP_WARN(_Potential::theLogger, "Warning! "<<__FUNCTION__<<"() is not yet implemented.");

      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");

      Tensor wlocal(0.0);
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      auto& particles = vectorization->particles;
      for(const auto& pair: *fixedpairList)
      {
        const auto& p1 = pair.first;
        const auto& p2 = pair.second;
      #if 0
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        Real3D r21;
        bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
        Real3D force;
        potential->computeColVarWeights(r21, bc);
        if(potential->_computeForce(force, r21)) {
          wlocal += Tensor(r21, force);
        }
      #endif
      }

      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
      w += wsum;
    }

    template < typename _Potential > inline void
    FixedPairListInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor& w, real z){
      LOG4ESPP_WARN(_Potential::theLogger, "Warning! "<<__FUNCTION__<<"() is not yet implemented.");

      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");

      Tensor wlocal(0.0);
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      auto& particles = vectorization->particles;
      for(const auto& pair: *fixedpairList)
      {
        const auto& p1 = pair.first;
        const auto& p2 = pair.second;
      #if 0
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        Real3D p1pos = p1.position();
        Real3D p2pos = p2.position();

        if(  (p1pos[2]>=z && p2pos[2]<=z) ||
             (p1pos[2]<=z && p2pos[2]>=z) ){
          Real3D r21;
          bc.getMinimumImageVectorBox(r21, p1pos, p2pos);
          Real3D force;
          potential->computeColVarWeights(r21, bc);
          if(potential->_computeForce(force, r21)) {
            wlocal += Tensor(r21, force);
          }
        }
      #endif
      }

      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
      w += wsum;
    }

    template < typename _Potential > inline void
    FixedPairListInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor *w, int n){
      LOG4ESPP_WARN(_Potential::theLogger, "Warning! "<<__FUNCTION__<<"() is not yet implemented.");

      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");

      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      Real3D Li = bc.getBoxL();
      Tensor *wlocal = new Tensor[n];
      for(int i=0; i<n; i++) wlocal[i] = Tensor(0.0);
      auto& particles = vectorization->particles;
      for(const auto& pair: *fixedpairList)
      {
        const auto& p1 = pair.first;
        const auto& p2 = pair.second;
      #if 0
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        Real3D p1pos = p1.position();
        Real3D p2pos = p2.position();

        int position1 = (int)( n * p1pos[2]/Li[2]);
        int position2 = (int)( n * p1pos[2]/Li[2]);

        int maxpos = std::max(position1, position2);
        int minpos = std::min(position1, position2);

        Real3D r21;
        bc.getMinimumImageVectorBox(r21, p1pos, p2pos);
        Real3D force;
        Tensor ww;
        if(potential->_computeForce(force, r21)) {
          ww = Tensor(r21, force);
        }

        int i = minpos + 1;
        while(i<=maxpos){
          wlocal[i] += ww;
          i++;
        }
      #endif
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
