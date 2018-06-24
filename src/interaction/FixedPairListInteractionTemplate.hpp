/*
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
#ifndef _INTERACTION_FIXEDPAIRLISTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDPAIRLISTINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "Interaction.hpp"
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

namespace espressopp {
  namespace interaction {
    template < typename _Potential >
    class FixedPairListInteractionTemplate: public Interaction, SystemAccess {

    protected:
      typedef _Potential Potential;

    public:
      FixedPairListInteractionTemplate
      (shared_ptr < System > system,
       shared_ptr < FixedPairList > _fixedpairList,
       shared_ptr < Potential > _potential)
        : SystemAccess(system), fixedpairList(_fixedpairList),
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
      setPotential(shared_ptr < Potential> _potential) {
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
      virtual int bondType() { return Pair; }

    protected:
      int ntypes;
      shared_ptr < FixedPairList > fixedpairList;
      shared_ptr < Potential > potential;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _Potential > inline void
    FixedPairListInteractionTemplate < _Potential >::addForces() {
      LOG4ESPP_INFO(_Potential::theLogger, "adding forces of FixedPairList");
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      real ltMaxBondSqr = fixedpairList->getLongtimeMaxBondSqr();
      for (FixedPairList::PairList::Iterator it(*fixedpairList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        Real3D dist;
        bc.getMinimumImageVectorBox(dist, p1.position(), p2.position());
        Real3D force;
        real d = dist.sqr();
        if (d > ltMaxBondSqr) {
        	fixedpairList->setLongtimeMaxBondSqr(d);
        	ltMaxBondSqr = d;
        }
        if(potential->_computeForce(force, dist)) {
          p1.force() += force;
          p2.force() -= force;
          LOG4ESPP_DEBUG(_Potential::theLogger, "p" << p1.id() << "(" << p1.position()[0] << "," << p1.position()[1] << "," << p1.position()[2] << ") "
        		                             << "p" << p2.id() << "(" << p2.position()[0] << "," << p2.position()[1] << "," << p2.position()[2] << ") "
        		                             << "dist=" << sqrt(dist*dist) << " "
        		                             << "force=(" << force[0] << "," << force[1] << "," << force[2] << ")" );
        }
      }
    }

    template < typename _Potential > inline real
    FixedPairListInteractionTemplate < _Potential >::
    computeEnergy() {

      LOG4ESPP_INFO(theLogger, "compute energy of the FixedPairList pairs");

      real e = 0.0;
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedPairList::PairList::Iterator it(*fixedpairList);
	   it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        Real3D r21;
        bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
        e += potential->_computeEnergy(r21);
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
      LOG4ESPP_INFO(theLogger, "compute the virial for the FixedPair List");

      real w = 0.0;
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedPairList::PairList::Iterator it(*fixedpairList);
           it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;

        Real3D r21;
        bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
        Real3D force;
        if(potential->_computeForce(force, r21)) {
          w += r21 * force;
        }
      }

      real wsum;
      boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
      return wsum;
    }

    template < typename _Potential > inline void
    FixedPairListInteractionTemplate < _Potential >::computeVirialTensor(Tensor& w){
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");

      Tensor wlocal(0.0);
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedPairList::PairList::Iterator it(*fixedpairList);
           it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        Real3D r21;
        bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
        Real3D force;
        if(potential->_computeForce(force, r21)) {
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
    computeVirialTensor(Tensor& w, real z){
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");

      Tensor wlocal(0.0);
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedPairList::PairList::Iterator it(*fixedpairList);
           it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        Real3D p1pos = p1.position();
        Real3D p2pos = p2.position();

        if(  (p1pos[2]>=z && p2pos[2]<=z) ||
             (p1pos[2]<=z && p2pos[2]>=z) ){
          Real3D r21;
          bc.getMinimumImageVectorBox(r21, p1pos, p2pos);
          Real3D force;
          if(potential->_computeForce(force, r21)) {
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
    computeVirialTensor(Tensor *w, int n){
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");

      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      Real3D Li = bc.getBoxL();
      Tensor *wlocal = new Tensor[n];
      for(int i=0; i<n; i++) wlocal[i] = Tensor(0.0);
      for (FixedPairList::PairList::Iterator it(*fixedpairList);
           it.isValid(); ++it) {
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
}
#endif
