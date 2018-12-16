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
#ifndef _INTERACTION_FIXEDTRIPLELISTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDTRIPLELISTINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "FixedTripleList.hpp"
#include "FixedTripleListAdress.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "types.hpp"

namespace espressopp {
  namespace interaction {
    template < typename _AngularPotential >
    class FixedTripleListInteractionTemplate : public Interaction, SystemAccess {

    protected:
      typedef _AngularPotential Potential;

    public:
      FixedTripleListInteractionTemplate
      (shared_ptr < System > _system,
       shared_ptr < FixedTripleList > _fixedtripleList,
       shared_ptr < Potential > _potential)
        : SystemAccess(_system), fixedtripleList(_fixedtripleList),
          potential(_potential)
      {
          if (! potential) {
                LOG4ESPP_ERROR(theLogger, "NULL potential");
          }
        //potentialArray = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());
      }

      virtual ~FixedTripleListInteractionTemplate() {};

      void
      setFixedTripleList(shared_ptr < FixedTripleList > _fixedtripleList) {
        fixedtripleList = _fixedtripleList;
      }

      shared_ptr < FixedTripleList > getFixedTripleList() {
        return fixedtripleList;
      }

      /*void
      setPotential(int type1, int type2, const Potential &potential) {
        potentialArray.at(type1, type2) = potential;
      }*/
      void
      setPotential(shared_ptr < Potential> _potential) {
         if (_potential) {
            potential = _potential;
         } else {
            LOG4ESPP_ERROR(theLogger, "NULL potential");
         }
      }

      /*Potential &getPotential(int type1, int type2) {
        return potentialArray.at(0, 0);
      }*/

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
      virtual int bondType() { return Angular; }

    protected:
      int ntypes;
      shared_ptr<FixedTripleList> fixedtripleList;
      //esutil::Array2D<Potential, esutil::enlarge> potentialArray;
      shared_ptr < Potential > potential;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _AngularPotential > inline void
    FixedTripleListInteractionTemplate <_AngularPotential>::
    addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by FixedTripleList");
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedTripleList::TripleList::Iterator it(*fixedtripleList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        Particle &p3 = *it->third;
        //const Potential &potential = getPotential(p1.type(), p2.type());
        Real3D dist12, dist32;
        bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
        bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
        Real3D force12, force32;
        potential->computeColVarWeights(dist12, dist32, bc);
        potential->_computeForce(force12, force32, dist12, dist32);
        p1.force() += force12;
        p2.force() -= force12 + force32;
        p3.force() += force32;
      }
    }

    template < typename _AngularPotential > inline real
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy of the triples");

      const bc::BC& bc = *getSystemRef().bc;
      real e = 0.0;
      for (FixedTripleList::TripleList::Iterator it(*fixedtripleList); it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        //const Potential &potential = getPotential(p1.type(), p2.type());
        Real3D dist12 = bc.getMinimumImageVector(p1.position(), p2.position());
        Real3D dist32 = bc.getMinimumImageVector(p3.position(), p2.position());
        potential->computeColVarWeights(dist12, dist32, bc);
        e += potential->_computeEnergy(dist12, dist32);
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
        std::cout << "Warning! At the moment computeVirialX in FixedTripleListInteractionTemplate does not work." << std::endl << "Therefore, the corresponding interactions won't be included in calculation." << std::endl;
    }

    template < typename _AngularPotential > inline real
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute scalar virial of the triples");

      const bc::BC& bc = *getSystemRef().bc;
      real w = 0.0;
      for (FixedTripleList::TripleList::Iterator it(*fixedtripleList); it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        //const Potential &potential = getPotential(p1.type(), p2.type());
        const espressopp::bc::BC& bc = *getSystemRef().bc;
        Real3D dist12, dist32;
        bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
        bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
        Real3D force12, force32;
        potential->computeColVarWeights(dist12, dist32, bc);
        potential->_computeForce(force12, force32, dist12, dist32);
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

      Tensor wlocal(0.0);
      const bc::BC& bc = *getSystemRef().bc;
      for (FixedTripleList::TripleList::Iterator it(*fixedtripleList); it.isValid(); ++it){
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        //const Potential &potential = getPotential(0, 0);
        Real3D r12, r32;
        bc.getMinimumImageVectorBox(r12, p1.position(), p2.position());
        bc.getMinimumImageVectorBox(r32, p3.position(), p2.position());
        Real3D force12, force32;
        potential->computeColVarWeights(r12, r32, bc);
        potential->_computeForce(force12, force32, r12, r32);
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

      /*
      Tensor wlocal(0.0);
      const bc::BC& bc = *getSystemRef().bc;
      for (FixedTripleList::TripleList::Iterator it(*fixedtripleList); it.isValid(); ++it){
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;

        Real3D p2pos = p2.position();

        if(  (p2pos[0]>xmin && p2pos[0]<xmax &&
              p2pos[1]>ymin && p2pos[1]<ymax &&
              p2pos[2]>zmin && p2pos[2]<zmax) ){
          Real3D r12, r32;
          bc.getMinimumImageVectorBox(r12, p1.position(), p2.position());
          bc.getMinimumImageVectorBox(r32, p3.position(), p2.position());
          Real3D force12, force32;
          potential->_computeForce(force12, force32, r12, r32);
          wlocal += Tensor(r12, force12) + Tensor(r32, force32);
        }
      }

      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal,6, (double*)&wsum, std::plus<double>());
      w += wsum;
       */
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
      /*real cutoff = 0.0;
      for (int i = 0; i < ntypes; i++) {
        for (int j = 0; j < ntypes; j++) {
          cutoff = std::max(cutoff, getPotential(i, j).getCutoff());
        }
      }*/
      return potential->getCutoff();
    }
  }
}
#endif
