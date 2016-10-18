/*
  Copyright (C) 2012,2013,2016
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
#ifndef _INTERACTION_FIXEDQUADRUPLEANGLELISTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDQUADRUPLEANGLELISTINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "FixedQuadrupleAngleList.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "types.hpp"

namespace espressopp {
  namespace interaction {
    template < typename _DihedralPotential >
    class FixedQuadrupleAngleListInteractionTemplate : public Interaction, SystemAccess {
    protected:
      typedef _DihedralPotential Potential;
    public:
      FixedQuadrupleAngleListInteractionTemplate
      (shared_ptr < System > _system,
       shared_ptr < FixedQuadrupleAngleList > _fixedQuadrupleAngleList,
       shared_ptr < Potential > _potential)
        : SystemAccess(_system), fixedQuadrupleAngleList(_fixedQuadrupleAngleList),
          potential(_potential)
      {
        if (! potential) {
          LOG4ESPP_ERROR(theLogger, "NULL potential");
        }
      }

      void
      setFixedQuadrupleAngleList(shared_ptr < FixedQuadrupleAngleList > _fixedQuadrupleAngleList) {
        fixedQuadrupleAngleList = _fixedQuadrupleAngleList;
      }

      virtual ~FixedQuadrupleAngleListInteractionTemplate() {};

      shared_ptr < FixedQuadrupleAngleList > getFixedQuadrupleAngleList() {
        return fixedQuadrupleAngleList;
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
      virtual void computeVirialX(std::vector<real> &p_xx_total, int bins); 
      virtual real computeVirial();
      virtual void computeVirialTensor(Tensor& w);
      virtual void computeVirialTensor(Tensor& w, real z);
      virtual void computeVirialTensor(Tensor *w, int n);
      virtual real getMaxCutoff();
      virtual int bondType() { return Dihedral; }

    protected:
      int ntypes;
      shared_ptr < FixedQuadrupleAngleList > fixedQuadrupleAngleList;
      shared_ptr < Potential > potential;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _DihedralPotential > inline void
    FixedQuadrupleAngleListInteractionTemplate < _DihedralPotential >::
    addForces() {

      LOG4ESPP_INFO(theLogger, "add forces computed by FixedQuadrupleAngleList");

      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions

      for (FixedQuadrupleAngleList::QuadrupleList::Iterator it(*fixedQuadrupleAngleList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        Particle &p3 = *it->third;
        Particle &p4 = *it->fourth;

        Real3D r21, r32, r43; // 

        bc.getMinimumImageVectorBox(r21, p2.position(), p1.position());
        bc.getMinimumImageVectorBox(r32, p3.position(), p2.position());
        bc.getMinimumImageVectorBox(r43, p4.position(), p3.position());

	    Real3D force1, force2, force3, force4;  // result forces
        
        real currentAngle = fixedQuadrupleAngleList->getAngle(p1.getId(), p2.getId(),
                p3.getId(), p4.getId());

	    potential->_computeForce(force1, force2, force3, force4,
                                r21, r32, r43, currentAngle);
        p1.force() += force1;
        p2.force() += force2; //p2.force() -= force2;
        p3.force() += force3;
        p4.force() += force4;
      }
    }

    template < typename _DihedralPotential >
    inline real
    FixedQuadrupleAngleListInteractionTemplate < _DihedralPotential >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy of the quadruples");

      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      real e = 0.0;
      for (FixedQuadrupleAngleList::QuadrupleList::Iterator it(*fixedQuadrupleAngleList); it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        const Particle &p4 = *it->fourth;

        Real3D r21, r32, r43; // 

        bc.getMinimumImageVectorBox(r21, p2.position(), p1.position());
        bc.getMinimumImageVectorBox(r32, p3.position(), p2.position());
        bc.getMinimumImageVectorBox(r43, p4.position(), p3.position());

        real currentAngle = fixedQuadrupleAngleList->getAngle(p1.getId(), p2.getId(), 
                p3.getId(), p4.getId());
        
        e += potential->_computeEnergy(r21, r32, r43, currentAngle);
      }
      real esum;
      boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
      return esum;
    }
    
    template < typename _DihedralPotential > inline real
    FixedQuadrupleAngleListInteractionTemplate < _DihedralPotential >::
    computeEnergyDeriv() {
      std::cout << "Warning! At the moment computeEnergyDeriv() in FixedQuadrupleAngleListInteractionTemplate does not work." << std::endl;
      return 0.0;
    }
    
    template < typename _DihedralPotential > inline real
    FixedQuadrupleAngleListInteractionTemplate < _DihedralPotential >::
    computeEnergyAA() {
      std::cout << "Warning! At the moment computeEnergyAA() in FixedQuadrupleAngleListInteractionTemplate does not work." << std::endl;
      return 0.0;
    }
    
    template < typename _DihedralPotential > inline real
    FixedQuadrupleAngleListInteractionTemplate < _DihedralPotential >::
    computeEnergyCG() {
      std::cout << "Warning! At the moment computeEnergyCG() in FixedQuadrupleAngleListInteractionTemplate does not work." << std::endl;
      return 0.0;
    }
    
    template < typename _DihedralPotential >
    inline void
    FixedQuadrupleAngleListInteractionTemplate < _DihedralPotential >::
    computeVirialX(std::vector<real> &p_xx_total, int bins) {
        std::cout << "Warning! At the moment computeVirialX in FixedQuadrupleAngleListInteractionTemplate does not work." << std::endl << "Therefore, the corresponding interactions won't be included in calculation." << std::endl;
    }

    template < typename _DihedralPotential >
    inline real
    FixedQuadrupleAngleListInteractionTemplate < _DihedralPotential >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute scalar virial of the quadruples");

      real w = 0.0;
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedQuadrupleAngleList::QuadrupleList::Iterator it(*fixedQuadrupleAngleList); it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        const Particle &p4 = *it->fourth;

        Real3D r21, r32, r43; 

        bc.getMinimumImageVectorBox(r21, p2.position(), p1.position());
        bc.getMinimumImageVectorBox(r32, p3.position(), p2.position());
        bc.getMinimumImageVectorBox(r43, p4.position(), p3.position());

        Real3D force1, force2, force3, force4;

        real currentAngle = fixedQuadrupleAngleList->getAngle(p1.getId(), p2.getId(), 
                p3.getId(), p4.getId());
        
        potential->_computeForce(force1, force2, force3, force4,
                                r21, r32, r43, currentAngle);

        // TODO: formulas are not correct yet?

        w += r21 * force1 + r32 * force2;
      }
      
      real wsum;
      boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
      return w;
    }

    template < typename _DihedralPotential >
    inline void
    FixedQuadrupleAngleListInteractionTemplate < _DihedralPotential >::
    computeVirialTensor(Tensor& w) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor of the quadruples");
    
      Tensor wlocal(0.0);
      const bc::BC& bc = *getSystemRef().bc;

      for (FixedQuadrupleAngleList::QuadrupleList::Iterator it(*fixedQuadrupleAngleList); it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        const Particle &p4 = *it->fourth;

        Real3D r21, r32, r43; 

        bc.getMinimumImageVectorBox(r21, p2.position(), p1.position());
        bc.getMinimumImageVectorBox(r32, p3.position(), p2.position());
        bc.getMinimumImageVectorBox(r43, p4.position(), p3.position());

        Real3D force1, force2, force3, force4;

        real currentAngle = fixedQuadrupleAngleList->getAngle(p1.getId(), p2.getId(), 
                p3.getId(), p4.getId());
        
        potential->_computeForce(force1, force2, force3, force4,
                                r21, r32, r43, currentAngle);

        // TODO: formulas are not correct yet

        wlocal += Tensor(r21, force1) - Tensor(r32, force2);
      }
      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
      w += wsum;
    }


    // TODO !!!!! This doesn't work
    template < typename _DihedralPotential >
    inline void
    FixedQuadrupleAngleListInteractionTemplate < _DihedralPotential >::
    computeVirialTensor(Tensor& w, real z) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor of the quadruples");
    
      Tensor wlocal(0.0);
      const bc::BC& bc = *getSystemRef().bc;
      
      std::cout<<"Warning!!! computeVirialTensor in specified volume doesn't work for "
              "FixedQuadrupleAngleListInteractionTemplate at the moment"<<std::endl;

      for (FixedQuadrupleAngleList::QuadrupleList::Iterator it(*fixedQuadrupleAngleList); it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        const Particle &p4 = *it->fourth;

        Real3D r21, r32, r43; 

        bc.getMinimumImageVectorBox(r21, p2.position(), p1.position());
        bc.getMinimumImageVectorBox(r32, p3.position(), p2.position());
        bc.getMinimumImageVectorBox(r43, p4.position(), p3.position());

        Real3D force1, force2, force3, force4;

        real currentAngle = fixedQuadrupleAngleList->getAngle(p1.getId(), p2.getId(), 
                p3.getId(), p4.getId());
        
        potential->_computeForce(force1, force2, force3, force4,
                                r21, r32, r43, currentAngle);

        // TODO: formulas are not correct yet

        wlocal += Tensor(r21, force1) - Tensor(r32, force2);
      }
      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
      w += wsum;
    }

    
    // TODO !!!!! This doesn't work
    template < typename _DihedralPotential >
    inline void
    FixedQuadrupleAngleListInteractionTemplate < _DihedralPotential >::
    computeVirialTensor(Tensor *w, int n) {
      std::cout<<"Warning!!! computeVirialTensor in specified volume doesn't work for "
              "FixedQuadrupleAngleListInteractionTemplate at the moment"<<std::endl;
    }
    
    
    template < typename _DihedralPotential >
    inline real
    FixedQuadrupleAngleListInteractionTemplate< _DihedralPotential >::
    getMaxCutoff() {
      return potential->getCutoff();
    }
  }
}
#endif
