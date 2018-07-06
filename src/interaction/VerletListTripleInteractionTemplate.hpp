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
#ifndef _INTERACTION_VERLETLISTTRIPLEINTERACTIONTEMPLATE_HPP
#define _INTERACTION_VERLETLISTTRIPLEINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "types.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "VerletListTriple.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "esutil/Array3D.hpp"

namespace espressopp {
  namespace interaction {
    template < typename _ThreeBodyPotential >
    class VerletListTripleInteractionTemplate : public Interaction, SystemAccess {

    protected:
      typedef _ThreeBodyPotential Potential;

    public:
      VerletListTripleInteractionTemplate
      (shared_ptr < System > _system,shared_ptr < VerletListTriple > _verletlisttriple)
      : SystemAccess(_system), verletListTriple(_verletlisttriple)
      {
        potentialArray = esutil::Array3D<Potential, esutil::enlarge>(0, 0, 0, Potential());
        ntypes = 0;
      }

      void
      setVerletListTriple(shared_ptr < VerletListTriple > _verletlisttriple) {
        verletListTriple = _verletlisttriple;
      }
      shared_ptr < VerletListTriple > getVerletListTriple() {
        return verletListTriple;
      }


      /* it will create the table for every possible combination
       * For example, type1 = 0, type2 = 1:
       * 000, 100, 001, 010, 101, 111
       */
      void
      setPotential(int type1, int type2, int type3, const Potential &potential) {
        ntypes = std::max(ntypes, std::max(type1+1, std::max(type2+1, type3+1)));

        potentialArray.at(type1, type2, type3) = potential;
      }

      Potential &getPotential(int type1, int type2, int type3) {
        return potentialArray.at(type1, type2, type3);
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
      shared_ptr<VerletListTriple> verletListTriple;
      esutil::Array3D<Potential, esutil::enlarge> potentialArray;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _ThreeBodyPotential > inline void
    VerletListTripleInteractionTemplate <_ThreeBodyPotential>::
    addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by VerletListTriple");
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions

      for (TripleList::Iterator it(verletListTriple->getTriples()); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second; // the main particle
        Particle &p3 = *it->third;
        Real3D r12, r32;
        bc.getMinimumImageVectorBox(r12, p1.position(), p2.position());
        bc.getMinimumImageVectorBox(r32, p3.position(), p2.position());

        int type1 = p1.type();
        int type2 = p2.type();
        int type3 = p3.type();
        const Potential &potential = getPotential(type1, type2, type3);

        Real3D force12(0.0,0.0,0.0), force32(0.0,0.0,0.0);

        if(potential._computeForce(force12, force32, r12, r32)){
          p1.force() += force12;
          p2.force() -= force12 + force32;
          p3.force() += force32;
        }
      }
    }

    template < typename _ThreeBodyPotential > inline real
    VerletListTripleInteractionTemplate < _ThreeBodyPotential >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy of the triples");

      const bc::BC& bc = *getSystemRef().bc;
      real e = 0.0;
      for (TripleList::Iterator it(verletListTriple->getTriples()); it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        Real3D r12 = bc.getMinimumImageVector(p1.position(), p2.position());
        Real3D r32 = bc.getMinimumImageVector(p3.position(), p2.position());

        int type1 = p1.type();
        int type2 = p2.type();
        int type3 = p3.type();
        const Potential &potential = getPotential(type1, type2, type3);

        e += potential._computeEnergy(r12, r32);
      }
      real esum;
      boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
      return esum;
    }

    template < typename _ThreeBodyPotential > inline real
    VerletListTripleInteractionTemplate < _ThreeBodyPotential >::
    computeEnergyDeriv() {
      std::cout << "Warning! At the moment computeEnergyDeriv() in VerletListTripleInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _ThreeBodyPotential > inline real
    VerletListTripleInteractionTemplate < _ThreeBodyPotential >::
    computeEnergyAA() {
      std::cout << "Warning! At the moment computeEnergyAA() in VerletListTripleInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _ThreeBodyPotential > inline real
    VerletListTripleInteractionTemplate < _ThreeBodyPotential >::
    computeEnergyAA(int atomtype) {
      std::cout << "Warning! At the moment computeEnergyAA(int atomtype) in VerletListTripleInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _ThreeBodyPotential > inline real
    VerletListTripleInteractionTemplate < _ThreeBodyPotential >::
    computeEnergyCG() {
      std::cout << "Warning! At the moment computeEnergyCG() in VerletListTripleInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _ThreeBodyPotential > inline real
    VerletListTripleInteractionTemplate < _ThreeBodyPotential >::
    computeEnergyCG(int atomtype) {
      std::cout << "Warning! At the moment computeEnergyCG(int atomtype) in VerletListTripleInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _ThreeBodyPotential >
    inline void
    VerletListTripleInteractionTemplate < _ThreeBodyPotential >::
    computeVirialX(std::vector<real> &p_xx_total, int bins) {
        std::cout << "Warning! At the moment computeVirialX in VerletListTripleInteractionTemplate does not work." << std::endl << "Therefore, the corresponding interactions won't be included in calculation." << std::endl;
    }

    template < typename _ThreeBodyPotential > inline real
    VerletListTripleInteractionTemplate < _ThreeBodyPotential >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute scalar virial of the triples");

      const bc::BC& bc = *getSystemRef().bc;
      real w = 0.0;
      for (TripleList::Iterator it(verletListTriple->getTriples()); it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        //const Potential &potential = getPotential(p1.type(), p2.type());
        const espressopp::bc::BC& bc = *getSystemRef().bc;
        Real3D dist12, dist32;
        bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
        bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());

        int type1 = p1.type();
        int type2 = p2.type();
        int type3 = p3.type();
        const Potential &potential = getPotential(type1, type2, type3);

        Real3D force12(0.0,0.0,0.0), force32(0.0,0.0,0.0);
        if(potential._computeForce(force12, force32, dist12, dist32)){
          w += dist12 * force12 + dist32 * force32;
        }
      }

      real wsum;
      boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
      return wsum;
    }

    template < typename _ThreeBodyPotential > inline void
    VerletListTripleInteractionTemplate < _ThreeBodyPotential >::
    computeVirialTensor(Tensor& w) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor of the triples");

      Tensor wlocal(0.0);
      const bc::BC& bc = *getSystemRef().bc;
      for (TripleList::Iterator it(verletListTriple->getTriples()); it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        Real3D r12, r32;
        bc.getMinimumImageVectorBox(r12, p1.position(), p2.position());
        bc.getMinimumImageVectorBox(r32, p3.position(), p2.position());

        int type1 = p1.type();
        int type2 = p2.type();
        int type3 = p3.type();
        const Potential &potential = getPotential(type1, type2, type3);

        Real3D force12(0.0,0.0,0.0), force32(0.0,0.0,0.0);
        if(potential._computeForce(force12, force32, r12, r32)){
          wlocal += Tensor(r12, force12) + Tensor(r32, force32);
        }
      }

      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
      w += wsum;
    }


    // Irvin-Kirkwood method
    template < typename _ThreeBodyPotential > inline void
    VerletListTripleInteractionTemplate < _ThreeBodyPotential >::
    computeVirialTensor(Tensor& w, real z) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor of the triples");

      std::cout << "At the moment IK computeVirialTensor for triples does'n work"<<std::endl;

      /*
      Tensor wlocal(0.0);
      const bc::BC& bc = *getSystemRef().bc;
      for (TripleList::Iterator it(verletListTriple->getTriples()); it.isValid(); ++it) {
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

          int type1 = p1.type();
          int type2 = p2.type();
          int type3 = p3.type();
          const Potential &potential = getPotential(type1, type2, type3);

          Real3D force12(0.0,0.0,0.0), force32(0.0,0.0,0.0);
          if(potential._computeForce(force12, force32, r12, r32)){
            wlocal += Tensor(r12, force12) + Tensor(r32, force32);
          }
        }
      }

      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
      w += wsum;
       */
    }
    // Irvin-Kirkwood method
    template < typename _ThreeBodyPotential > inline void
    VerletListTripleInteractionTemplate < _ThreeBodyPotential >::
    computeVirialTensor(Tensor *w, int n) {
      std::cout << "At the moment IK computeVirialTensor for triples does'n work"<<std::endl;
    }

    template < typename _ThreeBodyPotential >
    inline real
    VerletListTripleInteractionTemplate< _ThreeBodyPotential >::
    getMaxCutoff() {
      real cutoff = 0.0;
      for (int i = 0; i < ntypes; i++) {
        for (int j = 0; j < ntypes; j++) {
          for (int k = 0; k < ntypes; k++) {
            cutoff = std::max(cutoff, getPotential(i, j, k).getCutoff());
          }
        }
      }
      return cutoff;
    }
  }
}
#endif
