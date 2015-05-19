/*
  Copyright (C) 2012,2013,2014,2015
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
#ifndef _INTERACTION_FIXEDTRIPLELISTPIADRESSINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDTRIPLELISTPIADRESSINTERACTIONTEMPLATE_HPP

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
#include "FixedTupleListAdress.hpp"

namespace espressopp {
  namespace interaction {
    template < typename _AngularPotential >
    class FixedTripleListPIadressInteractionTemplate : public Interaction, SystemAccess {
        
    protected:
      typedef _AngularPotential Potential;
      
    public:
      FixedTripleListPIadressInteractionTemplate
      (shared_ptr < System > _system,
       shared_ptr < FixedTripleList > _fixedtripleList,
       shared_ptr <FixedTupleListAdress> _fixedtupleList,
       shared_ptr < Potential > _potential,
       int _ntrotter)
        : SystemAccess(_system), fixedtripleList(_fixedtripleList), fixedtupleList(_fixedtupleList),
          potential(_potential),  ntrotter(_ntrotter)
      {
          if (! potential) {
                LOG4ESPP_ERROR(theLogger, "NULL potential");
          }
        //potentialArray = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());
      }

      virtual ~FixedTripleListPIadressInteractionTemplate() {};

      void
      setFixedTripleList(shared_ptr < FixedTripleList > _fixedtripleList) {
        fixedtripleList = _fixedtripleList;
      }

      shared_ptr < FixedTripleList > getFixedTripleList() {
        return fixedtripleList;
      }
      
      void
      setFixedTupleList(shared_ptr<FixedTupleListAdress> _fixedtupleList) {
          fixedtupleList = _fixedtupleList;
      }
      
      void
      setNTrotter(int _ntrotter) {
          ntrotter = _ntrotter;
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
      virtual real computeEnergyAA();
      virtual real computeEnergyCG();      
      virtual void computeVirialX(std::vector<real> &p_xx_total, int bins); 
      virtual real computeVirial();
      virtual void computeVirialTensor(Tensor& w);
      virtual void computeVirialTensor(Tensor& w, real z);
      virtual void computeVirialTensor(Tensor *w, int n);
      virtual real getMaxCutoff();
      virtual int bondType() { return Angular; }

    protected:
      int ntypes;
      int ntrotter; // Trotter number
      shared_ptr<FixedTripleList> fixedtripleList;
      shared_ptr<FixedTupleListAdress> fixedtupleList;
      //esutil::Array2D<Potential, esutil::enlarge> potentialArray;
      shared_ptr < Potential > potential;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _AngularPotential > inline void
    FixedTripleListPIadressInteractionTemplate <_AngularPotential>::
    addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by FixedTripleList");
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedTripleList::TripleList::Iterator it(*fixedtripleList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        Particle &p3 = *it->third;
        
        //weights 
        real w1 = p1.lambda();               
        real w2 = p2.lambda();        
        real w3 = p3.lambda();        
        
        // Completely in classical region?
        if ( (w1 < 0.0000001) && (w2 < 0.0000001) && (w3 < 0.0000001) ) {
            Real3D dist12, dist32;
            bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
            bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
            Real3D force12, force32;
            potential->_computeForce(force12, force32, dist12, dist32);
            p1.force() += force12;
            p2.force() -= force12 + force32;
            p3.force() += force32;
        }
        // Otherwise...
        else{ 
            // Get the corresponding tuples
            FixedTupleListAdress::iterator it4;
            FixedTupleListAdress::iterator it5;
            FixedTupleListAdress::iterator it6;
            it4 = fixedtupleList->find(&p1);
            it5 = fixedtupleList->find(&p2);
            it6 = fixedtupleList->find(&p3);
            
            if (it4 != fixedtupleList->end() && it5 != fixedtupleList->end() && it6 != fixedtupleList->end()) {
            
                // Get the PI bead lists (i.e. the AdResS particles)
                std::vector<Particle*> atList1;
                std::vector<Particle*> atList2;
                std::vector<Particle*> atList3;
                atList1 = it4->second;
                atList2 = it5->second;
                atList3 = it6->second;
            
                // Iterate the two iterators in a parallel fashion
                std::vector<Particle*>::iterator itv2 = atList2.begin();
                std::vector<Particle*>::iterator itv3 = atList3.begin();
                for (std::vector<Particle*>::iterator itv = atList1.begin();
                     itv != atList1.end(); ++itv) {
                
                     // they should be the same length... Total Trotter beads the same everywhere in the system
                     if (itv2 == atList2.end() || itv3 == atList3.end()){
                         std::cout << "Tuplelists seem to have different lengths or not started properly. Corresponding to CG particle " << 
                         p1.id() << "\n";
                         exit(1);
                         return;
                     }

                     // Get the individual PI beads
                     Particle &p4 = **itv;
                     Particle &p5 = **itv2;
                     Particle &p6 = **itv3;

                     // the beads we get should have the same Trotter bead number to interact with each other
                     if (p4.pib() != p5.pib() || p5.pib() != p6.pib()){
                         std::cout << "Path Integral Beads numbers do not correspond in VerletListPIadressInteractionTemplate for particles " << 
                         p4.id() << " and " << p5.id() << " and " << p6.id() << "\n";
                         exit(1);
                         return;
                     }
                     
                     // Calculate forces
                     Real3D dist12, dist32;
                     bc.getMinimumImageVectorBox(dist12, p4.position(), p5.position());
                     bc.getMinimumImageVectorBox(dist32, p6.position(), p5.position());
                     Real3D force12, force32;
                     potential->_computeForce(force12, force32, dist12, dist32);
                     force12 *= 1.0/ntrotter;
                     force32 *= 1.0/ntrotter;
                     p4.force() += force12;
                     p5.force() -= force12 + force32;
                     p6.force() += force32;                
                     
                     //Iterate the second and third iterator
                     ++itv2;
                     ++itv3;
                     
                }               
            }
            else { // this should not happen
               std::cout << " one of the VP particles not found in tuples: " << p1.id() << "-" <<
                     p1.ghost() << ", " << p2.id() << "-" << p2.ghost() << ", " << p3.id() << "-" << p3.ghost();
               std::cout << " (" << p1.position() << ") (" << p2.position() << ") (" << p3.position() << ") \n";
               exit(1);
               return;
            }
        
        }

        
      }
    }

    template < typename _AngularPotential > inline real
    FixedTripleListPIadressInteractionTemplate < _AngularPotential >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy of the triples");

      const bc::BC& bc = *getSystemRef().bc;
      real e = 0.0;
      for (FixedTripleList::TripleList::Iterator it(*fixedtripleList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        Particle &p3 = *it->third;
        
        //weights 
        real w1 = p1.lambda();               
        real w2 = p2.lambda();        
        real w3 = p3.lambda();        
        
        // Completely in classical region?
        if ( (w1 < 0.0000001) && (w2 < 0.0000001) && (w3 < 0.0000001) ) {
            Real3D dist12 = bc.getMinimumImageVector(p1.position(), p2.position());
            Real3D dist32 = bc.getMinimumImageVector(p3.position(), p2.position());
            e += potential->_computeEnergy(dist12, dist32);
        }
        // Otherwise...
        else{  
            // Get the corresponding tuples
            FixedTupleListAdress::iterator it4;
            FixedTupleListAdress::iterator it5;
            FixedTupleListAdress::iterator it6;
            it4 = fixedtupleList->find(&p1);
            it5 = fixedtupleList->find(&p2);
            it6 = fixedtupleList->find(&p3);
            
            if (it4 != fixedtupleList->end() && it5 != fixedtupleList->end() && it6 != fixedtupleList->end()) {
            
                // Get the PI bead lists (i.e. the AdResS particles)
                std::vector<Particle*> atList1;
                std::vector<Particle*> atList2;
                std::vector<Particle*> atList3;
                atList1 = it4->second;
                atList2 = it5->second;
                atList3 = it6->second;
            
                // Iterate the two iterators in a parallel fashion
                std::vector<Particle*>::iterator itv2 = atList2.begin();
                std::vector<Particle*>::iterator itv3 = atList3.begin();
                for (std::vector<Particle*>::iterator itv = atList1.begin();
                     itv != atList1.end(); ++itv) {
                
                     // they should be the same length... Total Trotter beads the same everywhere in the system
                     if (itv2 == atList2.end() || itv3 == atList3.end()){
                         std::cout << "Tuplelists seem to have different lengths or not started properly. Corresponding to CG particle " << 
                         p1.id() << "\n";
                         exit(1);
                         return 0.0;
                     }

                     // Get the individual PI beads
                     Particle &p4 = **itv;
                     Particle &p5 = **itv2;
                     Particle &p6 = **itv3;

                     // the beads we get should have the same Trotter bead number to interact with each other
                     if (p4.pib() != p5.pib() || p5.pib() != p6.pib()){
                         std::cout << "Path Integral Beads numbers do not correspond in VerletListPIadressInteractionTemplate for particles " << 
                         p4.id() << " and " << p5.id() << " and " << p6.id() << "\n";
                         exit(1);
                         return 0.0;
                     }
                     
                     // Calculate energies
                     Real3D dist12 = bc.getMinimumImageVector(p4.position(), p5.position());
                     Real3D dist32 = bc.getMinimumImageVector(p6.position(), p5.position());
                     e += (1.0/ntrotter)*potential->_computeEnergy(dist12, dist32);
                     
                     // Iterate the second and third iterator
                     ++itv2;
                     ++itv3;
                }               
            }
            else { // this should not happen
               std::cout << " one of the VP particles not found in tuples: " << p1.id() << "-" <<
                     p1.ghost() << ", " << p2.id() << "-" << p2.ghost() << ", " << p3.id() << "-" << p3.ghost();
               std::cout << " (" << p1.position() << ") (" << p2.position() << ") (" << p3.position() << ") \n";
               exit(1);
               return 0.0;
            }
        
        }
       
        //const Potential &potential = getPotential(p1.type(), p2.type());
      }
      real esum;
      boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
      return esum;
    }
    
    template < typename _AngularPotential > inline real
    FixedTripleListPIadressInteractionTemplate < _AngularPotential >::
    computeEnergyAA() {
      std::cout << "Warning! At the moment computeEnergyAA() in FixedTripleListPIadressInteractionTemplate does not work." << std::endl;
      exit(1);
      return 0.0;
    }
    
    template < typename _AngularPotential > inline real
    FixedTripleListPIadressInteractionTemplate < _AngularPotential >::
    computeEnergyCG() {
      std::cout << "Warning! At the moment computeEnergyCG() in FixedTripleListPIadressInteractionTemplate does not work." << std::endl;
      exit(1);
      return 0.0;
    }
           
    template < typename _AngularPotential >
    inline void
    FixedTripleListPIadressInteractionTemplate < _AngularPotential >::
    computeVirialX(std::vector<real> &p_xx_total, int bins) {
        std::cout << "Warning! At the moment computeVirialX in FixedTripleListPIadressInteractionTemplate does not work." << std::endl;
        exit(1);
        return;
    }

    template < typename _AngularPotential > inline real
    FixedTripleListPIadressInteractionTemplate < _AngularPotential >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute scalar virial of the triples");
      std::cout << "Warning! At the moment computeVirial in FixedTripleListPIadressInteractionTemplate does not work." << std::endl;
      exit(1);
      return 0.0;
        
      /*const bc::BC& bc = *getSystemRef().bc;
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
        potential->_computeForce(force12, force32, dist12, dist32);
        w += dist12 * force12 + dist32 * force32;
      }
      real wsum;
      boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
      return wsum;*/
    }

    template < typename _AngularPotential > inline void
    FixedTripleListPIadressInteractionTemplate < _AngularPotential >::
    computeVirialTensor(Tensor& w) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor of the triples");
      std::cout << "Warning! At the moment computeVirialTensor in FixedTripleListPIadressInteractionTemplate does not work." << std::endl;
      exit(1);
      return;
      
      /*Tensor wlocal(0.0);
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
        potential->_computeForce(force12, force32, r12, r32);
        wlocal += Tensor(r12, force12) + Tensor(r32, force32);
      }
      
      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal,6, (double*)&wsum, std::plus<double>());
      w += wsum;*/
    }

    template < typename _AngularPotential > inline void
    FixedTripleListPIadressInteractionTemplate < _AngularPotential >::
    computeVirialTensor(Tensor& w, real z) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor of the triples");
      std::cout << "Warning! At the moment IK computeVirialTensor in FixedTripleListPIadressInteractionTemplate does not work." << std::endl;
      exit(1);
      return;
      
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
    FixedTripleListPIadressInteractionTemplate < _AngularPotential >::
    computeVirialTensor(Tensor *w, int n) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor of the triples");
      std::cout << "Warning! At the moment IK computeVirialTensor in FixedTripleListPIadressInteractionTemplate does not work." << std::endl;
      exit(1);
      return;
    }
    
    template < typename _AngularPotential >
    inline real
    FixedTripleListPIadressInteractionTemplate< _AngularPotential >::
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