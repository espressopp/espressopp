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
#ifndef _INTERACTION_VERLETLISTPIADRESSINTERACTIONTEMPLATE_HPP
#define _INTERACTION_VERLETLISTPIADRESSINTERACTIONTEMPLATE_HPP

//#include <typeinfo>

#include "System.hpp"
#include "bc/BC.hpp"

#include "types.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "VerletListAdress.hpp"
#include "FixedTupleListAdress.hpp"
#include "esutil/Array2D.hpp"
#include "SystemAccess.hpp"

namespace espressopp {
  namespace interaction {
    template < typename _PotentialQM, typename _PotentialCL >
    class VerletListPIadressInteractionTemplate: public Interaction {
    
    protected:
      typedef _PotentialQM PotentialQM;
      typedef _PotentialCL PotentialCL;
    
    public:
      VerletListPIadressInteractionTemplate
      (shared_ptr<VerletListAdress> _verletList, shared_ptr<FixedTupleListAdress> _fixedtupleList, int _ntrotter)
                : verletList(_verletList), fixedtupleList(_fixedtupleList), ntrotter(_ntrotter) {

          potentialArrayQM = esutil::Array2D<PotentialQM, esutil::enlarge>(0, 0, PotentialQM());
          potentialArrayCL = esutil::Array2D<PotentialCL, esutil::enlarge>(0, 0, PotentialCL());

          // AdResS stuff
          dhy = verletList->getHy();
          pidhy2 = M_PI/(dhy * 2.0);
          dex = verletList->getEx();
          dex2 = dex * dex;
          dexdhy = dex + verletList->getHy();
          dexdhy2 = dexdhy * dexdhy;
          
          ntypes = 0;
      }
                
      void
      setVerletList(shared_ptr < VerletListAdress > _verletList) {
        verletList = _verletList;
      }

      shared_ptr<VerletListAdress> getVerletList() {
        return verletList;
      }

      void
      setFixedTupleList(shared_ptr<FixedTupleListAdress> _fixedtupleList) {
          fixedtupleList = _fixedtupleList;
      }
      
      void
      setNTrotter(int _ntrotter) {
          ntrotter = _ntrotter;
      }

      void
      setPotentialQM(int type1, int type2, const PotentialQM &potential) {
          // typeX+1 because i<ntypes
          ntypes = std::max(ntypes, std::max(type1+1, type2+1));
        
          potentialArrayQM.at(type1, type2) = potential;
          if (type1 != type2) { // add potential in the other direction
             potentialArrayQM.at(type2, type1) = potential;
          }
      }

      void
      setPotentialCL(int type1, int type2, const PotentialCL &potential) {
          // typeX+1 because i<ntypes
          ntypes = std::max(ntypes, std::max(type1+1, type2+1));
        
         potentialArrayCL.at(type1, type2) = potential;
         if (type1 != type2) { // add potential in the other direction
             potentialArrayCL.at(type2, type1) = potential;
         }
      }

      PotentialQM &getPotentialQM(int type1, int type2) {
        return potentialArrayQM.at(type1, type2);
      }

      PotentialCL &getPotentialCL(int type1, int type2) {
        return potentialArrayCL.at(type1, type2);
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
      virtual int bondType() { return Nonbonded; }
      
    protected: 
      int ntypes;
      shared_ptr<VerletListAdress> verletList;
      shared_ptr<FixedTupleListAdress> fixedtupleList;
      esutil::Array2D<PotentialQM, esutil::enlarge> potentialArrayQM;
      esutil::Array2D<PotentialCL, esutil::enlarge> potentialArrayCL;

      // AdResS stuff
      real pidhy2; // pi / (dhy * 2)
      real dexdhy; // dex + dhy
      real dexdhy2; // dexdhy^2
      real dex;
      real dhy;
      real dex2; // dex^2
      int ntrotter; // Trotter number
      // Drift term WARNING - DOES NOT WORK CURRENTLY !!! Hence no need to set to set up...
      //std::map<Particle*, real> energydiff;  // Energydifference V_AA - V_CG map for particles in hybrid region for drift term calculation in H-AdResS
      std::set<Particle*> adrZone;  // Virtual particles in AdResS zone (HY and AT region)
      std::set<Particle*> cgZone;
      
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _PotentialQM, typename _PotentialCL > inline void
    VerletListPIadressInteractionTemplate < _PotentialQM, _PotentialCL >::
    addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by the Verlet List");
 
      // Get the cg and the adrZones
      std::set<Particle*> cgZone = verletList->getCGZone();
      std::set<Particle*> adrZone = verletList->getAdrZone();

      // Initialize the energy diff map to zero 
      // Drift term WARNING - DOES NOT WORK CURRENTLY !!! Hence no need to set to zero...
      /*for (std::set<Particle*>::iterator it=adrZone.begin();
                    it != adrZone.end(); ++it) {
                  	Particle &p = **it;
                  	// intitialize energy diff AA-CG
                  	energydiff[&p]=0.0;
      }*/
      

      // Pairs not inside the QM/Hybrid Zone (i.e. CL region)       

      // REMOVE FOR IDEAL GAS
      for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it) {

        // Get particles from pairlist
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
                  
        // Calculate forces in CL region
        const PotentialCL &potentialCL = getPotentialCL(p1.type(), p2.type());
        Real3D forcecl(0.0, 0.0, 0.0);
        if(potentialCL._computeForce(forcecl, p1, p2)) {     
             p1.force() += forcecl;
             p2.force() -= forcecl;  
        }
        
      }
      // REMOVE FOR IDEAL GAS
      
     
      // Pairs inside the QM/Hybrid Zone      
      for (PairList::Iterator it(verletList->getAdrPairs()); it.isValid(); ++it) {
         
         // for weights 
         real w1, w2;
         
         // Get particles from pairlist
         Particle &p1 = *it->first;
         Particle &p2 = *it->second;

         // Get their weights
         w1 = p1.lambda();               
         w2 = p2.lambda();        
         real w12 = (w1 + w2)/2.0;

         // Get the corresponding tuples
         FixedTupleListAdress::iterator it3;
         FixedTupleListAdress::iterator it4;
         it3 = fixedtupleList->find(&p1);
         it4 = fixedtupleList->find(&p2);

         if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

             // Get the PI bead lists (i.e. the AdResS particles)
             std::vector<Particle*> atList1;
             std::vector<Particle*> atList2;
             atList1 = it3->second;
             atList2 = it4->second;

             // Iterate the two iterators in a parallel fashion
             std::vector<Particle*>::iterator itv2 = atList2.begin();
             for (std::vector<Particle*>::iterator itv = atList1.begin();
                     itv != atList1.end(); ++itv) {

                 // they should be the same length... Total Trotter beads the same everywhere in the system
                 if (itv2 == atList2.end()){
                     std::cout << "Tuplelists seem to have different lengths or not started properly. Corresponding to CG particle " << 
                     p1.id() << "\n";
                     exit(1);
                     return;
                 }
                 
                 // Get the individual PI beads
                 Particle &p3 = **itv;
                 Particle &p4 = **itv2;
                 
                 // the beads we get should have the same Trotter bead number to interact with each other
                 if (p3.pib() != p4.pib()){
                     std::cout << "Path Integral Beads numbers do not correspond in VerletListPIadressInteractionTemplate for particles " << 
                     p3.id() << " and " << p4.id() << "\n";
                     exit(1);
                     return;
                 }
                 
                 // If at least one particle is in the hybrid region, calculate both forces.
                 if (w12 != 1.0) {
                 
                     // REMOVE FOR IDEAL GAS
                     // Calculate CL forces
                     const PotentialCL &potentialCL = getPotentialCL(p3.type(), p4.type());
                     Real3D forcecl(0.0, 0.0, 0.0);
                     if(potentialCL._computeForce(forcecl, p3, p4)) {
                         forcecl *= (1.0 - w12)/ntrotter;
                         p3.force() += forcecl;
                         p4.force() -= forcecl;  
                     }
                     // REMOVE FOR IDEAL GAS

                     // Calculate QM forces
                     const PotentialQM &potentialQM = getPotentialQM(p3.type(), p4.type());
                     Real3D forceqm(0.0, 0.0, 0.0);
                     if(potentialCL._computeForce(forceqm, p3, p4)) {
                         forceqm *= w12/ntrotter;
                         p3.force() += forceqm;
                         p4.force() -= forceqm;  
                     }
                     
                     // Drift term WARNING - DOES NOT WORK CURRENTLY !!!
                     /*if (w12 != 0.0) {   //at least one particle in hybrid region => need to do the energy calculation
                        real energyvpcl = potentialCL._computeEnergy(p3, p4);
                        if (w1 != 0.0) {   // if particle one is in hybrid region
                            energydiff[&p1] += energyvpcl;   // add CL energy for virtual particle 1
                        }
                        if (w2 != 0.0) {   // if particle two is in hybrid region
                            energydiff[&p2] += energyvpcl;   // add CL energy for virtual particle 2
                        }
                        
                        real energyvpqm = potentialQM._computeEnergy(p3, p4);
                        if (w1 != 0.0) {   // if particle one is in hybrid region
                            energydiff[&p1] -= energyvpqm;   // add CL energy for virtual particle 1
                        }
                        if (w2 != 0.0) {   // if particle two is in hybrid region
                            energydiff[&p2] -= energyvpqm;   // add CL energy for virtual particle 2
                        }
                        
                        // WARNING - DOES NOT WORK CURRENTLY !!!
                          std::cout << "No Drift terms in VerletListPIadressInteractionTemplate. This feature does currently not work.\n";
                          exit(1);
                          return;
                        // WARNING - DOES NOT WORK CURRENTLY !!!
                        
                    }*/
                     
                 }                 
                 else{   // both particles must be in the atomistic region now...

                     // Calculate QM forces
                     const PotentialQM &potentialQM = getPotentialQM(p3.type(), p4.type());
                     Real3D forceqm(0.0, 0.0, 0.0);
                     if(potentialQM._computeForce(forceqm, p3, p4)) {
                         forceqm *= 1.0/ntrotter;
                         p3.force() += forceqm;
                         p4.force() -= forceqm;  
                     }
                     
                 }
                 
                 //Iterate the second iterator
                 ++itv2;
                      
             }
             
         }
         else { // this should not happen
             std::cout << " one of the VP particles not found in tuples: " << p1.id() << "-" <<
                     p1.ghost() << ", " << p2.id() << "-" << p2.ghost();
             std::cout << " (" << p1.position() << ") (" << p2.position() << ")\n";
             exit(1);
             return;
         }
      
      }
      
      // Drift term WARNING - DOES NOT WORK CURRENTLY !!!
      // PI-AdResS - Drift Term application
      // Iterate over all particles in the hybrid region and calculate drift force
      
      /*for (std::set<Particle*>::iterator it=adrZone.begin();
        it != adrZone.end(); ++it) {   // Iterate over all particles
          Particle &vp = **it;
          real w = vp.lambda(); 
                  
          //if(w!=1.0 && w!=0.0){   //   only chose those in the hybrid region
          if(w<0.9999999 && w>0.0000001){   //   only chose those in the hybrid region
              // calculate distance to nearest adress particle or center
              std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
              Real3D pa = **it2; // position of adress particle
              Real3D mindriftforce(0.0, 0.0, 0.0);
              //Real3D mindriftforce = vp.position() - pa;                                                           // X SPLIT VS SPHERE CHANGE
              //real mindriftforce = vp.position()[0] - pa[0];                                                         // X SPLIT VS SPHERE CHANGE
              verletList->getSystem()->bc->getMinimumImageVector(mindriftforce, vp.position(), pa);
              real min1sq = mindriftforce[0]*mindriftforce[0]; // mindriftforce.sqr(); // set min1sq before loop
              ++it2;
              for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                   pa = **it2;
                   Real3D driftforce(0.0, 0.0, 0.0);
                   //Real3D driftforce = vp.position() - pa;                                                          // X SPLIT VS SPHERE CHANGE
                   //real driftforce = vp.position()[0] - pa[0];                                                    // X SPLIT VS SPHERE CHANGE
                   verletList->getSystem()->bc->getMinimumImageVector(driftforce, vp.position(), pa);
                   //real distsq1 = driftforce.sqr();//driftforce*driftforce; //driftforce.sqr();                     // X SPLIT VS SPHERE CHANGE
                   real distsq1 = driftforce[0]*driftforce[0];                                                          // X SPLIT VS SPHERE CHANGE
                   //std::cout << pa << " " << sqrt(distsq1) << "\n";
                   if (distsq1 < min1sq) {
                        min1sq = distsq1;
                        mindriftforce = driftforce;
                   }
              }
              min1sq = sqrt(min1sq);   // distance to nearest adress particle or center
              real mindriftforceX = (1.0/min1sq)*mindriftforce[0];  // normalized driftforce vector
              //mindriftforce *= weightderivative(min1sq);  // multiplication with derivative of the weighting function
              mindriftforceX *= 0.5;
              mindriftforceX *= energydiff.find(&vp)->second;   // get the energy differences which were calculated previously and put in drift force
                      
              vp.drift() += mindriftforceX ;//* vp.lambdaDeriv();    // USE ONLY LIKE THAT, IF DOING ITERATIVE FEC INCLUDING ITERATIVE PRESSURE FEC               
              
              mindriftforceX *= vp.lambdaDeriv();  // multiplication with derivative of the weighting function
              //vp.force() += mindriftforce;   // add drift force to virtual particles                                                                    // X SPLIT VS SPHERE CHANGE
              Real3D driftforceadd(mindriftforceX,0.0,0.0);                                                                                            // X SPLIT VS SPHERE CHANGE
              //Real3D driftforceadd(0.0,0.0,0.0);   
              vp.force() += driftforceadd;             // Goes in, if one wants to apply the "normal" drift force - also improve using [0] ...           // X SPLIT VS SPHERE CHANGE
              //std::cout << "Added Drift Force: " << driftforceadd << " for particle at pos(x).: " << vp.position()[0] << "\n";
              
          }
          
      }*/
      
      
      // Drift term WARNING - DOES NOT WORK CURRENTLY !!! Hence no need to clear any map...
      // energydiff.clear();  // clear the energy difference map                 
      
    }
    
    // Energy calculation does currently only work if integrator.run( ) (also with 0) and decompose have been executed before. This is due to the initialization of the tuples.
    template < typename _PotentialQM, typename _PotentialCL >
    inline real
    VerletListPIadressInteractionTemplate < _PotentialQM, _PotentialCL >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy of the Verlet list pairs");

      // Pairs not inside the QM/Hybrid Zone (i.e. CL region)
      
      // REMOVE FOR IDEAL GAS 
      real e = 0.0;        
      for (PairList::Iterator it(verletList->getPairs()); 
           it.isValid(); ++it) {
          Particle &p1 = *it->first;
          Particle &p2 = *it->second;
          int type1 = p1.type();
          int type2 = p2.type();
          const PotentialCL &potential = getPotentialCL(type1, type2);
          e += potential._computeEnergy(p1, p2);
      }
      // REMOVE FOR IDEAL GAS
           
      
      // Pairs inside the QM/Hybrid Zone
      for (PairList::Iterator it(verletList->getAdrPairs()); 
           it.isValid(); ++it) {
         Particle &p1 = *it->first;
         Particle &p2 = *it->second;
         
         // Get the weights
         real w1 = p1.lambda();
         real w2 = p2.lambda();
         real w12 = (w1 + w2)/2.0;
          
         // Get the corresponding tuples
         FixedTupleListAdress::iterator it3;
         FixedTupleListAdress::iterator it4;
         it3 = fixedtupleList->find(&p1);
         it4 = fixedtupleList->find(&p2);

         if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

             // Get the PI bead lists (i.e. the AdResS particles)
             std::vector<Particle*> atList1;
             std::vector<Particle*> atList2;
             atList1 = it3->second;
             atList2 = it4->second;

             // Iterate the two iterators in a parallel fashion
             std::vector<Particle*>::iterator itv2 = atList2.begin();
             for (std::vector<Particle*>::iterator itv = atList1.begin();
                     itv != atList1.end(); ++itv) {

                 // they should be the same length... Total Trotter beads the same everywhere in the system
                 if (itv2 == atList2.end()){
                     std::cout << "Tuplelists seem to have different lengths or not started properly. Corresponding to CG particle " << 
                     p1.id() << "\n";
                     exit(1);
                     return 0.0;
                 }
                 
                 // Get the individual PI beads
                 Particle &p3 = **itv;
                 Particle &p4 = **itv2;
                 
                 // the beads we get should have the same Trotter bead number to interact with each other
                 if (p3.pib()!= p4.pib()){
                     std::cout << "Path Integral Beads numbers do not correspond in VerletListPIadressInteractionTemplate for particles " << 
                     p3.id() << " and " << p4.id() << "\n";
                     exit(1);
                     return 0.0;
                 }
                 
                 // If at least one particle is in the hybrid region, calculate both forces.
                 if (w12 != 1.0) {
                 
                     // REMOVE FOR IDEAL GAS
                     // Calculate CL energy
                     const PotentialCL &potentialCL = getPotentialCL(p3.type(), p4.type());
                     e += ((1.0-w12)/ntrotter)*potentialCL._computeEnergy(p3, p4);
                     // REMOVE FOR IDEAL GAS
                     
                     // Calculate QM energy
                     const PotentialQM &potentialQM = getPotentialQM(p3.type(), p4.type());
                     e += (w12/ntrotter)*potentialQM._computeEnergy(p3, p4);
                     
                 }                 
                 else{   // both particles must be in the atomistic region now...

                     // Calculate QM energy
                     const PotentialQM &potential = getPotentialQM(p3.type(), p4.type());
                     e += (1.0/ntrotter)*potential._computeEnergy(p3, p4);
                     
                 }
                 
                 //Iterate the second iterator
                 ++itv2;
                      
             }
             
         }
         else { // this should not happen
             std::cout << " one of the VP particles not found in tuples: " << p1.id() << "-" <<
                     p1.ghost() << ", " << p2.id() << "-" << p2.ghost();
             std::cout << " (" << p1.position() << ") (" << p2.position() << ")\n";
             exit(1);
             return 0.0;
         }

      }
      real esum;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, e, esum, std::plus<real>());
      return esum;      
    }
    
    
    template < typename _PotentialQM, typename _PotentialCL >
    inline real
    VerletListPIadressInteractionTemplate < _PotentialQM, _PotentialCL >::
    computeEnergyAA() {
      LOG4ESPP_INFO(theLogger, "compute total AA energy of the Verlet list pairs");
      
      real e = 0.0;  

      for (PairList::Iterator it(verletList->getAdrPairs()); 
           it.isValid(); ++it) {
         Particle &p1 = *it->first;
         Particle &p2 = *it->second;
          
         // Get the corresponding tuples
         FixedTupleListAdress::iterator it3;
         FixedTupleListAdress::iterator it4;
         it3 = fixedtupleList->find(&p1);
         it4 = fixedtupleList->find(&p2);

         if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

             // Get the PI bead lists (i.e. the AdResS particles)
             std::vector<Particle*> atList1;
             std::vector<Particle*> atList2;
             atList1 = it3->second;
             atList2 = it4->second;

             // Iterate the two iterators in a parallel fashion
             std::vector<Particle*>::iterator itv2 = atList2.begin();
             for (std::vector<Particle*>::iterator itv = atList1.begin();
                     itv != atList1.end(); ++itv) {

                 // they should be the same length... Total Trotter beads the same everywhere in the system
                 if (itv2 == atList2.end()){
                     std::cout << "Tuplelists seem to have different lengths or not started properly. Corresponding to CG particle " << 
                     p1.id() << "\n";
                     exit(1);
                     return 0.0;
                 }
                 
                 // Get the individual PI beads
                 Particle &p3 = **itv;
                 Particle &p4 = **itv2;
                 
                 // the beads we get should have the same Trotter bead number to interact with each other
                 if (p3.pib() != p4.pib()){
                     std::cout << "Path Integral Beads numbers do not correspond in VerletListPIadressInteractionTemplate for particles " << 
                     p3.id() << " and " << p4.id() << "\n";
                     exit(1);
                     return 0.0;
                 }
                 

                 // Calculate QM energy
                 const PotentialQM &potential = getPotentialQM(p3.type(), p4.type());
                 e += (1.0/ntrotter)*potential._computeEnergy(p3, p4);                    
                 
                 //Iterate the second iterator
                 ++itv2;
                      
             }
             
         }
         else { // this should not happen
             std::cout << " one of the VP particles not found in tuples: " << p1.id() << "-" <<
                     p1.ghost() << ", " << p2.id() << "-" << p2.ghost();
             std::cout << " (" << p1.position() << ") (" << p2.position() << ")\n";
             exit(1);
             return 0.0;
         }

      }
      
      for (PairList::Iterator it(verletList->getPairs()); 
           it.isValid(); ++it) {
         Particle &p1 = *it->first;
         Particle &p2 = *it->second;
          
         // Get the corresponding tuples
         FixedTupleListAdress::iterator it3;
         FixedTupleListAdress::iterator it4;
         it3 = fixedtupleList->find(&p1);
         it4 = fixedtupleList->find(&p2);

         if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

             // Get the PI bead lists (i.e. the AdResS particles)
             std::vector<Particle*> atList1;
             std::vector<Particle*> atList2;
             atList1 = it3->second;
             atList2 = it4->second;

             // Iterate the two iterators in a parallel fashion
             std::vector<Particle*>::iterator itv2 = atList2.begin();
             for (std::vector<Particle*>::iterator itv = atList1.begin();
                     itv != atList1.end(); ++itv) {

                 // they should be the same length... Total Trotter beads the same everywhere in the system
                 if (itv2 == atList2.end()){
                     std::cout << "Tuplelists seem to have different lengths or not started properly. Corresponding to CG particle " << 
                     p1.id() << "\n";
                     exit(1);
                     return 0.0;
                 }
                 
                 // Get the individual PI beads
                 Particle &p3 = **itv;
                 Particle &p4 = **itv2;
                 
                 // the beads we get should have the same Trotter bead number to interact with each other
                 if (p3.pib() != p4.pib()){
                     std::cout << "Path Integral Beads numbers do not correspond in VerletListPIadressInteractionTemplate for particles " << 
                     p3.id() << " and " << p4.id() << "\n";
                     exit(1);
                     return 0.0;
                 }
                 

                 // Calculate QM energy
                 const PotentialQM &potential = getPotentialQM(p3.type(), p4.type());
                 e += (1.0/ntrotter)*potential._computeEnergy(p3, p4);                    
                 
                 //Iterate the second iterator
                 ++itv2;
                      
             }
             
         }
         else { // this should not happen
             std::cout << " one of the VP particles not found in tuples: " << p1.id() << "-" <<
                     p1.ghost() << ", " << p2.id() << "-" << p2.ghost();
             std::cout << " (" << p1.position() << ") (" << p2.position() << ")\n";
             exit(1);
             return 0.0;
         }

      }
       
      real esum;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, e, esum, std::plus<real>());
      return esum;      
    }
   
    
    template < typename _PotentialQM, typename _PotentialCL >
    inline real
    VerletListPIadressInteractionTemplate < _PotentialQM, _PotentialCL >::
    computeEnergyCG() {
      LOG4ESPP_INFO(theLogger, "compute total CG energy of the Verlet list pairs");
       
      real e = 0.0;  

      for (PairList::Iterator it(verletList->getAdrPairs()); 
           it.isValid(); ++it) {
         Particle &p1 = *it->first;
         Particle &p2 = *it->second;
          
         // Get the corresponding tuples
         FixedTupleListAdress::iterator it3;
         FixedTupleListAdress::iterator it4;
         it3 = fixedtupleList->find(&p1);
         it4 = fixedtupleList->find(&p2);

         if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

             // Get the PI bead lists (i.e. the AdResS particles)
             std::vector<Particle*> atList1;
             std::vector<Particle*> atList2;
             atList1 = it3->second;
             atList2 = it4->second;

             // Iterate the two iterators in a parallel fashion
             std::vector<Particle*>::iterator itv2 = atList2.begin();
             for (std::vector<Particle*>::iterator itv = atList1.begin();
                     itv != atList1.end(); ++itv) {

                 // they should be the same length... Total Trotter beads the same everywhere in the system
                 if (itv2 == atList2.end()){
                     std::cout << "Tuplelists seem to have different lengths or not started properly. Corresponding to CG particle " << 
                     p1.id() << "\n";
                     exit(1);
                     return 0.0;
                 }
                 
                 // Get the individual PI beads
                 Particle &p3 = **itv;
                 Particle &p4 = **itv2;
                 
                 // the beads we get should have the same Trotter bead number to interact with each other
                 if (p3.pib() != p4.pib()){
                     std::cout << "Path Integral Beads numbers do not correspond in VerletListPIadressInteractionTemplate for particles " << 
                     p3.id() << " and " << p4.id() << "\n";
                     exit(1);
                     return 0.0;
                 }
                 

                 // Calculate QM energy
                 const PotentialCL &potential = getPotentialCL(p3.type(), p4.type());
                 e += (1.0/ntrotter)*potential._computeEnergy(p3, p4);                    
                 
                 //Iterate the second iterator
                 ++itv2;
                      
             }
             
         }
         else { // this should not happen
             std::cout << " one of the VP particles not found in tuples: " << p1.id() << "-" <<
                     p1.ghost() << ", " << p2.id() << "-" << p2.ghost();
             std::cout << " (" << p1.position() << ") (" << p2.position() << ")\n";
             exit(1);
             return 0.0;
         }

      }
      
      for (PairList::Iterator it(verletList->getPairs()); 
           it.isValid(); ++it) {
         Particle &p1 = *it->first;
         Particle &p2 = *it->second;
          
         // Get the corresponding tuples
         FixedTupleListAdress::iterator it3;
         FixedTupleListAdress::iterator it4;
         it3 = fixedtupleList->find(&p1);
         it4 = fixedtupleList->find(&p2);

         if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

             // Get the PI bead lists (i.e. the AdResS particles)
             std::vector<Particle*> atList1;
             std::vector<Particle*> atList2;
             atList1 = it3->second;
             atList2 = it4->second;

             // Iterate the two iterators in a parallel fashion
             std::vector<Particle*>::iterator itv2 = atList2.begin();
             for (std::vector<Particle*>::iterator itv = atList1.begin();
                     itv != atList1.end(); ++itv) {

                 // they should be the same length... Total Trotter beads the same everywhere in the system
                 if (itv2 == atList2.end()){
                     std::cout << "Tuplelists seem to have different lengths or not started properly. Corresponding to CG particle " << 
                     p1.id() << "\n";
                     exit(1);
                     return 0.0;
                 }
                 
                 // Get the individual PI beads
                 Particle &p3 = **itv;
                 Particle &p4 = **itv2;
                 
                 // the beads we get should have the same Trotter bead number to interact with each other
                 if (p3.pib() != p4.pib()){
                     std::cout << "Path Integral Beads numbers do not correspond in VerletListPIadressInteractionTemplate for particles " << 
                     p3.id() << " and " << p4.id() << "\n";
                     exit(1);
                     return 0.0;
                 }
                 

                 // Calculate QM energy
                 const PotentialCL &potential = getPotentialCL(p3.type(), p4.type());
                 e += (1.0/ntrotter)*potential._computeEnergy(p3, p4);                    
                 
                 //Iterate the second iterator
                 ++itv2;
                      
             }
             
         }
         else { // this should not happen
             std::cout << " one of the VP particles not found in tuples: " << p1.id() << "-" <<
                     p1.ghost() << ", " << p2.id() << "-" << p2.ghost();
             std::cout << " (" << p1.position() << ") (" << p2.position() << ")\n";
             exit(1);
             return 0.0;
         }

      }
       
      real esum;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, e, esum, std::plus<real>());
      return esum;      
    }
    

    template < typename _PotentialQM, typename _PotentialCL > inline void
    VerletListPIadressInteractionTemplate < _PotentialQM, _PotentialCL >::
    computeVirialX(std::vector<real> &p_xx_total, int bins) {
      LOG4ESPP_INFO(theLogger, "compute virial p_xx of the pressure tensor slabwise");
      
      std::cout << "Warning! At the moment computeVirialX in VerletListPIadress does not work"<<std::endl;
      exit(1);
      return;
      
      /*std::set<Particle*> cgZone = verletList->getCGZone();
      for (std::set<Particle*>::iterator it=cgZone.begin();
              it != cgZone.end(); ++it) {

          Particle &vp = **it;

          FixedTupleListAdress::iterator it3;
          it3 = fixedtupleList->find(&vp);

          if (it3 != fixedtupleList->end()) {

              std::vector<Particle*> atList;
              atList = it3->second;

              // compute center of mass
              Real3D cmp(0.0, 0.0, 0.0); // center of mass position
              //Real3D cmv(0.0, 0.0, 0.0); // center of mass velocity
              //real M = vp.getMass(); // sum of mass of AT particles
              for (std::vector<Particle*>::iterator it2 = atList.begin();
                                   it2 != atList.end(); ++it2) {
                  Particle *at = *it2;
                  //Real3D d1 = at.position() - vp.position();
                  //Real3D d1;
                  //verletList->getSystem()->bc->getMinimumImageVectorBox(d1, at.position(), vp.position());
                  //cmp += at.mass() * d1;

                  cmp += at->mass() * at->position();
                  //cmv += at.mass() * at.velocity();
              }
              cmp /= vp.getMass();
              //cmv /= vp.getMass();
              //cmp += vp.position(); // cmp is a relative position
              //std::cout << " cmp M: "  << M << "\n\n";
              //std::cout << "  moving VP to " << cmp << ", velocitiy is " << cmv << "\n";

              // update (overwrite) the position and velocity of the VP
              vp.position() = cmp;
              //vp.velocity() = cmv;

          }
          else { // this should not happen
              std::cout << " VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
              std::cout << " (" << vp.position() << ")\n";
              exit(1);
              return;
          }
      }
      
      std::set<Particle*> adrZone = verletList->getAdrZone();
      for (std::set<Particle*>::iterator it=adrZone.begin();
              it != adrZone.end(); ++it) {

          Particle &vp = **it;

          FixedTupleListAdress::iterator it3;
          it3 = fixedtupleList->find(&vp);

          if (it3 != fixedtupleList->end()) {

              std::vector<Particle*> atList;
              atList = it3->second;

              // compute center of mass
              Real3D cmp(0.0, 0.0, 0.0); // center of mass position
              //Real3D cmv(0.0, 0.0, 0.0); // center of mass velocity
              //real M = vp.getMass(); // sum of mass of AT particles
              for (std::vector<Particle*>::iterator it2 = atList.begin();
                                   it2 != atList.end(); ++it2) {
                  Particle &at = **it2;
                  //Real3D d1 = at.position() - vp.position();
                  //Real3D d1;
                  //verletList->getSystem()->bc->getMinimumImageVectorBox(d1, at.position(), vp.position());
                  //cmp += at.mass() * d1;

                  cmp += at.mass() * at.position();
                  //cmv += at.mass() * at.velocity();
              }
              cmp /= vp.getMass();
              //cmv /= vp.getMass();
              //cmp += vp.position(); // cmp is a relative position
              //std::cout << " cmp M: "  << M << "\n\n";
              //std::cout << "  moving VP to " << cmp << ", velocitiy is " << cmv << "\n";

              // update (overwrite) the position and velocity of the VP
              vp.position() = cmp;
              //vp.velocity() = cmv;

          }
          else { // this should not happen
              std::cout << " VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
              std::cout << " (" << vp.position() << ")\n";
              exit(1);
              return;
          }
      }
      
      int i = 0;
      int bin1 = 0;
      int bin2 = 0;
      
      System& system = verletList->getSystemRef();
      Real3D Li = system.bc->getBoxL();
      real Delta_x = Li[0] / (real)bins;
      real Volume = Li[1] * Li[2] * Delta_x;
      size_t size = bins;
      std::vector <real> p_xx_local(size);      

      for (i = 0; i < bins; ++i)
        {
          p_xx_local[i] = 0.0;
        }
      
      for (PairList::Iterator it(verletList->getPairs());                
           it.isValid(); ++it) {                                         
        Particle &p1 = *it->first;                                       
        Particle &p2 = *it->second;                                      
        int type1 = p1.type();                                           
        int type2 = p2.type();
        const PotentialCG &potential = getPotentialCG(type1, type2);

        Real3D force(0.0, 0.0, 0.0);
        if(potential._computeForce(force, p1, p2)) {
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
      
      //if (KTI == false) {
      //makeWeights();
      //}
      
      for (PairList::Iterator it(verletList->getAdrPairs()); it.isValid(); ++it) {
         real w1, w2;
         // these are the two VP interacting
         Particle &p1 = *it->first;
         Particle &p2 = *it->second;
     
         Real3D dist = p1.position() - p2.position();
         w1 = p1.lambda();         
         w2 = p2.lambda();
         real w12 = (w1 + w2)/2.0;  // H-AdResS

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
         
         // force between VP particles
         int type1 = p1.type();
         int type2 = p2.type();
         const PotentialCG &potentialCG = getPotentialCG(type1, type2);
         Real3D forcevp(0.0, 0.0, 0.0);
                if (w12 != 1.0) { // calculate VP force if both VP are outside AT region (CG-HY, HY-HY)
                    if (potentialCG._computeForce(forcevp, p1, p2)) {
                          forcevp *= (1.0 - w12);
                          real vir_temp = 0.5 * dist[0] * forcevp[0];                          
                          p_xx_local.at(bin1) += vir_temp;
                          p_xx_local.at(bin2) += vir_temp;  

                    }
                }
       
         // force between AT particles
         if (w12 != 0.0) { // calculate AT force if both VP are outside CG region (HY-HY, HY-AT, AT-AT)
             FixedTupleListAdress::iterator it3;
             FixedTupleListAdress::iterator it4;
             it3 = fixedtupleList->find(&p1);
             it4 = fixedtupleList->find(&p2);

             //std::cout << "Interaction " << p1.id() << " - " << p2.id() << "\n";
             if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

                 std::vector<Particle*> atList1;
                 std::vector<Particle*> atList2;
                 atList1 = it3->second;
                 atList2 = it4->second;
                 
                 Real3D force_temp(0.0, 0.0, 0.0);
                 
                 for (std::vector<Particle*>::iterator itv = atList1.begin();
                         itv != atList1.end(); ++itv) {

                     Particle &p3 = **itv;

                     for (std::vector<Particle*>::iterator itv2 = atList2.begin();
                                          itv2 != atList2.end(); ++itv2) {

                         Particle &p4 = **itv2;

                         // AT forces
                         const PotentialAT &potentialAT = getPotentialAT(p3.type(), p4.type());
                         Real3D force(0.0, 0.0, 0.0);
                         if(potentialAT._computeForce(force, p3, p4)) {
                             force_temp += force; 
                         }
                     }
                 }
                 
                 force_temp *= w12;
                 real vir_temp = 0.5 * dist[0] * force_temp[0];                
                 p_xx_local.at(bin1) += vir_temp;
                 p_xx_local.at(bin2) += vir_temp;                                                                          
             }
             else { // this should not happen
                 std::cout << " one of the VP particles not found in tuples: " << p1.id() << "-" <<
                         p1.ghost() << ", " << p2.id() << "-" << p2.ghost();
                 std::cout << " (" << p1.position() << ") (" << p2.position() << ")\n";
                 exit(1);
                 return;
             }
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
          p_xx_total.at(i) += p_xx_sum.at(i);
      }*/
      
    }



    template < typename _PotentialQM, typename _PotentialCL > inline real
    VerletListPIadressInteractionTemplate < _PotentialQM, _PotentialCL >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute the virial for the Verlet List");
      
      std::cout << "Warning! At the moment IK computeVirial in VerletListPIadress does not work"<<std::endl;
      exit(1);
      return 0.0;
      
      /*std::set<Particle*> cgZone = verletList->getCGZone();
      for (std::set<Particle*>::iterator it=cgZone.begin();
              it != cgZone.end(); ++it) {

          Particle &vp = **it;

          FixedTupleListAdress::iterator it3;
          it3 = fixedtupleList->find(&vp);

          if (it3 != fixedtupleList->end()) {

              std::vector<Particle*> atList;
              atList = it3->second;

              // compute center of mass
              Real3D cmp(0.0, 0.0, 0.0); // center of mass position
              //Real3D cmv(0.0, 0.0, 0.0); // center of mass velocity
              //real M = vp.getMass(); // sum of mass of AT particles
              for (std::vector<Particle*>::iterator it2 = atList.begin();
                                   it2 != atList.end(); ++it2) {
                  Particle *at = *it2;
                  //Real3D d1 = at.position() - vp.position();
                  //Real3D d1;
                  //verletList->getSystem()->bc->getMinimumImageVectorBox(d1, at.position(), vp.position());
                  //cmp += at.mass() * d1;

                  cmp += at->mass() * at->position();
                  //cmv += at.mass() * at.velocity();
              }
              cmp /= vp.getMass();
              //cmv /= vp.getMass();
              //cmp += vp.position(); // cmp is a relative position
              //std::cout << " cmp M: "  << M << "\n\n";
              //std::cout << "  moving VP to " << cmp << ", velocitiy is " << cmv << "\n";

              // update (overwrite) the position and velocity of the VP
              vp.position() = cmp;
              //vp.velocity() = cmv;

          }
          else { // this should not happen
              std::cout << " VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
              std::cout << " (" << vp.position() << ")\n";
              exit(1);
              return 0.0;
          }
      }
      
      std::set<Particle*> adrZone = verletList->getAdrZone();
      for (std::set<Particle*>::iterator it=adrZone.begin();
              it != adrZone.end(); ++it) {

          Particle &vp = **it;

          FixedTupleListAdress::iterator it3;
          it3 = fixedtupleList->find(&vp);

          if (it3 != fixedtupleList->end()) {

              std::vector<Particle*> atList;
              atList = it3->second;

              // compute center of mass
              Real3D cmp(0.0, 0.0, 0.0); // center of mass position
              //Real3D cmv(0.0, 0.0, 0.0); // center of mass velocity
              //real M = vp.getMass(); // sum of mass of AT particles
              for (std::vector<Particle*>::iterator it2 = atList.begin();
                                   it2 != atList.end(); ++it2) {
                  Particle &at = **it2;
                  //Real3D d1 = at.position() - vp.position();
                  //Real3D d1;
                  //verletList->getSystem()->bc->getMinimumImageVectorBox(d1, at.position(), vp.position());
                  //cmp += at.mass() * d1;

                  cmp += at.mass() * at.position();
                  //cmv += at.mass() * at.velocity();
              }
              cmp /= vp.getMass();
              //cmv /= vp.getMass();
              //cmp += vp.position(); // cmp is a relative position
              //std::cout << " cmp M: "  << M << "\n\n";
              //std::cout << "  moving VP to " << cmp << ", velocitiy is " << cmv << "\n";

              // update (overwrite) the position and velocity of the VP
              vp.position() = cmp;
              //vp.velocity() = cmv;

          }
          else { // this should not happen
              std::cout << " VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
              std::cout << " (" << vp.position() << ")\n";
              exit(1);
              return 0.0;
          }
      }*/
      //const bc::BC& bc = *getSystemRef().bc;
      
      /*real w = 0.0; 
 
      for (PairList::Iterator it(verletList->getPairs());                
           it.isValid(); ++it) {*/
        /*std::cout << "Should not happen!\n";
        exit(1);
        return 0.0;*/
        /*Particle &p1 = *it->first;                                       
        Particle &p2 = *it->second;                                      
        int type1 = p1.type();                                           
        int type2 = p2.type();
        const PotentialCG &potential = getPotentialCG(type1, type2);

        Real3D force(0.0, 0.0, 0.0);
        if(potential._computeForce(force, p1, p2)) {
          Real3D dist = p1.position() - p2.position();
          //Real3D dist;
          //bc.getMinimumImageVectorBox(dist, p1.position(), p2.position());
          
          w += dist * force;
        }
      }
      
      //if (KTI == false) {
      //makeWeights();
      //}
      
      for (PairList::Iterator it(verletList->getAdrPairs()); it.isValid(); ++it) {
         real w1, w2;
         // these are the two VP interacting
         Particle &p1 = *it->first;
         Particle &p2 = *it->second;     

         w1 = p1.lambda();         
         w2 = p2.lambda();
         real w12 = (w1 + w2)/2.0;  // H-AdResS
         
         // force between VP particles
         int type1 = p1.type();
         int type2 = p2.type();
         const PotentialCG &potentialCG = getPotentialCG(type1, type2);
         Real3D forcevp(0.0, 0.0, 0.0);
                if (w12 != 1.0) {*/ // calculate VP force if both VP are outside AT region (CG-HY, HY-HY)
                    /*std::cout << "Should not happen!\n";
                    exit(1);
                    return 0.0;*/
                    /*if (potentialCG._computeForce(forcevp, p1, p2)) {
                          forcevp *= (1.0 - w12);
                          Real3D dist = p1.position() - p2.position();
                          //Real3D dist;
                          //bc.getMinimumImageVectorBox(dist, p1.position(), p2.position());
                          w += dist * forcevp;  

                    }
                }
       
         // force between AT particles
         if (w12 != 0.0) { // calculate AT force if both VP are outside CG region (HY-HY, HY-AT, AT-AT)
             FixedTupleListAdress::iterator it3;
             FixedTupleListAdress::iterator it4;
             it3 = fixedtupleList->find(&p1);
             it4 = fixedtupleList->find(&p2);

             //std::cout << "Interaction " << p1.id() << " - " << p2.id() << "\n";
             if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

                 std::vector<Particle*> atList1;
                 std::vector<Particle*> atList2;
                 atList1 = it3->second;
                 atList2 = it4->second;                
                 
                 for (std::vector<Particle*>::iterator itv = atList1.begin();
                         itv != atList1.end(); ++itv) {

                     Particle &p3 = **itv;

                     for (std::vector<Particle*>::iterator itv2 = atList2.begin();
                                          itv2 != atList2.end(); ++itv2) {

                         Particle &p4 = **itv2;

                         // AT forces
                         const PotentialAT &potentialAT = getPotentialAT(p3.type(), p4.type());
                         Real3D forceat(0.0, 0.0, 0.0);
                         if(potentialAT._computeForce(forceat, p3, p4)) {
                         forceat *= w12;
                         Real3D dist = p3.position() - p4.position();
                         //Real3D dist;
                         //bc.getMinimumImageVectorBox(dist, p3.position(), p3.position());*/
                         /*if (dist * forceat > 400.0){
                                std::cout << "IDs: id_p3=" << p3.id() << " and id_p4=" << p4.id() << "   #####   ";                             
                                std::cout << "Types: type_p3=" << p3.type() << " and type_p4=" << p4.type() << "   #####   ";   
                                std::cout << "dist: " << sqrt(dist*dist) << "   #####   ";
                                std::cout << "force: " << forceat << "   #####   ";
                                std::cout << "product: " << dist * forceat << "\n";
                                if (sqrt(dist*dist) > 1.0){
                                        std::cout << "DIST TOO BIG!\n";
                                        exit(1);
                                        return 0.0;
                                }
                         }*/
                         
                         /*w += dist * forceat; 
                         }
                     }
                 }
                 
             }
             else { // this should not happen
                 std::cout << " one of the VP particles not found in tuples: " << p1.id() << "-" <<
                         p1.ghost() << ", " << p2.id() << "-" << p2.ghost();
                 std::cout << " (" << p1.position() << ") (" << p2.position() << ")\n";
                 exit(1);
                 return 0.0;
             }
         }*/
      /*}      
      
      real wsum;
      wsum = 0.0;
      boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
      return wsum;      */
      
      /*real w = 0.0;           OLD STUFF
      for (PairList::Iterator it(verletList->getPairs());                
           it.isValid(); ++it) {                                         
        Particle &p1 = *it->first;                                       
        Particle &p2 = *it->second;                                      
        int type1 = p1.type();                                           
        int type2 = p2.type();
        const PotentialCG &potential = getPotentialCG(type1, type2);

        Real3D force(0.0, 0.0, 0.0);
        if(potential._computeForce(force, p1, p2)) {
          Real3D dist = p1.position() - p2.position();
          w = w + dist * force;
        }
      }

      for (PairList::Iterator it(verletList->getAdrPairs());
                 it.isValid(); ++it) {
              Particle &p1 = *it->first;
              Particle &p2 = *it->second;
              int type1 = p1.type();
              int type2 = p2.type();
              const PotentialCG &potential = getPotentialCG(type1, type2);

              Real3D force(0.0, 0.0, 0.0);
              if(potential._computeForce(force, p1, p2)) {
                Real3D dist = p1.position() - p2.position();
                w = w + dist * force;
              }
      }

      real wsum;
      boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
      return wsum;*/
    }

    template < typename _PotentialQM, typename _PotentialCL > inline void
    VerletListPIadressInteractionTemplate < _PotentialQM, _PotentialCL >::
    computeVirialTensor(Tensor& w) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the Verlet List");

      std::cout << "Warning! At the moment IK computeVirialTensor in VerletListPIadress does not work"<<std::endl;
      exit(1);
      return;
      
      /*Tensor wlocal(0.0);
      for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it){
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.type();
        int type2 = p2.type();
        const PotentialCG &potential = getPotentialCG(type1, type2);

        Real3D force(0.0, 0.0, 0.0);
        if(potential._computeForce(force, p1, p2)) {
          Real3D r21 = p1.position() - p2.position();
          wlocal += Tensor(r21, force);
        }
      }

      for (PairList::Iterator it(verletList->getAdrPairs()); it.isValid(); ++it){
              Particle &p1 = *it->first;
              Particle &p2 = *it->second;
              int type1 = p1.type();
              int type2 = p2.type();
              const PotentialCG &potential = getPotentialCG(type1, type2);

              Real3D force(0.0, 0.0, 0.0);
              if(potential._computeForce(force, p1, p2)) {
                Real3D r21 = p1.position() - p2.position();
                wlocal += Tensor(r21, force);
              }
      }

      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
      w += wsum;*/
    }
 
    template < typename _PotentialQM, typename _PotentialCL > inline void
    VerletListPIadressInteractionTemplate < _PotentialQM, _PotentialCL >::
    computeVirialTensor(Tensor& w, real z) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the Verlet List");

      std::cout << "Warning! At the moment IK computeVirialTensor in VerletListPIadress does not work"<<std::endl;
      exit(1);
      return;
      /*
      Tensor wlocal(0.0);
      for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it){
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.type();
        int type2 = p2.type();
        Real3D p1pos = p1.position();
        Real3D p2pos = p2.position();
        
        if(  (p1pos[0]>xmin && p1pos[0]<xmax && 
              p1pos[1]>ymin && p1pos[1]<ymax && 
              p1pos[2]>zmin && p1pos[2]<zmax) ||
             (p2pos[0]>xmin && p2pos[0]<xmax && 
              p2pos[1]>ymin && p2pos[1]<ymax && 
              p2pos[2]>zmin && p2pos[2]<zmax) ){
          const PotentialCG &potential = getPotentialCG(type1, type2);

          Real3D force(0.0, 0.0, 0.0);
          if(potential._computeForce(force, p1, p2)) {
            Real3D r21 = p1pos - p2pos;
            wlocal += Tensor(r21, force);
          }
        }
      }

      for (PairList::Iterator it(verletList->getAdrPairs()); it.isValid(); ++it){
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.type();
        int type2 = p2.type();
        Real3D p1pos = p1.position();
        Real3D p2pos = p2.position();
        if(  (p1pos[0]>xmin && p1pos[0]<xmax && 
              p1pos[1]>ymin && p1pos[1]<ymax && 
              p1pos[2]>zmin && p1pos[2]<zmax) ||
             (p2pos[0]>xmin && p2pos[0]<xmax && 
              p2pos[1]>ymin && p2pos[1]<ymax && 
              p2pos[2]>zmin && p2pos[2]<zmax) ){
          const PotentialCG &potential = getPotentialCG(type1, type2);

          Real3D force(0.0, 0.0, 0.0);
          if(potential._computeForce(force, p1, p2)) {
            Real3D r21 = p1pos - p2pos;
            wlocal += Tensor(r21, force);
          }
        }
      }

      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
      w += wsum;
       */
    }

    
    template < typename _PotentialQM, typename _PotentialCL > inline void
    VerletListPIadressInteractionTemplate < _PotentialQM, _PotentialCL >::
    computeVirialTensor(Tensor *w, int n) {
      std::cout << "Warning! At the moment IK computeVirialTensor in VerletListPIadress does not work"<<std::endl;
      exit(1);
      return;
    }
    
    template < typename _PotentialQM, typename _PotentialCL >
    inline real
    VerletListPIadressInteractionTemplate< _PotentialQM, _PotentialCL >::
    getMaxCutoff() {
      real cutoff = 0.0;
      for (int i = 0; i < ntypes; i++) {
        for (int j = 0; j < ntypes; j++) {
          cutoff = std::max(cutoff, getPotentialCL(i, j).getCutoff());
        }
      }
      return cutoff;
    }
  }
}
#endif
