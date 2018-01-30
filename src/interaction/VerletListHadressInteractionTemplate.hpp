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
#ifndef _INTERACTION_VERLETLISTHADRESSINTERACTIONTEMPLATE_HPP
#define _INTERACTION_VERLETLISTHADRESSINTERACTIONTEMPLATE_HPP

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
    template < typename _PotentialAT, typename _PotentialCG >
    class VerletListHadressInteractionTemplate: public Interaction {

    protected:
      typedef _PotentialAT PotentialAT;
      typedef _PotentialCG PotentialCG;

    public:
      VerletListHadressInteractionTemplate
      (shared_ptr<VerletListAdress> _verletList, shared_ptr<FixedTupleListAdress> _fixedtupleList)
                : verletList(_verletList), fixedtupleList(_fixedtupleList) {

          potentialArrayAT = esutil::Array2D<PotentialAT, esutil::enlarge>(0, 0, PotentialAT());
          potentialArrayCG = esutil::Array2D<PotentialCG, esutil::enlarge>(0, 0, PotentialCG());

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
      setPotentialAT(int type1, int type2, const PotentialAT &potential) {
          // typeX+1 because i<ntypes
          ntypes = std::max(ntypes, std::max(type1+1, type2+1));

          potentialArrayAT.at(type1, type2) = potential;
          if (type1 != type2) { // add potential in the other direction
             potentialArrayAT.at(type2, type1) = potential;
          }
      }

      void
      setPotentialCG(int type1, int type2, const PotentialCG &potential) {
          // typeX+1 because i<ntypes
          ntypes = std::max(ntypes, std::max(type1+1, type2+1));

         potentialArrayCG.at(type1, type2) = potential;
         if (type1 != type2) { // add potential in the other direction
             potentialArrayCG.at(type2, type1) = potential;
         }
      }

      PotentialAT &getPotentialAT(int type1, int type2) {
        return potentialArrayAT.at(type1, type2);
      }

      PotentialCG &getPotentialCG(int type1, int type2) {
        return potentialArrayCG.at(type1, type2);
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
      virtual int bondType() { return Nonbonded; }

    protected:
      int ntypes;
      shared_ptr<VerletListAdress> verletList;
      shared_ptr<FixedTupleListAdress> fixedtupleList;
      esutil::Array2D<PotentialAT, esutil::enlarge> potentialArrayAT;
      esutil::Array2D<PotentialCG, esutil::enlarge> potentialArrayCG;

      // AdResS stuff
      real pidhy2; // pi / (dhy * 2)
      real dexdhy; // dex + dhy
      real dexdhy2; // dexdhy^2
      real dex;
      real dhy;
      real dex2; // dex^2
      std::map<Particle*, real> energydiff;  // Energydifference V_AA - V_CG map for particles in hybrid region for drift term calculation in H-AdResS
      std::set<Particle*> adrZone;  // Virtual particles in AdResS zone (HY and AT region)
      std::set<Particle*> cgZone;

    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _PotentialAT, typename _PotentialCG > inline void
    VerletListHadressInteractionTemplate < _PotentialAT, _PotentialCG >::
    addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by the Verlet List");

      std::set<Particle*> cgZone = verletList->getCGZone();
      std::set<Particle*> adrZone = verletList->getAdrZone();

      for (std::set<Particle*>::iterator it=adrZone.begin();
                    it != adrZone.end(); ++it) {
                  	Particle &p = **it;
                  	// intitialize energy diff AA-CG
                  	energydiff[&p]=0.0;
      }


      // Pairs not inside the AdResS Zone (CG region)
      // COMMENT FOR IDEAL GAS
      for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it) {

        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.type();
        int type2 = p2.type();

        const PotentialCG &potentialCG = getPotentialCG(type1, type2);

        Real3D force(0.0, 0.0, 0.0);

        // CG forces
        if(potentialCG._computeForce(force, p1, p2)) {
          p1.force() += force;
          p2.force() -= force;
        }
      }
      // REMOVE FOR IDEAL GAS

      // Compute forces (AT and VP) of Pairs inside AdResS zone
      for (PairList::Iterator it(verletList->getAdrPairs()); it.isValid(); ++it) {
         real w1, w2;
         // these are the two VP interacting
         Particle &p1 = *it->first;
         Particle &p2 = *it->second;

         w1 = p1.lambda();
         w2 = p2.lambda();

         real w12 = (w1 + w2)/2.0;  // H-AdResS



         // REMOVE FOR IDEAL GAS
         // force between VP particles
         int type1 = p1.type();
         int type2 = p2.type();
         const PotentialCG &potentialCG = getPotentialCG(type1, type2);
         Real3D forcevp(0.0, 0.0, 0.0);
                if (w12 != 1.0) { // calculate VP force if both VP are outside AT region (CG-HY, HY-HY)
                    if (potentialCG._computeForce(forcevp, p1, p2)) {
                        forcevp *= (1.0 - w12);
                        p1.force() += forcevp;
                        p2.force() -= forcevp;
                    }

                    // H-AdResS - Drift Term part 1
                    // Compute CG energies of particles in the hybrid and store and add up in map energydiff
                    if (w12 != 0.0) {   //at least one particle in hybrid region => need to do the energy calculation
                        real energyvp = potentialCG._computeEnergy(p1, p2);
                        if (w1 != 0.0) {   // if particle one is in hybrid region
                            energydiff[&p1] += energyvp;   // add CG energy for virtual particle 1
                        }
                        if (w2 != 0.0) {   // if particle two is in hybrid region
                            energydiff[&p2] += energyvp;   // add CG energy for virtual particle 2
                        }
                    }

                }
         // REMOVE FOR IDEAL GAS

         // force between AT particles
         if (w12 != 0.0) { // calculate AT force if both VP are outside CG region (HY-HY, HY-AT, AT-AT)
             FixedTupleListAdress::iterator it3;
             FixedTupleListAdress::iterator it4;
             it3 = fixedtupleList->find(&p1);
             it4 = fixedtupleList->find(&p2);

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
                         Real3D force(0.0, 0.0, 0.0);
                         if(potentialAT._computeForce(force, p3, p4)) {
                             force *= w12;
                             p3.force() += force;
                             p4.force() -= force;
                         }

                         // H-AdResS - Drift Term part 2
                         // Compute AT energies of particles in the hybrid and store and subtract in map energydiff
                         if(w12!=1.0){   //at least one particle in hybrid region => need to do the energy calculation
                             real energyat = potentialAT._computeEnergy(p3, p4);
                             if(w1!=1.0){   // if particle one is in hybrid region
                                    energydiff[&p1] -= energyat;   // subtract AT energy for virtual particle 1
                             }
                             if(w2!=1.0){   // if particle two is in hybrid region
                                    energydiff[&p2] -= energyat;   // subtract AT energy for virtual particle 2
                             }
                         }

                     }

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
      }

      // H-AdResS - Drift Term part 3
      // Iterate over all particles in the hybrid region and calculate drift force
      for (std::set<Particle*>::iterator it=adrZone.begin();
        it != adrZone.end(); ++it) {   // Iterate over all particles
          Particle &vp = **it;
          real w = vp.lambda();

          if(w<0.9999999 && w>0.0000001){   //   only chose those in the hybrid region
              // calculate distance to nearest adress particle or center
              std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
              Real3D pa = **it2; // position of adress particle
              Real3D mindriftforce(0.0, 0.0, 0.0);
              verletList->getSystem()->bc->getMinimumImageVector(mindriftforce, vp.position(), pa);
              real min1sq = 0.0;
              if(verletList->getAdrRegionType()){
                min1sq = mindriftforce.sqr();
              }
              else{
                min1sq = mindriftforce[0]*mindriftforce[0];
              }
              ++it2;
              for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                   pa = **it2;
                   Real3D driftforce(0.0, 0.0, 0.0);
                   verletList->getSystem()->bc->getMinimumImageVector(driftforce, vp.position(), pa);
                   real distsq1 = 0.0;
                   if(verletList->getAdrRegionType()){
                     distsq1 = driftforce.sqr();
                   }
                   else{
                     distsq1 = driftforce[0]*driftforce[0];
                   }
                   if (distsq1 < min1sq) {
                        min1sq = distsq1;
                        mindriftforce = driftforce;
                   }
              }
              min1sq = sqrt(min1sq);   // distance to nearest adress particle or center

              if(verletList->getAdrRegionType()){
                mindriftforce = (1.0/min1sq)*mindriftforce; // normalized driftforce vector
                mindriftforce *= (0.5 * energydiff.find(&vp)->second); // get the energy differences which were calculated previously and put in drift force
                mindriftforce *= vp.lambdaDeriv();
                vp.force() += mindriftforce;
              }
              else{
                real mindriftforceX = (1.0/min1sq)*mindriftforce[0]; // normalized driftforce vector
                mindriftforceX *= (0.5 * energydiff.find(&vp)->second); // get the energy differences which were calculated previously and put in drift force
                mindriftforceX *= vp.lambdaDeriv();
                Real3D driftforceadd(mindriftforceX,0.0,0.0);
                vp.force() += driftforceadd;
              }
              vp.drift() += 0.5 * energydiff.find(&vp)->second;
          }

      }

      energydiff.clear();  // clear the energy difference map
    }

    // Energy calculation does currently only work if integrator.run( ) (also with 0) and decompose have been executed before. This is due to the initialization of the tuples.
    template < typename _PotentialAT, typename _PotentialCG >
    inline real
    VerletListHadressInteractionTemplate < _PotentialAT, _PotentialCG >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy of the Verlet list pairs");

      // REMOVE FOR IDEAL GAS
      real e = 0.0;
      for (PairList::Iterator it(verletList->getPairs());
           it.isValid(); ++it) {
          Particle &p1 = *it->first;
          Particle &p2 = *it->second;
          int type1 = p1.type();
          int type2 = p2.type();
          const PotentialCG &potential = getPotentialCG(type1, type2);
          e += potential._computeEnergy(p1, p2);
      }
      // REMOVE FOR IDEAL GAS

      for (PairList::Iterator it(verletList->getAdrPairs());
           it.isValid(); ++it) {
          Particle &p1 = *it->first;
          Particle &p2 = *it->second;
          real w1 = p1.lambda();
          real w2 = p2.lambda();
          real w12 = (w1 + w2)/2.0;

          // REMOVE FOR IDEAL GAS
          int type1 = p1.type();
          int type2 = p2.type();
          const PotentialCG &potentialCG = getPotentialCG(type1, type2);
          e += (1.0-w12)*potentialCG._computeEnergy(p1, p2);
          // REMOVE FOR IDEAL GAS

          FixedTupleListAdress::iterator it3;
          FixedTupleListAdress::iterator it4;
          it3 = fixedtupleList->find(&p1);
          it4 = fixedtupleList->find(&p2);

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

                      // AT energies
                      const PotentialAT &potentialAT = getPotentialAT(p3.type(), p4.type());
                      e += w12*potentialAT._computeEnergy(p3, p4);

                  }
              }
          }
      }

      real esum;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, e, esum, std::plus<real>());
      return esum;
    }

    template < typename _PotentialAT, typename _PotentialCG >
    inline real
    VerletListHadressInteractionTemplate < _PotentialAT, _PotentialCG >::
    computeEnergyDeriv() {
      std::cout << "Warning! At the moment computeEnergyDeriv() in VerletListHadressInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _PotentialAT, typename _PotentialCG >
    inline real
    VerletListHadressInteractionTemplate < _PotentialAT, _PotentialCG >::
    computeEnergyAA() {
      LOG4ESPP_INFO(theLogger, "compute total AA energy of the Verlet list pairs");

      real e = 0.0;
      for (PairList::Iterator it(verletList->getPairs());
           it.isValid(); ++it) {
          Particle &p1 = *it->first;
          Particle &p2 = *it->second;
          FixedTupleListAdress::iterator it3;
          FixedTupleListAdress::iterator it4;
          it3 = fixedtupleList->find(&p1);
          it4 = fixedtupleList->find(&p2);

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

                      // AT energies
                      const PotentialAT &potentialAT = getPotentialAT(p3.type(), p4.type());
                      e += potentialAT._computeEnergy(p3, p4);

                  }
              }
          }
      }

      for (PairList::Iterator it(verletList->getAdrPairs());
           it.isValid(); ++it) {
          Particle &p1 = *it->first;
          Particle &p2 = *it->second;
          FixedTupleListAdress::iterator it3;
          FixedTupleListAdress::iterator it4;
          it3 = fixedtupleList->find(&p1);
          it4 = fixedtupleList->find(&p2);

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

                      // AT energies
                      const PotentialAT &potentialAT = getPotentialAT(p3.type(), p4.type());
                      e += potentialAT._computeEnergy(p3, p4);

                  }
              }
          }
      }

      real esum;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, e, esum, std::plus<real>());
      return esum;
    }

    template < typename _PotentialAT, typename _PotentialCG > inline real
    VerletListHadressInteractionTemplate < _PotentialAT, _PotentialCG >::
    computeEnergyAA(int atomtype) {
      std::cout << "Warning! At the moment computeEnergyAA(int atomtype) in VerletListHadressInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _PotentialAT, typename _PotentialCG >
    inline real
    VerletListHadressInteractionTemplate < _PotentialAT, _PotentialCG >::
    computeEnergyCG() {
      LOG4ESPP_INFO(theLogger, "compute total CG energy of the Verlet list pairs");

      real e = 0.0;
      for (PairList::Iterator it(verletList->getPairs());
           it.isValid(); ++it) {
          Particle &p1 = *it->first;
          Particle &p2 = *it->second;
          int type1 = p1.type();
          int type2 = p2.type();
          const PotentialCG &potential = getPotentialCG(type1, type2);
          e += potential._computeEnergy(p1, p2);
      }

      for (PairList::Iterator it(verletList->getAdrPairs());
           it.isValid(); ++it) {
          Particle &p1 = *it->first;
          Particle &p2 = *it->second;
          int type1 = p1.type();
          int type2 = p2.type();
          const PotentialCG &potentialCG = getPotentialCG(type1, type2);
          e += potentialCG._computeEnergy(p1, p2);
      }

      real esum;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, e, esum, std::plus<real>());
      return esum;
    }

    template < typename _PotentialAT, typename _PotentialCG > inline real
    VerletListHadressInteractionTemplate < _PotentialAT, _PotentialCG >::
    computeEnergyCG(int atomtype) {
      std::cout << "Warning! At the moment computeEnergyCG(int atomtype) in VerletListHadressInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _PotentialAT, typename _PotentialCG > inline void
    VerletListHadressInteractionTemplate < _PotentialAT, _PotentialCG >::
    computeVirialX(std::vector<real> &p_xx_total, int bins) {
      LOG4ESPP_INFO(theLogger, "compute virial p_xx of the pressure tensor slabwise");

      std::set<Particle*> cgZone = verletList->getCGZone();
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
              for (std::vector<Particle*>::iterator it2 = atList.begin();
                                   it2 != atList.end(); ++it2) {
                  Particle *at = *it2;
                  cmp += at->mass() * at->position();
              }
              cmp /= vp.getMass();

              // update (overwrite) the position and velocity of the VP
              vp.position() = cmp;
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
              for (std::vector<Particle*>::iterator it2 = atList.begin();
                                   it2 != atList.end(); ++it2) {
                  Particle &at = **it2;
                  cmp += at.mass() * at.position();
              }
              cmp /= vp.getMass();

              // update (overwrite) the position and velocity of the VP
              vp.position() = cmp;
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
      }

    }



    template < typename _PotentialAT, typename _PotentialCG > inline real
    VerletListHadressInteractionTemplate < _PotentialAT, _PotentialCG >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute the virial for the Verlet List");

      real w = 0.0;

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
          //Real3D dist;
          //bc.getMinimumImageVectorBox(dist, p1.position(), p2.position());

          w += dist * force;
        }
      }

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
                if (w12 != 1.0) { // calculate VP force if both VP are outside AT region (CG-HY, HY-HY)
                    if (potentialCG._computeForce(forcevp, p1, p2)) {
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
                         //bc.getMinimumImageVectorBox(dist, p3.position(), p3.position());
                         w += dist * forceat;
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
         }
      }

      real wsum;
      wsum = 0.0;
      boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
      return wsum;
    }

    template < typename _PotentialAT, typename _PotentialCG > inline void
    VerletListHadressInteractionTemplate < _PotentialAT, _PotentialCG >::
    computeVirialTensor(Tensor& w) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the Verlet List");

      Tensor wlocal(0.0);
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
      w += wsum;
    }

    template < typename _PotentialAT, typename _PotentialCG > inline void
    VerletListHadressInteractionTemplate < _PotentialAT, _PotentialCG >::
    computeVirialTensor(Tensor& w, real z) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the Verlet List");

      std::cout << "Warning! At the moment IK computeVirialTensor in VerletListHAdress does'n work"<<std::endl;
    }


    template < typename _PotentialAT, typename _PotentialCG > inline void
    VerletListHadressInteractionTemplate < _PotentialAT, _PotentialCG >::
    computeVirialTensor(Tensor *w, int n) {
      std::cout << "Warning! At the moment IK computeVirialTensor in VerletListHAdress does'n work"<<std::endl;
    }

    template < typename _PotentialAT, typename _PotentialCG >
    inline real
    VerletListHadressInteractionTemplate< _PotentialAT, _PotentialCG >::
    getMaxCutoff() {
      real cutoff = 0.0;
      for (int i = 0; i < ntypes; i++) {
        for (int j = 0; j < ntypes; j++) {
          cutoff = std::max(cutoff, getPotentialCG(i, j).getCutoff());
        }
      }
      return cutoff;
    }
  }
}
#endif
