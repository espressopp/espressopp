/*
  Copyright (C) 2018
      Max Planck Institute for Polymer Research

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
#ifndef _INTERACTION_VERLETLISTHADRESSCGINTERACTIONTEMPLATE_HPP
#define _INTERACTION_VERLETLISTHADRESSCGINTERACTIONTEMPLATE_HPP

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
#include <boost/unordered_map.hpp>

namespace espressopp {
  namespace interaction {
    template < typename _Potential >
    class VerletListHadressCGInteractionTemplate: public Interaction {

      protected:
        typedef _Potential Potential;

      public:
        VerletListHadressCGInteractionTemplate
        (shared_ptr<VerletListAdress> _verletList, shared_ptr<FixedTupleListAdress> _fixedtupleList)
          : verletList(_verletList), fixedtupleList(_fixedtupleList) {

          potentialArray = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());

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
        setPotential(int type1, int type2, const Potential &potential) {
          // typeX+1 because i<ntypes
          ntypes = std::max(ntypes, std::max(type1+1, type2+1));

          potentialArray.at(type1, type2) = potential;
          if (type1 != type2) { // add potential in the other direction
            potentialArray.at(type2, type1) = potential;
          }
        }

        Potential &getPotential(int type1, int type2) {
          return potentialArray.at(type1, type2);
        }

        shared_ptr<Potential> getPotentialPtr(int type1, int type2) {
          return  make_shared<Potential>(potentialArray.at(type1, type2));
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
        virtual int bondType() {
          return NonbondedSlow;
        }

      protected:
        int ntypes;
        shared_ptr<VerletListAdress> verletList;
        shared_ptr<FixedTupleListAdress> fixedtupleList;
        esutil::Array2D<Potential, esutil::enlarge> potentialArray;

        // AdResS stuff
        real pidhy2; // pi / (dhy * 2)
        real dexdhy; // dex + dhy
        real dexdhy2; // dexdhy^2
        real dex;
        real dhy;
        real dex2; // dex^2
        boost::unordered_map<Particle*, real> energydiff;  // Energydifference V_AA - V_CG map for particles in hybrid region for drift term calculation in H-AdResS
        std::set<Particle*> adrZone;  // Virtual particles in AdResS zone (HY and AT region)
        std::set<Particle*> cgZone;

    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _Potential > inline void
    VerletListHadressCGInteractionTemplate < _Potential >::
    addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by the Verlet List");

      std::set<Particle*> cgZone = verletList->getCGZone();
      std::set<Particle*> adrZone = verletList->getAdrZone();

      for (std::set<Particle*>::iterator it=adrZone.begin();
           it != adrZone.end(); ++it) {
        Particle &p = **it;
        // intitialize energy diff AA-CG
        if (p.lambda()<0.9999999 && p.lambda()>0.0000001) {
          energydiff[&p]=0.0;
        }
        // energydiff[&p]=0.0;
      }


      // Pairs not inside the AdResS Zone (CG region)
      for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it) {

        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.type();
        int type2 = p2.type();

        const Potential &potential = getPotential(type1, type2);
        Real3D force(0.0, 0.0, 0.0);
        if(potential._computeForce(force, p1, p2)) {
          p1.force() += force;
          p2.force() -= force;
        }
      }

      // Compute forces (AT and VP) of Pairs inside AdResS zone
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
        const Potential &potential = getPotential(type1, type2);
        Real3D forcevp(0.0, 0.0, 0.0);
        if (w12 != 1.0) { // calculate VP force if both VP are outside AT region (CG-HY, HY-HY)
          if (potential._computeForce(forcevp, p1, p2)) {
            forcevp *= (1.0 - w12);
            p1.force() += forcevp;
            p2.force() -= forcevp;
          }

          // H-AdResS - Drift Term part 1
          // Compute CG energies of particles in the hybrid and store and add up in map energydiff
          if (w12 != 0.0) {   //at least one particle in hybrid region => need to do the energy calculation
            real energyvp = potential._computeEnergy(p1, p2);
            if (w1 != 0.0) {   // if particle one is in hybrid region
              energydiff[&p1] += energyvp;   // add CG energy for virtual particle 1
            }
            if (w2 != 0.0) {   // if particle two is in hybrid region
              energydiff[&p2] += energyvp;   // add CG energy for virtual particle 2
            }
          }
        }
      }

      // H-AdResS - Drift Term part 3
      // Iterate over all particles in the hybrid region and calculate drift force
      for (std::set<Particle*>::iterator it=adrZone.begin();
           it != adrZone.end(); ++it) {   // Iterate over all particles
        Particle &vp = **it;
        real w = vp.lambda();

        if(w<0.9999999 && w>0.0000001) {  //   only chose those in the hybrid region
          // calculate distance to nearest adress particle or center
          std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
          Real3D pa = **it2; // position of adress particle
          Real3D mindriftforce(0.0, 0.0, 0.0);
          verletList->getSystem()->bc->getMinimumImageVector(mindriftforce, vp.position(), pa);
          real min1sq = 0.0;
          if(verletList->getAdrRegionType()) {
            min1sq = mindriftforce.sqr();
          }
          else {
            min1sq = mindriftforce[0]*mindriftforce[0];
          }
          ++it2;
          for (; it2 != verletList->getAdrPositions().end(); ++it2) {
            pa = **it2;
            Real3D driftforce(0.0, 0.0, 0.0);
            verletList->getSystem()->bc->getMinimumImageVector(driftforce, vp.position(), pa);
            real distsq1 = 0.0;
            if(verletList->getAdrRegionType()) {
              distsq1 = driftforce.sqr();
            }
            else {
              distsq1 = driftforce[0]*driftforce[0];
            }
            if (distsq1 < min1sq) {
              min1sq = distsq1;
              mindriftforce = driftforce;
            }
          }
          min1sq = sqrt(min1sq);   // distance to nearest adress particle or center

          if(verletList->getAdrRegionType()) {
            mindriftforce = (1.0/min1sq)*mindriftforce; // normalized driftforce vector
            // mindriftforce *= (0.5 * energydiff.find(&vp)->second); // get the energy differences which were calculated previously and put in drift force
            mindriftforce *= (0.5 * energydiff[&vp]);
            mindriftforce *= vp.lambdaDeriv();
            vp.force() += mindriftforce;
          }
          else {
            real mindriftforceX = (1.0/min1sq)*mindriftforce[0]; // normalized driftforce vector
            // mindriftforceX *= (0.5 * energydiff.find(&vp)->second); // get the energy differences which were calculated previously and put in drift force
            mindriftforceX *= (0.5 * energydiff[&vp]);
            mindriftforceX *= vp.lambdaDeriv();
            Real3D driftforceadd(mindriftforceX,0.0,0.0);
            vp.force() += driftforceadd;
          }
          vp.drift() += 0.5 * energydiff.find(&vp)->second;
        }

      }

      // energydiff.clear();  // clear the energy difference map
    }

    // Energy calculation does currently only work if integrator.run( ) (also with 0) and decompose have been executed before. This is due to the initialization of the tuples.
    template < typename _Potential >
    inline real
    VerletListHadressCGInteractionTemplate < _Potential >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy of the Verlet list pairs");

      real e = 0.0;
      for (PairList::Iterator it(verletList->getPairs());
           it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.type();
        int type2 = p2.type();
        const Potential &potential = getPotential(type1, type2);
        e += potential._computeEnergy(p1, p2);
      }

      for (PairList::Iterator it(verletList->getAdrPairs());
           it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        real w1 = p1.lambda();
        real w2 = p2.lambda();
        real w12 = (w1 + w2)/2.0;
        int type1 = p1.type();
        int type2 = p2.type();
        const Potential &potential = getPotential(type1, type2);
        e += (1.0-w12)*potential._computeEnergy(p1, p2);
      }

      real esum;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, e, esum, std::plus<real>());
      return esum;
    }

    template < typename _Potential >
    inline real
    VerletListHadressCGInteractionTemplate < _Potential >::
    computeEnergyDeriv() {
      std::cout << "Warning! At the moment computeEnergyDeriv in VerletListHadressCGInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline real
    VerletListHadressCGInteractionTemplate < _Potential >::
    computeEnergyAA() {
      std::cout << "Warning! At the moment computeEnergyAA in VerletListHadressCGInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline real
    VerletListHadressCGInteractionTemplate < _Potential >::
    computeEnergyAA(int atomtype) {
      std::cout << "Warning! At the moment computeEnergyAA(int atomtype) in VerletListHadressCGInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline real
    VerletListHadressCGInteractionTemplate < _Potential >::
    computeEnergyCG() {
      std::cout << "Warning! At the moment computeEnergyCG in VerletListHadressCGInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline real
    VerletListHadressCGInteractionTemplate < _Potential >::
    computeEnergyCG(int atomtype) {
      std::cout << "Warning! At the moment computeEnergyCG(int atomtype) in VerletListHadressCGInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline void
    VerletListHadressCGInteractionTemplate < _Potential >::
    computeVirialX(std::vector<real> &p_xx_total, int bins) {
      std::cout << "Warning! At the moment computeVirialX in VerletListHadressCGInteractionTemplate does'n work"<<std::endl;
    }

    template < typename _Potential > inline real
    VerletListHadressCGInteractionTemplate < _Potential >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute the virial for the Verlet List");
      std::cout << "Warning! At the moment computeVirial in VerletListHadressCGInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline void
    VerletListHadressCGInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor& w) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the Verlet List");
      std::cout << "Warning! At the moment computeVirialTensor in VerletListHadressCGInteractionTemplate does'n work"<<std::endl;
    }

    template < typename _Potential > inline void
    VerletListHadressCGInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor& w, real z) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the Verlet List");
      std::cout << "Warning! At the moment computeVirialTensor in VerletListHadressCGInteractionTemplate does'n work"<<std::endl;
    }

    template < typename _Potential > inline void
    VerletListHadressCGInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor *w, int n) {
      std::cout << "Warning! At the moment computeVirialTensor in VerletListHadressCGInteractionTemplate does'n work"<<std::endl;
    }

    template < typename _Potential >
    inline real
    VerletListHadressCGInteractionTemplate< _Potential >::
    getMaxCutoff() {
      real cutoff = 0.0;
      for (int i = 0; i < ntypes; i++) {
        for (int j = 0; j < ntypes; j++) {
          cutoff = std::max(cutoff, getPotential(i, j).getCutoff());
        }
      }
      return cutoff;
    }
  }
}
#endif
