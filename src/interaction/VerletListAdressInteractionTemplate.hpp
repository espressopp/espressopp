// ESPP_CLASS 
#ifndef _INTERACTION_VERLETLISTADRESSINTERACTIONTEMPLATE_HPP
#define _INTERACTION_VERLETLISTADRESSINTERACTIONTEMPLATE_HPP

//#include <typeinfo>

#include "System.hpp"
#include "bc/BC.hpp"

#include "types.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "VerletListAdress.hpp"
#include "FixedTupleList.hpp"
#include "esutil/Array2D.hpp"

namespace espresso {
  namespace interaction {
    template < typename _PotentialAT, typename _PotentialCG >
    class VerletListAdressInteractionTemplate: public Interaction {
    
    protected:
      typedef _PotentialAT PotentialAT;
      typedef _PotentialCG PotentialCG;
    
    public:
      VerletListAdressInteractionTemplate
      (shared_ptr<VerletListAdress> _verletList, shared_ptr<FixedTupleList> _fixedtupleList)
                : verletList(_verletList), fixedtupleList(_fixedtupleList) {

          potentialArrayAT = esutil::Array2D<PotentialAT, esutil::enlarge>(0, 0, PotentialAT());
          potentialArrayCG = esutil::Array2D<PotentialCG, esutil::enlarge>(0, 0, PotentialCG());

          // AdResS stuff
          pidhy2 = M_PI/(verletList->getHy() * 2);
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
      setFixedTupleList(shared_ptr<FixedTupleList> _fixedtupleList) {
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
      virtual real computeVirial();
      virtual void computeVirialTensor(Tensor& w);
      virtual real getMaxCutoff();
      virtual int bondType() { return Nonbonded; }

    protected:
      int ntypes;
      shared_ptr<VerletListAdress> verletList;
      shared_ptr<FixedTupleList> fixedtupleList;
      esutil::Array2D<PotentialAT, esutil::enlarge> potentialArrayAT;
      esutil::Array2D<PotentialCG, esutil::enlarge> potentialArrayCG;

      // AdResS stuff
      real pidhy2; // pi / (dhy * 2)
      real dexdhy; // dex + dhy
      real dexdhy2; // dexdhy^2
      real dex;
      real dex2; // dex^2
      std::map<Particle*, real> weights;

    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _PotentialAT, typename _PotentialCG > inline void
    VerletListAdressInteractionTemplate < _PotentialAT, _PotentialCG >::
    addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by the Verlet List");

      // Pairs not inside the AdResS Zone (CG region)
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

      // Loop over CG particles and overwrite AT forces and velocity.
      // This makes the AT particles move along with CG particles.
      std::set<Particle*> cgZone = verletList->getCGZone();
      for (std::set<Particle*>::iterator it=cgZone.begin();
                    it != cgZone.end(); ++it) {

            Particle &vp = **it;

            FixedTupleList::iterator it3;
            it3 = fixedtupleList->find(&vp);

            if (it3 != fixedtupleList->end()) {

                std::vector<Particle*> atList1;
                atList1 = it3->second;

                Real3D vpfm = vp.force() / vp.getMass();
                for (std::vector<Particle*>::iterator itv = atList1.begin();
                        itv != atList1.end(); ++itv) {
                    Particle &at = **itv;
                    at.velocity() = vp.velocity(); // overwrite velocity
                    at.force() += at.mass() * vpfm;
                }

            }
            else { // this should not happen
                std::cout << " VP particle not found in tuples: " << vp.id() << "-" << vp.ghost();
                exit(1);
                return;
            }
      }

      // Compute center of mass and weights for virtual particles in Adress zone (HY and AT region).
      std::set<Particle*> adrZone = verletList->getAdrZone();
      for (std::set<Particle*>::iterator it=adrZone.begin();
              it != adrZone.end(); ++it) {

          Particle &vp = **it;

          FixedTupleList::iterator it3;
          it3 = fixedtupleList->find(&vp);

          if (it3 != fixedtupleList->end()) {

              std::vector<Particle*> atList;
              atList = it3->second;

              // compute center of mass
              Real3D cmp(0.0, 0.0, 0.0); // center of mass position
              Real3D cmv(0.0, 0.0, 0.0); // center of mass velocity
              //real M = vp.getMass(); // sum of mass of AT particles
              for (std::vector<Particle*>::iterator it2 = atList.begin();
                                   it2 != atList.end(); ++it2) {
                  Particle &at = **it2;
                  //Real3D d1 = at.position() - vp.position();
                  //Real3D d1;
                  //verletList->getSystem()->bc->getMinimumImageVectorBox(d1, at.position(), vp.position());
                  //cmp += at.mass() * d1;

                  cmp += at.mass() * at.position();
                  cmv += at.mass() * at.velocity();
              }
              cmp /= vp.getMass();
              cmv /= vp.getMass();
              //cmp += vp.position(); // cmp is a relative position
              //std::cout << " cmp M: "  << M << "\n\n";
              //std::cout << "  moving VP to " << cmp << ", velocitiy is " << cmv << "\n";

              // update (overwrite) the position and velocity of the VP
              vp.position() = cmp;
              vp.velocity() = cmv;

              // calculate distance to nearest adress particle or center
              std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
              Real3D pa = **it2; // position of adress particle
              Real3D d1 = vp.position() - pa;
              real min1sq = d1.sqr(); // set min1sq before loop
              ++it2;
              for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                   pa = **it2;
                   d1 = vp.position() - pa;
                   real distsq1 = d1.sqr();
                   //std::cout << pa << " " << sqrt(distsq1) << "\n";
                   if (distsq1 < min1sq) min1sq = distsq1;
              }

              //real min1 = sqrt(min1sq);
              //std::cout << vp.id() << " min: " << min1 << "\n";
              //std::cout << vp.id() << " dex: " << dex << "\n";
              //std::cout << vp.id() << " dex+dhy: " << dexdhy << "\n";

              // calculate weight and write it in the map
              real w;
              if (dex2 > min1sq) w = 1;
              else if (dexdhy2 < min1sq) w = 0;
              else {
                   w = cos(pidhy2 * (sqrt(min1sq) - dex));
                   w *= w;
              }

              weights.insert(std::make_pair(&vp, w));

              //if (w1 == 1 || w2 == 1) std::cout << p1.id() << " ";
              //std::cout << vp.id() << " weight: " << w << "\n";
          }
          else { // this should not happen
              std::cout << " VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
              std::cout << " (" << vp.position() << ")\n";
              exit(1);
              return;
          }
      }


      // Compute forces (AT and VP) of Pairs inside AdResS zone
      for (PairList::Iterator it(verletList->getAdrPairs()); it.isValid(); ++it) {

         // these are the two VP interacting
         Particle &p1 = *it->first;
         Particle &p2 = *it->second;

         // read weights
         real w1 = weights.find(&p1)->second;
         real w2 = weights.find(&p2)->second;
         real w12 = w1 * w2;

         // force between VP particles
         int type1 = p1.type();
         int type2 = p2.type();
         const PotentialCG &potentialCG = getPotentialCG(type1, type2);
         Real3D forcevp(0.0, 0.0, 0.0);
         if (w12 != 1) { // calculate VP force if both VP are outside AT region (CG-HY, HY-HY)
             if(potentialCG._computeForce(forcevp, p1, p2)) {
                 forcevp *= (1 - w12);
                 p1.force() += forcevp;
                 p2.force() -= forcevp;
             }
         }
         /*
         else {
             std::cout << "skipping VP forces...\n";
         }*/

         // force between AT particles
         if (w12 != 0) { // calculate AT force if both VP are outside CG region (HY-HY, HY-AT, AT-AT)
             FixedTupleList::iterator it3;
             FixedTupleList::iterator it4;
             it3 = fixedtupleList->find(&p1);
             it4 = fixedtupleList->find(&p2);

             //std::cout << "Interaction " << p1.id() << " - " << p2.id() << "\n";
             if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

                 std::vector<Particle*> atList1;
                 std::vector<Particle*> atList2;
                 atList1 = it3->second;
                 atList2 = it4->second;

                 //std::cout << "AT forces ...\n";
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

      weights.clear();

      // distribute forces from VP to AT (HY and AT region)
      for (std::set<Particle*>::iterator it=adrZone.begin();
                it != adrZone.end(); ++it) {

        Particle &vp = **it;

        FixedTupleList::iterator it3;
        it3 = fixedtupleList->find(&vp);

        if (it3 != fixedtupleList->end()) {

            std::vector<Particle*> atList;
            atList = it3->second;

            // update force of AT particles belonging to a VP
            Real3D vpfm = vp.force() / vp.getMass();
            for (std::vector<Particle*>::iterator it2 = atList.begin();
                                 it2 != atList.end(); ++it2) {
                Particle &at = **it2;
                at.force() += at.mass() * vpfm;
            }
        }
        else { // this should not happen
            std::cout << " particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
            std::cout << " (" << vp.position() << ")\n";
            exit(1);
            return;
        }
      }
    }
    
    template < typename _PotentialAT, typename _PotentialCG >
    inline real
    VerletListAdressInteractionTemplate < _PotentialAT, _PotentialCG >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy of the Verlet list pairs");

      //std::cout << "compute energy of the Verlet list pairs" << "\n";
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

      //std::cout << "compute energy of the AdReSs pairs" << "\n";
      for (PairList::Iterator it(verletList->getAdrPairs());
             it.isValid(); ++it) {
              Particle &p1 = *it->first;
              Particle &p2 = *it->second;
              int type1 = p1.type();
              int type2 = p2.type();
              const PotentialCG &potential = getPotentialCG(type1, type2);
              e += potential._computeEnergy(p1, p2);
      }

      real esum;
      boost::mpi::reduce(*getVerletList()->getSystem()->comm, e, esum, std::plus<real>(), 0);
      return esum;
    }





    template < typename _PotentialAT, typename _PotentialCG > inline real
    VerletListAdressInteractionTemplate < _PotentialAT, _PotentialCG >::
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
      boost::mpi::reduce(*mpiWorld, w, wsum, std::plus<real>(), 0);
      return wsum; 
    }

    template < typename _PotentialAT, typename _PotentialCG > inline void
    VerletListAdressInteractionTemplate < _PotentialAT, _PotentialCG >::
    computeVirialTensor(Tensor& w) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the Verlet List");

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
          w += Tensor(dist, force);
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
                w += Tensor(dist, force);
              }
      }

    }
 
    template < typename _PotentialAT, typename _PotentialCG >
    inline real
    VerletListAdressInteractionTemplate< _PotentialAT, _PotentialCG >::
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
