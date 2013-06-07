// ESPP_CLASS 
#ifndef _INTERACTION_VERLETLISTHADRESSINTERACTIONTEMPLATE_HPP
#define _INTERACTION_VERLETLISTHADRESSINTERACTIONTEMPLATE_HPP

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
    class VerletListHadressInteractionTemplate: public Interaction {
    
    protected:
      typedef _PotentialAT PotentialAT;
      typedef _PotentialCG PotentialCG;
    
    public:
      VerletListHadressInteractionTemplate
      (shared_ptr<VerletListAdress> _verletList, shared_ptr<FixedTupleList> _fixedtupleList)
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
      virtual void computeVirialTensor(Tensor& w, real z);
      virtual void computeVirialTensor(Tensor *w, int n);
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
      real dhy;
      real dex2; // dex^2
      std::map<Particle*, real> weights;
      std::map<Particle*, real> energydiff;  // Energydifference V_AA - V_CG map for particles in hybrid region for drift term calculation in H-AdResS
      std::set<Particle*> adrZone;  // Virtual particles in AdResS zone (HY and AT region)

      // AdResS Weighting function
      real weight(real distanceSqr){
          if (dex2 > distanceSqr) return 1.0;
          else if (dexdhy2 < distanceSqr) return 0.0;
          else {
              real argument = sqrt(distanceSqr) - dex;
              return 1.0-(30.0/(pow(dhy, 5.0)))*(1.0/5.0*pow(argument, 5.0)-dhy/2.0*pow(argument, 4.0)+1.0/3.0*pow(argument, 3.0)*dhy*dhy);
              //return pow(cos(pidhy2 * argument),2.0); // for cosine squared weighting function
          }
      }
      real weightderivative(real distance){
          real argument = distance - dex;
          return -(30.0/(pow(dhy, 5.0)))*(pow(argument, 4.0)-2.0*dhy*pow(argument, 3.0)+argument*argument*dhy*dhy);
          //return -pidhy2 * 2.0 * cos(pidhy2*argument) * sin(pidhy2*argument); // for cosine squared weighting function
      }

      // Compute center of mass and set the weights for virtual particles in AdResS zone (HY and AT region).
      void makeWeights(){
          adrZone = verletList->getAdrZone();
          for (std::set<Particle*>::iterator it=adrZone.begin();
                  it != adrZone.end(); ++it) {

              Particle &vp = **it;

              FixedTupleList::iterator it3;
              it3 = fixedtupleList->find(&vp);

              if (it3 != fixedtupleList->end()) {

                  std::vector<Particle*> atList;
                  atList = it3->second;

                  // Compute center of mass
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
                  //Real3D d1 = vp.position() - pa;                                               // X SPLIT VS SPHERE CHANGE
                  real d1 = vp.position()[0] - pa[0];                                             // X SPLIT VS SPHERE CHANGE
                  real min1sq = d1*d1; // d1.sqr(); // set min1sq before loop                     // X SPLIT VS SPHERE CHANGE
                  ++it2;
                  for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                       pa = **it2;
                       //d1 = vp.position() - pa;                                                 // X SPLIT VS SPHERE CHANGE
                       d1 = vp.position()[0] - pa[0];                                             // X SPLIT VS SPHERE CHANGE
                       real distsq1 = d1*d1; //d1.sqr();
                       //std::cout << pa << " " << sqrt(distsq1) << "\n";
                       if (distsq1 < min1sq) min1sq = distsq1;
                  }
                  
                  real w = weight(min1sq);
                  
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
        }
      
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _PotentialAT, typename _PotentialCG > inline void
    VerletListHadressInteractionTemplate < _PotentialAT, _PotentialCG >::
    addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by the Verlet List");

                  
      // Update CG particle position according to center of masses of AT particles (CG region).
      
      // Note (Karsten): This is a different approach compared to Force-AdResS. In
      // Force-AdResS, we overwrite intra-molecular rotations and vibrations in the CG zone. This leads to failures in the kinetic energy. However, in Force-AdResS there is no energy
      // conservation anyway. In Force-AdResS we calculate CG forces/velocities and distribute them to AT particles. 
      // In contrast, in H-AdResS, we calculate AT forces from intra-molecular interactions and inter-molecular center-of-mass interactions and just update the positions
      // of the center-of-mass CG particles.
      std::set<Particle*> cgZone = verletList->getCGZone();
      for (std::set<Particle*>::iterator it=cgZone.begin();
              it != cgZone.end(); ++it) {

          Particle &vp = **it;

          FixedTupleList::iterator it3;
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
      // std::set<Particle*> cgZone = verletList->getCGZone();
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
                    // at.velocity() = vp.velocity(); // overwrite velocity
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
      //std::set<Particle*> adrZone = verletList->getAdrZone();
      makeWeights();
      
      // Compute forces (AT and VP) of Pairs inside AdResS zone
      int count = 0;
      for (PairList::Iterator it(verletList->getAdrPairs()); it.isValid(); ++it) {

         // these are the two VP interacting
         Particle &p1 = *it->first;
         Particle &p2 = *it->second;

         // read weights
         real w1 = weights.find(&p1)->second;
         real w2 = weights.find(&p2)->second;
         real w12 = (w1 + w2)/2.0;  // H-AdResS

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
         /*
         else {
             std::cout << "skipping VP forces...\n";
         }*/

         // force between AT particles
         if (w12 != 0.0) { // calculate AT force if both VP are outside CG region (HY-HY, HY-AT, AT-AT)
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
                            
                             if(force[0]!=0.0 || force[1]!=0.0 || force[2]!=0.0){count +=1;} 
                             
                             force *= w12;
                             p3.force() += force;
                             p4.force() -= force;
                         }
                         
                         // H-AdResS - Drift Term part 2
                         // Compute AT energies of particles in the hybrid and store and subtract in map energydiff
                         if(w12!=1.0){   //at least one particle in hybrid region => need to do the energy calculation
                             real energyat = potentialAT._computeEnergy(p3, p4);   
                             if(w1!=1.0){   // if particle one is in hybrid region
                                    //energydiff.find(&p1)->second -= energyat;
                                    energydiff[&p1] -= energyat;   // subtract AT energy for virtual particle 1
                             }
                             if(w2!=1.0){   // if particle two is in hybrid region
                                    //energydiff.find(&p2)->second -= energyat;
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
            std::cout << "Count Adr " << count << std::endl;
      
      // H-AdResS - Drift Term part 3
      // Iterate over all particles in the hybrid region and calculate drift force
      for (std::set<Particle*>::iterator it=adrZone.begin();
        it != adrZone.end(); ++it) {   // Iterate over all particles
          Particle &vp = **it;
          real w = weights.find(&vp)->second;
                  
          if(w!=1.0 && w!=0.0){   //   only chose those in the hybrid region
              
              // calculate distance to nearest adress particle or center
              std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
              Real3D pa = **it2; // position of adress particle
              //Real3D mindriftforce = vp.position() - pa;                                                 // X SPLIT VS SPHERE CHANGE
              real mindriftforce = vp.position()[0] - pa[0];                                               // X SPLIT VS SPHERE CHANGE
              real min1sq = mindriftforce*mindriftforce; // mindriftforce.sqr(); // set min1sq before loop
              ++it2;
              for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                   pa = **it2;
                   //Real3D driftforce = vp.position() - pa;                                                 // X SPLIT VS SPHERE CHANGE
                   real driftforce = vp.position()[0] - pa[0];                                               // X SPLIT VS SPHERE CHANGE
                   real distsq1 = driftforce*driftforce; //driftforce.sqr();
                   //std::cout << pa << " " << sqrt(distsq1) << "\n";
                   if (distsq1 < min1sq) {
                        min1sq = distsq1;
                        mindriftforce = driftforce;
                   }
              }
              min1sq = sqrt(min1sq);   // distance to nearest adress particle or center
              mindriftforce = (1.0/min1sq)*mindriftforce;  // normalized driftforce vector
              mindriftforce *= weightderivative(min1sq);  // multiplication with derivative of the weighting function
              mindriftforce *= energydiff.find(&vp)->second;   // get the energy differences which were calculated previously and put in drift force
              //vp.force() += mindriftforce;   // add drift force to virtual particles                                                           // X SPLIT VS SPHERE CHANGE
              Real3D driftforceadd(mindriftforce,0.0,0.0);                                                                                       // X SPLIT VS SPHERE CHANGE
              vp.force() += driftforceadd;                                                                                                       // X SPLIT VS SPHERE CHANGE
          }
          
      }
      
      energydiff.clear();  // clear the energy difference map
      
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
                //std::cout << "Force of atomistic particle (AdResS sim.) with id " << at.id() << " is: " << at.force() << "\n";
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
    
    // Energy calculation does currently only work if integrator.run( ) (also with 0) and decompose have been executed before. This is due to the initialization of the tuples.
    template < typename _PotentialAT, typename _PotentialCG >
    inline real
    VerletListHadressInteractionTemplate < _PotentialAT, _PotentialCG >::
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
          //std::cout << "Energy calculation CG region done" << "\n";
      }
      //std::cout << "Energy CG region:" << e << "\n";
      makeWeights();
      int counter = 0;
      for (PairList::Iterator it(verletList->getAdrPairs()); 
           it.isValid(); ++it) {
          counter += 1;
          Particle &p1 = *it->first;
          Particle &p2 = *it->second;                           
          real w1 = weights.find(&p1)->second;
          real w2 = weights.find(&p2)->second;
          real w12 = (w1 + w2)/2.0;
          int type1 = p1.type();
          int type2 = p2.type();
          const PotentialCG &potentialCG = getPotentialCG(type1, type2);
          e += (1.0-w12)*potentialCG._computeEnergy(p1, p2);
          //std::cout << "CG Energy calculation AT/HY region done:" << e << "\n";
          
          FixedTupleList::iterator it3;
          FixedTupleList::iterator it4;
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
                      //std::cout << "Nonbonded AT energy: " << e << "\n";
                      //counter += 1;
                      //std::cout << "AT Energy calculation AT/HY region done, w12 = " << w12 << " and w1 = " << w1 << " and w2 = " << w2 << " okay...\n";
                  }                  
              }              
          }          
      }
      //std::cout << "Energy AT + CG region:" << e << "\n";
      std::cout << "Total number of pairs in calculation (AdResS):" << counter << "\n";        
      real esum;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, e, esum, std::plus<real>());
      return esum;      
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
      boost::mpi::all_reduce(*mpiWorld, wlocal, wsum, std::plus<Tensor>());
      w += wsum;
    }
 
    template < typename _PotentialAT, typename _PotentialCG > inline void
    VerletListHadressInteractionTemplate < _PotentialAT, _PotentialCG >::
    computeVirialTensor(Tensor& w, real z) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the Verlet List");

      std::cout << "Warning! At the moment IK computeVirialTensor in VerletListHAdress does'n work"<<std::endl;
      
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
      boost::mpi::all_reduce(*mpiWorld, wlocal, wsum, std::plus<Tensor>());
      w += wsum;
       */
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
