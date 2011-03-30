// ESPP_CLASS 
#ifndef _INTERACTION_VERLETLISTADRESSINTERACTIONTEMPLATE_HPP
#define _INTERACTION_VERLETLISTADRESSINTERACTIONTEMPLATE_HPP

//#include <typeinfo>


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
    template < typename _Potential >
    class VerletListAdressInteractionTemplate: public Interaction {
    
    protected:
      typedef _Potential Potential;
    
    public:
      VerletListAdressInteractionTemplate
          (shared_ptr<VerletListAdress> _verletList)
          : verletList(_verletList) {
          potentialArray = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());

          // AdResS stuff
          pidhy2 = M_PI/(verletList->getHy() * 2);
          dex = verletList->getEx();
          dexdhy = dex + verletList->getHy();
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
      setPotential(int type1, int type2, const Potential &potential) {
          potentialArray.at(type1, type2) = potential;
          if (type1 != type2) { // add potential in the other direction
             potentialArray.at(type2, type1) = potential;
          }
      }

      Potential &getPotential(int type1, int type2) {
        return potentialArray.at(type1, type2);
      }

      virtual void addForces();
      virtual real computeEnergy();
      virtual real computeVirial();
      virtual void computeVirialTensor(Tensor& w);
      virtual real getMaxCutoff();

    protected:
      int ntypes;
      shared_ptr<VerletListAdress> verletList;
      shared_ptr<FixedTupleList> fixedtupleList;
      esutil::Array2D<Potential, esutil::enlarge> potentialArray;

      // AdResS stuff
      real pidhy2; // pi / (dhy * 2)
      real dexdhy; // dex + dhy
      real dex;

    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _Potential > inline void
    VerletListAdressInteractionTemplate < _Potential >::
    addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by the Verlet List");

      //std::cout << "add forces computed by the Verlet List" << "\n";
      for (PairList::Iterator it(verletList->getPairs()); 
	   it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.type();
        int type2 = p2.type();
        const Potential &potential = getPotential(type1, type2);

        Real3D force(0.0, 0.0, 0.0);

        if(potential._computeForce(force, p1, p2)) {
          p1.force() += force;
          p2.force() -= force;

          Real3D force1(0.0, 0.0, 0.0);
          Real3D force2(0.0, 0.0, 0.0);
          force1 = force * p1.getMass();
          force2 = force * p2.getMass();

          //std::cout << p1.id() << "-" << p2.id() << " force: (" << force << ")\n";

          // iterate through atomistic particles in fixedtuplelist
          // and add them the proportional amount of force (depending
          // on their mass) to make them move along with their CG particles
          FixedTupleList::iterator it3;
          FixedTupleList::iterator it4;
          it3 = fixedtupleList->find(p1.id());
          it4 = fixedtupleList->find(p2.id());

          if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

              std::vector<Particle*> atList1;
              std::vector<Particle*> atList2;
              atList1 = it3->second;
              atList2 = it4->second;

              for (std::vector<Particle*>::iterator itv = atList1.begin();
                      itv != atList1.end(); ++itv) {
                  Particle &p3 = **itv;
                  p3.force() += force1 * (1/p3.getMass());
              }

              for (std::vector<Particle*>::iterator itv2 = atList2.begin();
                      itv2 != atList2.end(); ++itv2) {
                  Particle &p4 = **itv2;
                  p4.force() -= force2 * (1/p4.getMass());
              }
              //std::cout << "\n";
          }

        }
      }

      // adress TODO
      //std::cout << "add forces computed by the AdResS List" << "\n";
      for (PairList::Iterator it(verletList->getAdrPairs()); it.isValid(); ++it) {
         Particle &p1 = *it->first;
         Particle &p2 = *it->second;
         int type1 = p1.type();
         int type2 = p2.type();
         const Potential &potential = getPotential(type1, type2);


         // calculate distance to nearest atomistic particle
         std::vector<Real3D*>::iterator it2 = verletList->getAtmPositions().begin();
         Real3D pa = **it2; // position of atomistic particle
         Real3D d1 = p1.position() - pa;
         Real3D d2 = p2.position() - pa;
         real distsq1 = d1.sqr();
         real distsq2 = d2.sqr();
         real min1 = distsq1;
         real min2 = distsq2;
         for (; it2 != verletList->getAtmPositions().end(); ++it2) {
             pa = **it2;
             d1 = p1.position() - pa;
             d2 = p2.position() - pa;
             distsq1 = d1.sqr();
             distsq2 = d2.sqr();
             if (distsq1 < min1) min1 = distsq1;
             if (distsq2 < min2) min2 = distsq2;
         }
         //std::cout << "("<< p1.id() << ", " << p2.id() << " min: " << sqrt(min1) << ", " << sqrt(min2) << ") ";

         // calculate weight
         real w1, w2;
         if (dex > min1) w1 = 1;
         else if (dexdhy < min1) w1 = 0;
         else {
             w1 = cos(pidhy2 * (min1 - dex));
             w1 *= w1;
         }

         if (dex > min2) w2 = 1;
         else if (dexdhy < min2) w2 = 0;
         else {
             w2 = cos(pidhy2 * (min2 - dex));
             w2 *= w2;
         }

         //if (w1 == 1 || w2 == 1) std::cout << p1.id() << " ";
         //std::cout << p1.id() << ": " << w1 << "\n";
         //std::cout << p2.id() << ": " << w2 << "\n";

         // iterate through atomistic particles in fixedtuplelist
         FixedTupleList::iterator it3;
         FixedTupleList::iterator it4;
         it3 = fixedtupleList->find(p1.id());
         it4 = fixedtupleList->find(p2.id());

         if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

             std::vector<Particle*> atList1;
             std::vector<Particle*> atList2;
             atList1 = it3->second;
             atList2 = it4->second;

             //std::cout << p1.id() << "-" << p2.id() << ":\n";
             for (std::vector<Particle*>::iterator itv = atList1.begin();
                     itv != atList1.end(); ++itv) {
                 for (std::vector<Particle*>::iterator itv2 = atList2.begin();
                                      itv2 != atList2.end(); ++itv2) {

                     Particle &p3 = **itv;
                     Particle &p4 = **itv2;

                     //std::cout << p3.id()  << " " << p4.id() << "\n";

                     const Potential &potential2 =
                             getPotential(p3.type(),p4.type());

                     // AT forces
                     Real3D force(0.0, 0.0, 0.0);
                     if(potential._computeForce(force, p3, p4)) {
                        p3.force() += force;
                        p4.force() -= force;
                     }
                 }
             }
             //std::cout << "\n";
         }



         // CG forces
         Real3D force(0.0, 0.0, 0.0);
         if(potential._computeForce(force, p1, p2)) {
            p1.force() += force;
            p2.force() -= force;
         }
      }
      //std::cout << "\n\n";
    }
    
    template < typename _Potential >
    inline real
    VerletListAdressInteractionTemplate < _Potential >::
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
        const Potential &potential = getPotential(type1, type2);
        e += potential._computeEnergy(p1, p2);
      }

      // adress TODO
      //std::cout << "compute energy of the AdReSs list pairs" << "\n";
      for (PairList::Iterator it(verletList->getAdrPairs());
             it.isValid(); ++it) {
              Particle &p1 = *it->first;
              Particle &p2 = *it->second;
              int type1 = p1.type();
              int type2 = p2.type();
              const Potential &potential = getPotential(type1, type2);
              e += potential._computeEnergy(p1, p2);
      }

      real esum;
      boost::mpi::reduce(*getVerletList()->getSystem()->comm, e, esum, std::plus<real>(), 0);
      return esum;
    }





    template < typename _Potential > inline real
    VerletListAdressInteractionTemplate < _Potential >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute the virial for the Verlet List");
      

      real w = 0.0;
      for (PairList::Iterator it(verletList->getPairs());                
           it.isValid(); ++it) {                                         
        Particle &p1 = *it->first;                                       
        Particle &p2 = *it->second;                                      
        int type1 = p1.type();                                           
        int type2 = p2.type();
        const Potential &potential = getPotential(type1, type2);

        Real3D force(0.0, 0.0, 0.0);
        if(potential._computeForce(force, p1, p2)) {
          Real3D dist = p1.position() - p2.position();
          w = w + dist * force;
        }
      }

      //adress TODO
      for (PairList::Iterator it(verletList->getAdrPairs());
                 it.isValid(); ++it) {
              Particle &p1 = *it->first;
              Particle &p2 = *it->second;
              int type1 = p1.type();
              int type2 = p2.type();
              const Potential &potential = getPotential(type1, type2);

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

    template < typename _Potential > inline void
    VerletListAdressInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor& w) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the Verlet List");

      for (PairList::Iterator it(verletList->getPairs());
           it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.type();
        int type2 = p2.type();
        const Potential &potential = getPotential(type1, type2);

        Real3D force(0.0, 0.0, 0.0);
        if(potential._computeForce(force, p1, p2)) {
          Real3D dist = p1.position() - p2.position();
          w += Tensor(dist, force);
        }
      }

      // adress TODO
      for (PairList::Iterator it(verletList->getAdrPairs());
                 it.isValid(); ++it) {
              Particle &p1 = *it->first;
              Particle &p2 = *it->second;
              int type1 = p1.type();
              int type2 = p2.type();
              const Potential &potential = getPotential(type1, type2);

              Real3D force(0.0, 0.0, 0.0);
              if(potential._computeForce(force, p1, p2)) {
                Real3D dist = p1.position() - p2.position();
                w += Tensor(dist, force);
              }
      }

    }
 
    template < typename _Potential >
    inline real
    VerletListAdressInteractionTemplate< _Potential >::
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
