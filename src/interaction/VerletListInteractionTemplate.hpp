// ESPP_CLASS 
#ifndef _INTERACTION_VERLETLISTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_VERLETLISTINTERACTIONTEMPLATE_HPP

#include "types.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "VerletList.hpp"
#include "esutil/Array2D.hpp"

namespace espresso {
  namespace interaction {
    template < typename _Potential >
    class VerletListInteractionTemplate: public Interaction {
    
    protected:
      typedef _Potential Potential;
    
    public:
      VerletListInteractionTemplate (shared_ptr < VerletList > _verletList)
        : verletList(_verletList) {
        potentialArray = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());

        pidhy2 = M_PI/(verletList->getHy() * 2);
        dex = verletList->getEx();
        dexdhy = dex + verletList->getHy();
      }

      void
      setVerletList(shared_ptr < VerletList > _verletList) {
        verletList = _verletList;
      }

      shared_ptr < VerletList > getVerletList() {
        return verletList;
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
      shared_ptr < VerletList > verletList;
      esutil::Array2D<Potential, esutil::enlarge> potentialArray;

      // adress stuff
      real pidhy2; // pi / (dhy * 2)
      real dexdhy; // dhy + dex
      real dex;

    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _Potential > inline void
    VerletListInteractionTemplate < _Potential >::
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
        }
      }

      // adress TODO
      //std::cout << "add forces computed by the AdResS List" << "\n";
      for (PairList::Iterator it(verletList->getAdrPairs());
         it.isValid(); ++it) {
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
         //else if (dexdhy < min1) w1 = 0;
         else {
             w1 = cos(pidhy2 * (min1 - dex));
             w1 *= w1;
         }

         if (dex > min2) w2 = 1;
         //else if (dexdhy < min2) w2 = 0;
         else {
             w2 = cos(pidhy2 * (min2 - dex));
             w2 *= w2;
         }

         //if (w1 == 1 || w2 == 1) std::cout << p1.id() << " ";
         //std::cout << p1.id() << ": " << w1 << "\n";
         //std::cout << p2.id() << ": " << w2 << "\n";

         Real3D force(0.0, 0.0, 0.0);
         if(potential._computeForce(force, p1, p2)) {
            p1.force() += force;
            p2.force() -= force;
         }
      }
      std::cout << "\n\n";

    }
    
    template < typename _Potential >
    inline real
    VerletListInteractionTemplate < _Potential >::
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
    VerletListInteractionTemplate < _Potential >::
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
    VerletListInteractionTemplate < _Potential >::
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
    VerletListInteractionTemplate< _Potential >::
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
