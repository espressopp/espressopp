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

namespace espresso {
  namespace interaction {
    template < typename _ThreeBodyPotential >
    class VerletListTripleInteractionTemplate : public Interaction, SystemAccess {
        
    protected:
      typedef _ThreeBodyPotential Potential;
      
    public:
      VerletListTripleInteractionTemplate
      (shared_ptr < System > _system,
       shared_ptr < VerletListTriple > _verletlisttriple,
       shared_ptr < Potential > _potential)
        : SystemAccess(_system), verletListTriple(_verletlisttriple),
          potential(_potential)
      {
          if (! potential) {
                LOG4ESPP_ERROR(theLogger, "NULL potential");
          }
        //potentialArray = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());
      }

      void
      setVerletListTriple(shared_ptr < VerletListTriple > _verletlisttriple) {
        verletListTriple = _verletlisttriple;
      }

      shared_ptr < VerletListTriple > getVerletListTriple() {
        return verletListTriple;
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
      virtual real computeVirial();
      virtual void computeVirialTensor(Tensor& w);
      virtual real getMaxCutoff();
      virtual int bondType() { return Angular; }

    protected:
      int ntypes;
      shared_ptr<VerletListTriple> verletListTriple;
      //esutil::Array2D<Potential, esutil::enlarge> potentialArray;
      shared_ptr < Potential > potential;
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
        //const Potential &potential = getPotential(p1.type(), p2.type());
        Real3D dist12, dist32;
        //bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
        //bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
        bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
        bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
        Real3D force12, force32;
        if(potential->_computeForce(force12, force32, dist12, dist32)){
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
      //for (VerletListTriple::TripleList::Iterator it(*verletListTriple); it.isValid(); ++it) {
      for (TripleList::Iterator it(verletListTriple->getTriples()); it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        //const Potential &potential = getPotential(p1.type(), p2.type());
        Real3D dist12 = bc.getMinimumImageVector(p1.position(), p2.position());
        Real3D dist32 = bc.getMinimumImageVector(p3.position(), p2.position());
        e += potential->_computeEnergy(dist12, dist32);
      }
      real esum;
      //boost::mpi::reduce(*mpiWorld, e, esum, std::plus<real>(), 0);
      boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
      return esum;
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
        const espresso::bc::BC& bc = *getSystemRef().bc;
        Real3D dist12, dist32;
        bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
        bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
        Real3D force12, force32;
        if(potential->_computeForce(force12, force32, dist12, dist32)){
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

      const bc::BC& bc = *getSystemRef().bc;
      for (TripleList::Iterator it(verletListTriple->getTriples()); it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        //const Potential &potential = getPotential(0, 0);
        Real3D dist12, dist32;
        bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
        bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
        Real3D force12, force32;
        if(potential->_computeForce(force12, force32, dist12, dist32)){
          w += Tensor(dist12, force12) + Tensor(dist32, force32);
        }
      }
    }

    template < typename _ThreeBodyPotential >
    inline real
    VerletListTripleInteractionTemplate< _ThreeBodyPotential >::
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
