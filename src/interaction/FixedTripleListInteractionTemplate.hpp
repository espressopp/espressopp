// ESPP_CLASS
#ifndef _INTERACTION_FIXEDTRIPLELISTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDTRIPLELISTINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "types.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "FixedTripleList.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"

namespace espresso {
  namespace interaction {
    template < typename _AngularPotential >
    class FixedTripleListInteractionTemplate : public Interaction, SystemAccess {
        
    protected:
      typedef _AngularPotential Potential;
      
    public:
      FixedTripleListInteractionTemplate
      (shared_ptr < System > _system,
       shared_ptr < FixedTripleList > _fixedtripleList)
        : SystemAccess(_system), fixedtripleList(_fixedtripleList) 
      {
        potentialArray = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());
      }

      void
      setFixedTripleList(shared_ptr < FixedTripleList > _fixedtripleList) {
        fixedtripleList = _fixedtripleList;
      }

      shared_ptr < FixedTripleList > getFixedTripleList() {
        return fixedtripleList;
      }

      void
      setPotential(int type1, int type2, const Potential &potential) {
        potentialArray.at(type1, type2) = potential;
      }

      Potential &getPotential(int type1, int type2) {
        return potentialArray.at(0, 0);
      }

      virtual void addForces();
      virtual real computeEnergy();
      virtual real computeVirial();
      virtual void computeVirialTensor(Tensor& w);
      virtual real getMaxCutoff();

    protected:
      int ntypes;
      shared_ptr < FixedTripleList > fixedtripleList;
      esutil::Array2D<Potential, esutil::enlarge> potentialArray;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _AngularPotential > inline void
    FixedTripleListInteractionTemplate < _AngularPotential >::
    addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by FixedTripleList");
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedTripleList::Iterator it(*fixedtripleList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        Particle &p3 = *it->third;
        const Potential &potential = getPotential(p1.type(), p2.type());
        Real3D dist12, dist32;
        bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
        bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
        Real3D force12, force32;
        potential._computeForce(force12, force32, dist12, dist32);
        p1.force() += force12;
        p2.force() -= force12 + force32;
        p3.force() += force32;
      }
    }

    template < typename _AngularPotential >
    inline real
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy of the triples");

      const bc::BC& bc = *getSystemRef().bc;
      real e = 0.0;
      for (FixedTripleList::Iterator it(*fixedtripleList); it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        const Potential &potential = getPotential(p1.type(), p2.type());
        Real3D dist12 = bc.getMinimumImageVector(p1.position(), p2.position());
        Real3D dist32 = bc.getMinimumImageVector(p3.position(), p2.position());
        e += potential._computeEnergy(dist12, dist32);
      }
      real esum;
      boost::mpi::reduce(*mpiWorld, e, esum, std::plus<real>(), 0);
      return esum;
    }

    template < typename _AngularPotential >
    inline real
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute scalar virial of the triples");

      const bc::BC& bc = *getSystemRef().bc;
      real w = 0.0;
      for (FixedTripleList::Iterator it(*fixedtripleList); it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        const Potential &potential = getPotential(p1.type(), p2.type());
        const espresso::bc::BC& bc = *getSystemRef().bc;
        Real3D dist12, dist32;
        bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
        bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
        Real3D force12, force32;
        potential._computeForce(force12, force32, dist12, dist32);
        w += dist12 * force12 + dist32 * force32;
      }
      return w;
    }

    template < typename _AngularPotential >
    inline void
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeVirialTensor(Tensor& w) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor of the triples");

      const bc::BC& bc = *getSystemRef().bc;
      for (FixedTripleList::Iterator it(*fixedtripleList); it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        const Potential &potential = getPotential(0, 0);
        Real3D dist12, dist32;
        bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
        bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
        Real3D force12, force32;
        potential._computeForce(force12, force32, dist12, dist32);
        w += Tensor(dist12, force12) + Tensor(dist32, force32);
      }
    }

    template < typename _AngularPotential >
    inline real
    FixedTripleListInteractionTemplate< _AngularPotential >::
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
