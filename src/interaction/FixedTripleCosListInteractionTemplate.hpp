// ESPP_CLASS
#ifndef _INTERACTION_FIXEDTRIPLECOSLISTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDTRIPLECOSLISTINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "types.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "FixedTripleCosList.hpp"
//#include "FixedTripleListAdress.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"

/*
 * Temporary interaction template for angular potential based on fixed triple list
 * with fixed unique angles (cosines)
 */

namespace espresso {
  namespace interaction {
    template < typename _AngularPotential >
    class FixedTripleCosListInteractionTemplate : public Interaction, SystemAccess {
        
    protected:
      typedef _AngularPotential Potential;
      
    public:
      FixedTripleCosListInteractionTemplate
      (shared_ptr < System > _system,
       shared_ptr < FixedTripleCosList > _fixedtripleList,
       shared_ptr < Potential > _potential)
        : SystemAccess(_system), fixedtripleList(_fixedtripleList), potential(_potential)
      {
        if(! potential){
          LOG4ESPP_ERROR(theLogger, "NULL potential");
        }
      }

      virtual ~FixedTripleCosListInteractionTemplate() {};

      void
      setFixedTripleCosList(shared_ptr < FixedTripleCosList > _fixedtripleList) {
        fixedtripleList = _fixedtripleList;
      }

      shared_ptr < FixedTripleCosList > 
      getFixedTripleList() {
        return fixedtripleList;
      }

      void
      setPotential(shared_ptr < Potential> _potential) {
        if (_potential) {
          potential = _potential;
        }
        else{
          LOG4ESPP_ERROR(theLogger, "NULL potential");
        }
      }

      shared_ptr < Potential >
      getPotential() {
        return potential;
      }

      virtual void addForces();
      virtual real computeEnergy();
      virtual real computeVirial();
      virtual void computeVirialTensor(Tensor& w);
      virtual void computeVirialTensor(Tensor& w, real xmin, real xmax,
          real ymin, real ymax, real zmin, real zmax);
      virtual real getMaxCutoff();
      virtual int bondType() { return Angular; }

    protected:
      int ntypes;
      shared_ptr<FixedTripleCosList> fixedtripleList;
      shared_ptr < Potential > potential;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _AngularPotential > inline void
    FixedTripleCosListInteractionTemplate <_AngularPotential>::
    addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by FixedTripleCosList");
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedTripleCosList::TripleList::Iterator it(*fixedtripleList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        Particle &p3 = *it->third;
        //const Potential &potential = getPotential(p1.type(), p2.type());
        Real3D dist12, dist32;
        bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
        bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
        Real3D force12, force32;
        
        real currentCos = fixedtripleList->getCos(p1.getId(), p2.getId(), p3.getId());
        
        potential->_computeForce(force12, force32, dist12, dist32, currentCos);
        p1.force() += force12;
        p2.force() -= force12 + force32;
        p3.force() += force32;
      }
    }

    template < typename _AngularPotential > inline real
    FixedTripleCosListInteractionTemplate < _AngularPotential >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy of the triples");

      const bc::BC& bc = *getSystemRef().bc;
      real e = 0.0;
      for (FixedTripleCosList::TripleList::Iterator it(*fixedtripleList); it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        //const Potential &potential = getPotential(p1.type(), p2.type());
        Real3D dist12 = bc.getMinimumImageVector(p1.position(), p2.position());
        Real3D dist32 = bc.getMinimumImageVector(p3.position(), p2.position());
        
        real currentCos = fixedtripleList->getCos(p1.getId(), p2.getId(), p3.getId());
        
        e += potential->_computeEnergy(dist12, dist32, currentCos);
      }
      real esum;
      boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
      return esum;
    }

    template < typename _AngularPotential > inline real
    FixedTripleCosListInteractionTemplate < _AngularPotential >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute scalar virial of the triples");

      const bc::BC& bc = *getSystemRef().bc;
      real w = 0.0;
      for (FixedTripleCosList::TripleList::Iterator it(*fixedtripleList); it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        //const Potential &potential = getPotential(p1.type(), p2.type());
        const espresso::bc::BC& bc = *getSystemRef().bc;
        Real3D dist12, dist32;
        bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
        bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
        Real3D force12, force32;

        real currentCos = fixedtripleList->getCos(p1.getId(), p2.getId(), p3.getId());
        
        potential->_computeForce(force12, force32, dist12, dist32, currentCos);
        w += dist12 * force12 + dist32 * force32;
      }
      real wsum;
      boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
      return wsum;
    }

    template < typename _AngularPotential > inline void
    FixedTripleCosListInteractionTemplate < _AngularPotential >::
    computeVirialTensor(Tensor& w) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor of the triples");

      Tensor wlocal(0.0);
      const bc::BC& bc = *getSystemRef().bc;
      for (FixedTripleCosList::TripleList::Iterator it(*fixedtripleList); it.isValid(); ++it){
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        //const Potential &potential = getPotential(0, 0);
        Real3D r12, r32;
        bc.getMinimumImageVectorBox(r12, p1.position(), p2.position());
        bc.getMinimumImageVectorBox(r32, p3.position(), p2.position());
        Real3D force12, force32;

        real currentCos = fixedtripleList->getCos(p1.getId(), p2.getId(), p3.getId());
        
        potential->_computeForce(force12, force32, r12, r32, currentCos);
        wlocal += Tensor(r12, force12) + Tensor(r32, force32);
      }
      
      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, wlocal, wsum, std::plus<Tensor>());
      w += wsum;
    }

    // compute the pressure tensor localized between xmin, xmax, ymin, ymax, zmin, zmax
    // TODO physics should be checked
    // now it calculates if the 2 particle is inside the volume
    template < typename _AngularPotential > inline void
    FixedTripleCosListInteractionTemplate < _AngularPotential >::
    computeVirialTensor(Tensor& w,
            real xmin, real xmax, real ymin, real ymax, real zmin, real zmax) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor of the triples");

      Tensor wlocal(0.0);
      const bc::BC& bc = *getSystemRef().bc;
      for (FixedTripleCosList::TripleList::Iterator it(*fixedtripleList); it.isValid(); ++it){
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

          real currentCos = fixedtripleList->getCos(p1.getId(), p2.getId(), p3.getId());
        
          potential->_computeForce(force12, force32, r12, r32, currentCos);
          wlocal += Tensor(r12, force12) + Tensor(r32, force32);
        }
      }
      
      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, wlocal, wsum, std::plus<Tensor>());
      w += wsum;
    }
    
    template < typename _AngularPotential >
    inline real
    FixedTripleCosListInteractionTemplate< _AngularPotential >::
    getMaxCutoff() {
      return potential->getCutoff();
    }
  }
}
#endif
