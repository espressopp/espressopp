// ESPP_CLASS
#ifndef _INTERACTION_FIXEDPAIRDISTLISTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDPAIRDISTLISTINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "FixedPairDistList.hpp"
#include "FixedPairListAdress.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "types.hpp"

namespace espresso {
  namespace interaction {
    template < typename _Potential >
    class FixedPairDistListInteractionTemplate: public Interaction, SystemAccess {
        
    protected:
      typedef _Potential Potential;
      
    public:
      FixedPairDistListInteractionTemplate
      (shared_ptr < System > system,
       shared_ptr < FixedPairDistList > _fixedPairDistList,
       shared_ptr < Potential > _potential)
        : SystemAccess(system), fixedPairDistList(_fixedPairDistList),
          potential(_potential)
      {
        if (! potential) {
          LOG4ESPP_ERROR(theLogger, "NULL potential");
        }
      }

      virtual ~FixedPairDistListInteractionTemplate() {};

      void
      setFixedPairList(shared_ptr < FixedPairDistList > _fixedPairDistList) {
        fixedPairDistList = _fixedPairDistList;
      }

      shared_ptr < FixedPairDistList > getFixedPairList() {
        return fixedPairDistList;
      }

      void
      setPotential(shared_ptr < Potential> _potential) {
        if (_potential) {
          potential = _potential;
        } else {
          LOG4ESPP_ERROR(theLogger, "NULL potential");
        }
      }

      shared_ptr < Potential > getPotential() {
        return potential;
      }

      virtual void addForces();
      virtual real computeEnergy();
      virtual void computeVirialX(std::vector<real> &p_xx_total, int bins); 
      virtual real computeVirial();
      virtual void computeVirialTensor(Tensor& w);
      virtual void computeVirialTensor(Tensor& w, real z);
      virtual void computeVirialTensor(Tensor *w, int n);
      virtual real getMaxCutoff();
      virtual int bondType() { return Pair; }

    protected:
      int ntypes;
      shared_ptr < FixedPairDistList > fixedPairDistList;
      shared_ptr < Potential > potential;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _Potential > inline void
    FixedPairDistListInteractionTemplate < _Potential >::addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by the FixedPair List");
      const bc::BC& bc = *getSystemRef().bc;
      for (FixedPairDistList::PairList::Iterator it(*fixedPairDistList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        Real3D r21(0.0,0.0,0.0);
        bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
        Real3D force;
        
        real currentDist = fixedPairDistList->getDist(p1.getId(), p2.getId());
        
        if(potential->_computeForce(force, r21, currentDist)) {
          p1.force() += force;
          p2.force() -= force;
        }
      }
    }
    
    template < typename _Potential > inline real
    FixedPairDistListInteractionTemplate < _Potential >::
    computeEnergy() {

      LOG4ESPP_INFO(theLogger, "compute energy of the FixedPair list pairs");

      real e = 0.0;
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedPairDistList::PairList::Iterator it(*fixedPairDistList);
	   it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        Real3D r21(0.0,0.0,0.0);
        bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
        real currentDist = fixedPairDistList->getDist(p1.getId(), p2.getId());
        
        e += potential->_computeEnergy(r21, currentDist);
      }
      real esum;
      boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
      return esum;
    }
    
    template < typename _Potential >
    inline void
    FixedPairDistListInteractionTemplate < _Potential >::
    computeVirialX(std::vector<real> &p_xx_total, int bins) {
        std::cout << "Warning! At the moment computeVirialX in FixedPairDistListInteractionTemplate does not work." << std::endl << "Therefore, the corresponding interactions won't be included in calculation." << std::endl;
    }   

    template < typename _Potential > inline real
    FixedPairDistListInteractionTemplate < _Potential >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute the virial for the FixedPair List");
      
      real w = 0.0;
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedPairDistList::PairList::Iterator it(*fixedPairDistList);
           it.isValid(); ++it) {                                         
        const Particle &p1 = *it->first;                                       
        const Particle &p2 = *it->second;                                      

        Real3D r21(0.0,0.0,0.0);
        bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
        real currentDist = fixedPairDistList->getDist(p1.getId(), p2.getId());
        
        Real3D force;
        if(potential->_computeForce(force, r21, currentDist)) {
          w += r21 * force;
        }
      }
      
      real wsum;
      boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
      return wsum;
    }

    template < typename _Potential > inline void
    FixedPairDistListInteractionTemplate < _Potential >::computeVirialTensor(Tensor& w){
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");

      Tensor wlocal(0.0);
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedPairDistList::PairList::Iterator it(*fixedPairDistList);
           it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        Real3D r21(0.0,0.0,0.0);
        bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
        real currentDist = fixedPairDistList->getDist(p1.getId(), p2.getId());
        
        Real3D force;
        if(potential->_computeForce(force, r21, currentDist)) { 
          wlocal += Tensor(r21, force);
        }
      }
      
      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, wlocal, wsum, std::plus<Tensor>());
      w += wsum;
    }

    template < typename _Potential > inline void
    FixedPairDistListInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor& w, real z){
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");

      Tensor wlocal(0.0);
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedPairDistList::PairList::Iterator it(*fixedPairDistList);
           it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        Real3D p1pos = p1.position();
        Real3D p2pos = p2.position();
        
        if(  (p1pos[2]>=z && p2pos[2]<=z) ||
             (p1pos[2]<=z && p2pos[2]>=z) ){
          Real3D r21(0.0,0.0,0.0);
          bc.getMinimumImageVectorBox(r21, p1pos, p2pos);
          real currentDist = fixedPairDistList->getDist(p1.getId(), p2.getId());
        
          Real3D force;
          if(potential->_computeForce(force, r21, currentDist)) { 
            wlocal += Tensor(r21, force);
          }
        }
      }
      
      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, wlocal, wsum, std::plus<Tensor>());
      w += wsum;
    }
    
    template < typename _Potential > inline void
    FixedPairDistListInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor *w, int n){
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");

      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      Real3D Li = bc.getBoxL();
      Tensor *wlocal = new Tensor[n];
      for(int i=0; i<n; i++) wlocal[i] = Tensor(0.0);
      for (FixedPairDistList::PairList::Iterator it(*fixedPairDistList);
           it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        Real3D p1pos = p1.position();
        Real3D p2pos = p2.position();
        
        int position1 = (int)( n * p1pos[2]/Li[2]);
        int position2 = (int)( n * p1pos[2]/Li[2]);
        
        int maxpos = std::max(position1, position2);
        int minpos = std::min(position1, position2); 
        
        Real3D r21(0.0,0.0,0.0);
        bc.getMinimumImageVectorBox(r21, p1pos, p2pos);
        real currentDist = fixedPairDistList->getDist(p1.getId(), p2.getId());

        Real3D force;
        Tensor ww;
        if(potential->_computeForce(force, r21, currentDist)) { 
          ww = Tensor(r21, force);
        }
        
        int i = minpos + 1;
        while(i<=maxpos){
          wlocal[i] += ww;
          i++;
        }
      }
      
      // reduce over all CPUs
      Tensor *wsum = new Tensor[n];
      boost::mpi::all_reduce(*mpiWorld, wlocal, n, wsum, std::plus<Tensor>());
      
      for(int j=0; j<n; j++){
        w[j] += wsum[j];
      }

      delete [] wsum;
      delete [] wlocal;
    }
    
    template < typename _Potential >
    inline real
    FixedPairDistListInteractionTemplate< _Potential >::
    getMaxCutoff() {
      return potential->getCutoff();
    }
  }
}
#endif
