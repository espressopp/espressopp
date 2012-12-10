// ESPP_CLASS 
#ifndef _INTERACTION_VERLETLISTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_VERLETLISTINTERACTIONTEMPLATE_HPP

//#include <typeinfo>


#include "types.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "VerletList.hpp"
#include "FixedTupleList.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"

namespace espresso {
  namespace interaction {
    template < typename _Potential >
    class VerletListInteractionTemplate: public Interaction {
    
    protected:
      typedef _Potential Potential;
    
    public:
      VerletListInteractionTemplate
          (shared_ptr<VerletList> _verletList)
          : verletList(_verletList) {
    	  potentialArray = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());
    	  // potentialArray = esutil::Array2D< shared_ptr<Potential>, esutil::enlarge>(0, 0, shared_ptr<Potential>());
          
        ntypes = 0;
      }

      virtual ~VerletListInteractionTemplate() {};

      void
      setVerletList(shared_ptr < VerletList > _verletList) {
        verletList = _verletList;
      }

      shared_ptr<VerletList> getVerletList() {
        return verletList;
      }

      void
      setPotential(int type1, int type2, const Potential &potential) {
       //setPotential(int type1, int type2, shared_ptr<Potential> potential) {
        // typeX+1 because i<ntypes
        ntypes = std::max(ntypes, std::max(type1+1, type2+1));
        
        potentialArray.at(type1, type2) = potential;
        if (type1 != type2) { // add potential in the other direction
           potentialArray.at(type2, type1) = potential;
        }
      }

      Potential &getPotential(int type1, int type2) {
      // shared_ptr<Potential> getPotential(int type1, int type2) {
        return potentialArray.at(type1, type2);
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
      shared_ptr<VerletList> verletList;
      esutil::Array2D<Potential, esutil::enlarge> potentialArray;
      //esutil::Array2D<shared_ptr<Potential>, esutil::enlarge> potentialArray;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _Potential > inline void
    VerletListInteractionTemplate < _Potential >::
    addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by the Verlet List");

      for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.type();
        int type2 = p2.type();
        const Potential &potential = getPotential(type1, type2);
        // shared_ptr<Potential> potential = getPotential(type1, type2);

        Real3D force(0.0, 0.0, 0.0);
        if(potential._computeForce(force, p1, p2)) {
        //if(potential->_computeForce(force, p1, p2)) {
          p1.force() += force;
          p2.force() -= force;
          LOG4ESPP_TRACE(theLogger, "id1=" << p1.id() << " id2=" << p2.id() << " force=" << force);
        }
      }
    }
    
    template < typename _Potential >
    inline real
    VerletListInteractionTemplate < _Potential >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy of the Verlet list pairs");

      real e = 0.0;
      real es = 0.0;
      for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.type();
        int type2 = p2.type();
        const Potential &potential = getPotential(type1, type2);
        // shared_ptr<Potential> potential = getPotential(type1, type2);
        e   = potential._computeEnergy(p1, p2);
        // e   = potential->_computeEnergy(p1, p2);
        es += e;
        LOG4ESPP_TRACE(theLogger, "id1=" << p1.id() << " id2=" << p2.id() << " potential energy=" << e);
      }

      // reduce over all CPUs
      real esum;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, es, esum, std::plus<real>());
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
        // shared_ptr<Potential> potential = getPotential(type1, type2);

        Real3D force(0.0, 0.0, 0.0);
        if(potential._computeForce(force, p1, p2)) {
        // if(potential->_computeForce(force, p1, p2)) {
          Real3D r21 = p1.position() - p2.position();
          w = w + r21 * force;
        }
      }

      // reduce over all CPUs
      real wsum;
      boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
      return wsum; 
    }

    template < typename _Potential > inline void
    VerletListInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor& w) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the Verlet List");

      Tensor wlocal(0.0);
      for (PairList::Iterator it(verletList->getPairs());
           it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.type();
        int type2 = p2.type();
        const Potential &potential = getPotential(type1, type2);
        // shared_ptr<Potential> potential = getPotential(type1, type2);

        Real3D force(0.0, 0.0, 0.0);
        if(potential._computeForce(force, p1, p2)) {
        // if(potential->_computeForce(force, p1, p2)) {
          Real3D r21 = p1.position() - p2.position();
          wlocal += Tensor(r21, force);
        }
      }
      
      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, wlocal, wsum, std::plus<Tensor>());
      w += wsum;
    }

    
    template < typename _Potential > inline void
    VerletListInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor& w, real z) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the Verlet List");

      System& system = verletList->getSystemRef();
      Real3D Li = system.bc->getBoxL();
      
      real rc_cutoff = verletList->getVerletCutoff();
      
      // boundaries should be taken into account
      bool ghost_layer = false;
      real zghost = -100;
      if(z<rc_cutoff){
        zghost = z + Li[2];
        ghost_layer = true;
      }
      else if(z>=Li[2]-rc_cutoff){
        zghost = z - Li[2];
        ghost_layer = true;
      }
      
      int count =0;
      Tensor wlocal(0.0);
      for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        Real3D p1pos = p1.position();
        Real3D p2pos = p2.position();

        
        if( (p1pos[2]>z && p2pos[2]<z) || 
            (p1pos[2]<z && p2pos[2]>z) ||
                (ghost_layer &&
                ((p1pos[2]>zghost && p2pos[2]<zghost) || 
                 (p1pos[2]<zghost && p2pos[2]>zghost))
                ) 
          ){
         
        /*
        if( (p1pos[2]-z>0.0 && p2pos[2]-z<0.0) || 
            (p1pos[2]-z<0.0 && p2pos[2]-z>0.0)
          ){*/
          int type1 = p1.type();
          int type2 = p2.type();
          const Potential &potential = getPotential(type1, type2);

          Real3D force(0.0, 0.0, 0.0);
          if(potential._computeForce(force, p1, p2)) {
            Real3D r21 = p1pos - p2pos;
            //wlocal += Tensor(r21, force) / r21.abs();
            Tensor ttt = Tensor(r21, force) / r21.abs();
            ttt[2] /= 2.0;
            
            wlocal += ttt;
            
            /*
            if(r21.abs()<1.0)
              std::cout<<"cpu: "<< system.comm->rank() << "  p1= "<< p1pos[2] << "  p2= "<< p2pos[2]
                      << "  r21= " << r21.abs() << "  z= " << z << "  wloc= "<< ttt[0] <<std::endl;
             **/
             
          }
          count ++;
          //if(count == 10000) exit(1);
        }
      }
      
      /*
      std::cout<<"cpu: "<< system.comm->rank() << "  count= "<< count
              << "  z= " << z << "  wloc= "<< wlocal[0] <<std::endl;
      **/
      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, wlocal, wsum, std::plus<Tensor>());
      w += wsum;
    }
    
    
    // TODO it doesn't work yet
    template < typename _Potential > inline void
    VerletListInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor *w, int n) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the Verlet List");

      System& system = verletList->getSystemRef();
      Real3D Li = system.bc->getBoxL();
      
      Tensor wlocal[n];
      for (PairList::Iterator it(verletList->getPairs());
           it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.type();
        int type2 = p2.type();
        Real3D p1pos = p1.position();
        Real3D p2pos = p2.position();
        
        int position1 = (int)( n * p1pos[2]/Li[2]);
        int position2 = (int)( n * p2pos[2]/Li[2]);
        
        int maxpos = std::max(position1, position2);
        int minpos = std::min(position1, position2); 
        
        const Potential &potential = getPotential(type1, type2);

        Real3D force(0.0, 0.0, 0.0);
        Tensor ww;
        if(potential._computeForce(force, p1, p2)) {
          Real3D r21 = p1pos - p2pos;
          ww = Tensor(r21, force);
        }
        
        int i = minpos + 1;
        while(i<=maxpos){
          wlocal[i] += ww;
          i++;
        }
      }
      
      // reduce over all CPUs
      Tensor wsum[n];
      boost::mpi::all_reduce(*mpiWorld, wlocal, n, wsum, std::plus<Tensor>());
      
      for(int j=0; j<n; j++){
        w[j] += wsum[j];
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
            // cutoff = std::max(cutoff, getPotential(i, j)->getCutoff());
        }
      }
      return cutoff;
    }
  }
}
#endif
