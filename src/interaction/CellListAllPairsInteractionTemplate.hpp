// ESPP_CLASS
#ifndef _INTERACTION_CELLLISTALLPAIRSINTERACTIONTEMPLATE_HPP
#define _INTERACTION_CELLLISTALLPAIRSINTERACTIONTEMPLATE_HPP

#include "types.hpp"
#include "Tensor.hpp"
#include "Interaction.hpp"
#include "storage/Storage.hpp"
#include "esutil/Array2D.hpp"
#include "iterator/CellListAllPairsIterator.hpp"

namespace espresso {
  namespace interaction {
    template < typename _Potential >
    class CellListAllPairsInteractionTemplate
        : public Interaction {
    protected:
      typedef _Potential Potential;
    public:
      CellListAllPairsInteractionTemplate
      (shared_ptr < storage::Storage > _storage)
      : storage(_storage){
        potentialArray = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());
        ntypes=0;
      }

      void
      setPotential(int type1, int type2, const Potential &potential) {
        // typeX+1 because i<ntypes
        ntypes = std::max(ntypes, std::max(type1+1, type2+1));
        
        potentialArray.at(type1, type2) = potential;
        
        // @TODO should it be here???
//        if (type1 != type2) { // add potential in the other direction
//           potentialArray.at(type2, type1) = potential;
//        }
        
      }

      Potential &getPotential(int type1, int type2) {
        return potentialArray(type1, type2);
      }

      virtual void addForces();
      virtual real computeEnergy();
      virtual real computeVirial();
      virtual void computeVirialTensor(Tensor& wij);
      virtual void computeVirialTensor(Tensor& w, real xmin, real xmax,
          real ymin, real ymax, real zmin, real zmax);
      virtual real getMaxCutoff();
      virtual int bondType() { return Nonbonded; }

    protected:
      int ntypes;
      esutil::Array2D< Potential, esutil::enlarge > potentialArray;
      shared_ptr< storage::Storage > storage;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _Potential > inline void
    CellListAllPairsInteractionTemplate < _Potential >::
    addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed for all pairs in the cell lists");

      for (iterator::CellListAllPairsIterator it(storage->getRealCells()); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        const Potential &potential = getPotential(p1.type(), p2.type());

        Real3D force(0.0, 0.0, 0.0);
        if(potential._computeForce(force, p1, p2)) {
          p1.force() += force;
          p2.force() -= force;
        }
      }
    }

    template < typename _Potential > inline real
    CellListAllPairsInteractionTemplate < _Potential >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy by the Verlet List");

      real e = 0.0;
      for (iterator::CellListAllPairsIterator it(storage->getRealCells()); it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Potential &potential = getPotential(p1.type(), p2.type());
        e += potential._computeEnergy(p1, p2);
      }

      // reduce over all CPUs
      real esum;
      boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
      return esum;
    }

    template < typename _Potential > inline real 
    CellListAllPairsInteractionTemplate < _Potential >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "computed virial for all pairs in the cell lists");
     
      real w = 0.0;
      for (iterator::CellListAllPairsIterator it(storage->getRealCells());
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
      
      // reduce over all CPUs
      real wsum;
      boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
      return wsum; 
    }

    template < typename _Potential > inline void
    CellListAllPairsInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor& wij) {
      LOG4ESPP_INFO(theLogger, "computed virial tensor for all pairs in the cell lists");

      Tensor wlocal(0.0);
      for (iterator::CellListAllPairsIterator it(storage->getRealCells());
           it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Potential &potential = getPotential(p1.type(), p2.type());

        Real3D force(0.0, 0.0, 0.0);
        if(potential._computeForce(force, p1, p2)) {
          Real3D dist = p1.position() - p2.position();
          wlocal += Tensor(dist, force);
        }
      }
      
      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, wlocal, wsum, std::plus<Tensor>());
      wij += wsum;
    }
    

    // compute the pressure tensor localized between xmin, xmax, ymin, ymax, zmin, zmax
    // TODO (vit) physics should be checked
    template < typename _Potential > inline void
    CellListAllPairsInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor& wij,
            real xmin, real xmax, real ymin, real ymax, real zmin, real zmax) {
      LOG4ESPP_INFO(theLogger, "computed virial tensor for all pairs in the cell lists");

      Tensor wlocal(0.0);
      for (iterator::CellListAllPairsIterator it(storage->getRealCells());
           it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        Real3D p1pos = p1.position();
        Real3D p2pos = p2.position();
        
        if(  (p1pos[0]>xmin && p1pos[0]<xmax && 
              p1pos[1]>ymin && p1pos[1]<ymax && 
              p1pos[2]>zmin && p1pos[2]<zmax) ||
             (p2pos[0]>xmin && p2pos[0]<xmax && 
              p2pos[1]>ymin && p2pos[1]<ymax && 
              p2pos[2]>zmin && p2pos[2]<zmax) ){
          const Potential &potential = getPotential(p1.type(), p2.type());

          Real3D force(0.0, 0.0, 0.0);
          if(potential._computeForce(force, p1, p2)) {
            Real3D dist = p1pos - p2pos;
            wlocal += Tensor(dist, force);
          }
        }
      }
      
      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, wlocal, wsum, std::plus<Tensor>());
      wij += wsum;
    }

    template < typename _Potential >
    inline real
    CellListAllPairsInteractionTemplate < _Potential >::getMaxCutoff() {
      real cutoff = 0.0;
      for(int i = 0; i < ntypes; i++) {
        for(int j = 0; j < ntypes; j++) {
          cutoff = std::max(cutoff, getPotential(i, j).getCutoff());
        }
      }
      return cutoff;
    }
  }
}


#endif
