/*
  Copyright (C) 2012,2013,2014,2015,2016,2017,2018
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI

  This file is part of ESPResSo++.

  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// ESPP_CLASS
#ifndef _INTERACTION_CELLLISTALLPAIRSINTERACTIONTEMPLATE_HPP
#define _INTERACTION_CELLLISTALLPAIRSINTERACTIONTEMPLATE_HPP

#include "types.hpp"
#include "Tensor.hpp"
#include "Interaction.hpp"
#include "storage/Storage.hpp"
#include "esutil/Array2D.hpp"
#include "iterator/CellListAllPairsIterator.hpp"

namespace espressopp {
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
      virtual real computeEnergyDeriv();
      virtual real computeEnergyAA();
      virtual real computeEnergyCG();
      virtual real computeEnergyAA(int atomtype);
      virtual real computeEnergyCG(int atomtype);
      virtual void computeVirialX(std::vector<real> &p_xx_total, int bins);
      virtual real computeVirial();
      virtual void computeVirialTensor(Tensor& wij);
      virtual void computeVirialTensor(Tensor& w, real z);
      virtual void computeVirialTensor(Tensor *w, int n);
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
    computeEnergyDeriv() {
      std::cout << "Warning! At the moment computeEnergyDeriv() in CellListAllPairsInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline real
    CellListAllPairsInteractionTemplate < _Potential >::
    computeEnergyAA() {
      std::cout << "Warning! At the moment computeEnergyAA() in CellListAllPairsInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline real
    CellListAllPairsInteractionTemplate < _Potential >::
    computeEnergyAA(int atomtype) {
      std::cout << "Warning! At the moment computeEnergyAA(int atomtype) in CellListAllPairsInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline real
    CellListAllPairsInteractionTemplate < _Potential >::
    computeEnergyCG() {
      std::cout << "Warning! At the moment computeEnergyCG() in CellListAllPairsInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline real
    CellListAllPairsInteractionTemplate < _Potential >::
    computeEnergyCG(int atomtype) {
      std::cout << "Warning! At the moment computeEnergyCG(int atomtype) in CellListAllPairsInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential >
    inline void
    CellListAllPairsInteractionTemplate < _Potential >::
    computeVirialX(std::vector<real> &p_xx_total, int bins) {
        std::cout << "Warning! At the moment computeVirialX in CellListAllPairsInteractionTemplate does not work." << std::endl << "Therefore, the corresponding interactions won't be included in calculation." << std::endl;
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
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
      wij += wsum;
    }

    template < typename _Potential > inline void
    CellListAllPairsInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor& wij, real z) {
      LOG4ESPP_INFO(theLogger, "computed virial tensor for all pairs in the cell lists");

      Tensor wlocal(0.0);
      const bc::BC& bc = *storage->getSystemRef().bc;  // boundary conditions
      for (iterator::CellListAllPairsIterator it(storage->getRealCells());
           it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        Real3D p1pos = p1.position();
        Real3D p2pos = p2.position();

        if(  (p1pos[2]>=z && p2pos[2]<=z) ||
             (p1pos[2]<=z && p2pos[2]>=z) ){
          const Potential &potential = getPotential(p1.type(), p2.type());

          Real3D force(0.0, 0.0, 0.0);
          if(potential._computeForce(force, p1, p2)) {
            Real3D r21;
            bc.getMinimumImageVectorBox(r21, p1pos, p2pos);
            wlocal += Tensor(r21, force);
          }
        }
      }

      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
      wij += wsum;
    }


    template < typename _Potential > inline void
    CellListAllPairsInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor *wij, int n) {
      LOG4ESPP_INFO(theLogger, "computed virial tensor for all pairs in the cell lists");

      const bc::BC& bc = *storage->getSystemRef().bc;  // boundary conditions
      Real3D Li = bc.getBoxL();
      Tensor *wlocal = new Tensor[n];
      for(int i=0; i<n; i++) wlocal[i] = Tensor(0.0);
      for (iterator::CellListAllPairsIterator it(storage->getRealCells());
           it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        Real3D p1pos = p1.position();
        Real3D p2pos = p2.position();

        int position1 = (int)( n * p1pos[2]/Li[2]);
        int position2 = (int)( n * p1pos[2]/Li[2]);

        int maxpos = std::max(position1, position2);
        int minpos = std::min(position1, position2);

        const Potential &potential = getPotential(p1.type(), p2.type());

        Real3D force(0.0, 0.0, 0.0);
        Tensor ww;
        if(potential._computeForce(force, p1, p2)) {
          Real3D r21;
          bc.getMinimumImageVectorBox(r21, p1pos, p2pos);
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
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, n, (double*)&wsum, std::plus<double>());

      for(int j=0; j<n; j++){
        wij[j] += wsum[j];
      }

      delete [] wsum;
      delete [] wlocal;
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
