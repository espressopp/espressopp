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
#ifndef _INTERACTION_VERLETLISTVSPHEREINTERACTIONTEMPLATE_HPP
#define _INTERACTION_VERLETLISTVSPHEREINTERACTIONTEMPLATE_HPP

//#include <typeinfo>


#include "types.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "VerletList.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"

#include "storage/Storage.hpp"

namespace espressopp {
  namespace interaction {
    template < typename _Potential >
    class VerletListVSphereInteractionTemplate: public Interaction {

    protected:
      typedef _Potential Potential;

    public:
      VerletListVSphereInteractionTemplate
          (shared_ptr<VerletList> _verletList)
          : verletList(_verletList) {
    	  potentialArray    = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());
        ntypes = 0;
      }

      virtual ~VerletListVSphereInteractionTemplate() {};

      void
      setVerletList(shared_ptr < VerletList > _verletList) {
        verletList = _verletList;
      }

      shared_ptr<VerletList> getVerletList() {
        return verletList;
      }

      void
      setPotential(int type1, int type2, const Potential &potential) {
        // typeX+1 because i<ntypes
        ntypes = std::max(ntypes, std::max(type1+1, type2+1));
        potentialArray.at(type1, type2) = potential;
        if (type1 != type2) { // add potential in the other direction
           potentialArray.at(type2, type1) = potential;
        }
      }

      // this is used in the innermost force-loop
      Potential &getPotential(int type1, int type2) {
        return potentialArray.at(type1, type2);
      }

      // this is mainly used to access the potential from Python (e.g. to change parameters of the potential)
      shared_ptr<Potential> getPotentialPtr(int type1, int type2) {
    	return  make_shared<Potential>(potentialArray.at(type1, type2));
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
      virtual void computeVirialTensor(Tensor& w);
      virtual void computeVirialTensor(Tensor& w, real z);
      virtual void computeVirialTensor(Tensor *w, int n);
      virtual real getMaxCutoff();
      virtual int bondType() { return Nonbonded; }

    protected:
      int ntypes;
      shared_ptr<VerletList> verletList;
      esutil::Array2D<Potential, esutil::enlarge> potentialArray;
      esutil::Array2D<shared_ptr<Potential>, esutil::enlarge> potentialArrayPtr;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _Potential > inline void
    VerletListVSphereInteractionTemplate < _Potential >::
    addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by the Verlet List");

      for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.type();
        int type2 = p2.type();
        const Potential &potential = getPotential(type1, type2);
        // shared_ptr<Potential> potential = getPotential(type1, type2);

        Real3D force(0.0);
        real fsi=0.0, fsj=0.0;
        if(potential._computeForce(force, fsi, fsj, p1, p2)) {
          p1.force() += force;
          p2.force() -= force;

          LOG4ESPP_TRACE(theLogger, "id1=" << p1.id() << " id2=" << p2.id() << " force=" << force);
        }
      }
    }

    template < typename _Potential >
    inline real
    VerletListVSphereInteractionTemplate < _Potential >::
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
    VerletListVSphereInteractionTemplate < _Potential >::
    computeEnergyDeriv() {
      std::cout << "Warning! At the moment computeEnergyDeriv() in VerletListVSphereInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline real
    VerletListVSphereInteractionTemplate < _Potential >::
    computeEnergyAA() {
      std::cout << "Warning! At the moment computeEnergyAA() in VerletListVSphereInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline real
    VerletListVSphereInteractionTemplate < _Potential >::
    computeEnergyAA(int atomtype) {
      std::cout << "Warning! At the moment computeEnergyAA(int atomtype) in VerletListVSphereInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline real
    VerletListVSphereInteractionTemplate < _Potential >::
    computeEnergyCG() {
      std::cout << "Warning! At the moment computeEnergyCG() in VerletListVSphereInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline real
    VerletListVSphereInteractionTemplate < _Potential >::
    computeEnergyCG(int atomtype) {
      std::cout << "Warning! At the moment computeEnergyCG(int atomtype) in VerletListVSphereInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential >
    inline void
    VerletListVSphereInteractionTemplate < _Potential >::
    computeVirialX(std::vector<real> &p_xx_total, int bins) {
        std::cout << "Warning! At the moment computeVirialX in VerletListInteractionTemplate does not work." << std::endl << "Therefore, the corresponding interactions won't be included in calculation." << std::endl;
    }

    template < typename _Potential > inline real
    VerletListVSphereInteractionTemplate < _Potential >::
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
        real fsi, fsj;
        if(potential._computeForce(force, fsi, fsj, p1, p2)) {
        // if(potential->_computeForce(force, p1, p2)) {
        //TODO think of how to incorporate sigmaij-force into virial calculation
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
    VerletListVSphereInteractionTemplate < _Potential >::
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
        real fsi, fsj;
        if(potential._computeForce(force, fsi, fsj, p1, p2)) {
        // if(potential->_computeForce(force, p1, p2)) {
        //TODO think of how to incorporate sigmaij-force into virial calculation
          Real3D r21 = p1.position() - p2.position();
          wlocal += Tensor(r21, force);
        }
      }

      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
      w += wsum;
    }

    // local pressure tensor for layer, plane is defined by z coordinate
    template < typename _Potential > inline void
    VerletListVSphereInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor& w, real z) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the Verlet List");

      System& system = verletList->getSystemRef();
      Real3D Li = system.bc->getBoxL();

      real rc_cutoff = verletList->getVerletCutoff();

      // boundaries should be taken into account
      bool ghost_layer = false;
      real zghost = -100.0;
      if(z<rc_cutoff){
        zghost = z + Li[2];
        ghost_layer = true;
      }
      else if(z>=Li[2]-rc_cutoff){
        zghost = z - Li[2];
        ghost_layer = true;
      }

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
          int type1 = p1.type();
          int type2 = p2.type();
          const Potential &potential = getPotential(type1, type2);

          Real3D force(0.0, 0.0, 0.0);
          real fsi, fsj;
          if(potential._computeForce(force, fsi, fsj, p1, p2)) {
          //TODO think of how to incorporate sigmaij-force into virial calculation
            Real3D r21 = p1pos - p2pos;
            wlocal += Tensor(r21, force) / fabs(r21[2]);
          }
        }
      }

      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
      w += wsum;
    }

    // it will calculate the pressure in 'n' layers along Z axis
    // the first layer has coordinate 0.0 the last - (Lz - Lz/n)
    template < typename _Potential > inline void
    VerletListVSphereInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor *w, int n) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the Verlet List");

      System& system = verletList->getSystemRef();
      Real3D Li = system.bc->getBoxL();

      real z_dist = Li[2] / float(n);  // distance between two layers
      Tensor *wlocal = new Tensor[n];
      for(int i=0; i<n; i++) wlocal[i] = Tensor(0.0);
      for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.type();
        int type2 = p2.type();
        Real3D p1pos = p1.position();
        Real3D p2pos = p2.position();

        const Potential &potential = getPotential(type1, type2);

        Real3D force(0.0, 0.0, 0.0);
        Tensor ww;
        real fsi, fsj;
        if(potential._computeForce(force, fsi, fsj, p1, p2)) {
          //TODO think of how to incorporate sigmaij-force into virial calculation
          Real3D r21 = p1pos - p2pos;
          ww = Tensor(r21, force) / fabs(r21[2]);

          int position1 = (int)( p1pos[2]/z_dist );
          int position2 = (int)( p2pos[2]/z_dist );

          int maxpos = std::max(position1, position2);
          int minpos = std::min(position1, position2);

          // boundaries should be taken into account
          bool boundaries1 = false;
          bool boundaries2 = false;
          if(minpos < 0){
            minpos += n;
            boundaries1 =true;
          }
          if(maxpos >=n){
            maxpos -= n;
            boundaries2 =true;
          }

          if(boundaries1 || boundaries2){
            for(int i = 0; i<=maxpos; i++){
              wlocal[i] += ww;
            }
            for(int i = minpos+1; i<n; i++){
              wlocal[i] += ww;
            }
          }
          else{
            for(int i = minpos+1; i<=maxpos; i++){
              wlocal[i] += ww;
            }
          }
        }
      }

      // reduce over all CPUs
      Tensor *wsum = new Tensor[n];
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, n, (double*)&wsum, std::plus<double>());

      for(int j=0; j<n; j++){
        w[j] += wsum[j];
      }

      delete [] wsum;
      delete [] wlocal;
    }

    template < typename _Potential >
    inline real
    VerletListVSphereInteractionTemplate< _Potential >::
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
