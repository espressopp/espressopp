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
#ifndef _INTERACTION_FIXEDPAIRLISTTYPESINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDPAIRLISTTYPESINTERACTIONTEMPLATE_HPP

//#include <typeinfo>


#include "types.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "FixedPairList.hpp"
#include "FixedPairListAdress.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "types.hpp"

#include "storage/Storage.hpp"

//In FixedPairListInteractionTemplate.hpp, one potential is used for the entire FixedPairList.
//Here, the potential used for each FixedPair comes from an array of potentials, as a function of the types of the two particles in the FixedPair (i.e. same as how it works in VerletListInteractionTemplate.hpp)
//This is needed e.g. for non-bonded fixed-pair interactions, such as 1-4 LJ and Coulomb interactions in protein forcefields.

//Unlike in FixedPairListInteractionTemplate.hpp, no ltMaxBondSqr variable is calculated here. (The FixedPairListInteractionTemplate.hpp variable seems to be unused.)


namespace espressopp {
  namespace interaction {
    template < typename _Potential >
    class FixedPairListTypesInteractionTemplate: public Interaction, SystemAccess {

    protected:
      typedef _Potential Potential;

    public:
      FixedPairListTypesInteractionTemplate
          (shared_ptr < System > system,
           shared_ptr<FixedPairList> _fixedpairList)
          : SystemAccess(system), fixedpairList(_fixedpairList) {
    	  potentialArray    = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());
        ntypes = 0;
      }

      virtual ~FixedPairListTypesInteractionTemplate() {};

      void
      setFixedPairList(shared_ptr < FixedPairList > _fixedpairList) {
        fixedpairList = _fixedpairList;
      }

      shared_ptr < FixedPairList > getFixedPairList() {
        return fixedpairList;
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
      virtual int bondType() { return Pair; }

    protected:
      int ntypes;
      shared_ptr < FixedPairList > fixedpairList;
      esutil::Array2D<Potential, esutil::enlarge> potentialArray;
      esutil::Array2D<shared_ptr<Potential>, esutil::enlarge> potentialArrayPtr;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _Potential > inline void
    FixedPairListTypesInteractionTemplate < _Potential >::addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by the FixedPair List");
      const bc::BC& bc = *getSystemRef().bc;

      for (FixedPairList::PairList::Iterator it(*fixedpairList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.type();
        int type2 = p2.type();
        const Potential &potential = getPotential(type1, type2);
        // shared_ptr<Potential> potential = getPotential(type1, type2);

        Real3D force(0.0);
        //if(potential._computeForce(force, p1, p2)) {
        ////if(potential->_computeForce(force, p1, p2)) {
        //  p1.force() += force;
        //  p2.force() -= force;
        //  LOG4ESPP_TRACE(theLogger, "id1=" << p1.id() << " id2=" << p2.id() << " force=" << force);
        //}
        Real3D dist;
        bc.getMinimumImageVectorBox(dist, p1.position(), p2.position());
        if(potential._computeForce(force, p1, p2, dist)) {
          p1.force() += force;
          p2.force() -= force;
        }
      }
    }

    template < typename _Potential >
    inline real
    FixedPairListTypesInteractionTemplate < _Potential >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy of the FixedPair list pairs");

      real e = 0.0;
      real es = 0.0;
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedPairList::PairList::Iterator it(*fixedpairList);
           it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.type();
        int type2 = p2.type();
        const Potential &potential = getPotential(type1, type2);
        // shared_ptr<Potential> potential = getPotential(type1, type2);
        Real3D r21;
        bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
        //e   = potential._computeEnergy(p1, p2);
        // e   = potential->_computeEnergy(p1, p2);
        e = potential._computeEnergy(p1,p2,r21);
        es += e;
        LOG4ESPP_TRACE(theLogger, "id1=" << p1.id() << " id2=" << p2.id() << " potential energy=" << e);
        //std::cout << "id1=" << p1.id() << " id2=" << p2.id() << " potential energy=" << e << std::endl;
      }

      // reduce over all CPUs
      real esum;
      boost::mpi::all_reduce(*mpiWorld, es, esum, std::plus<real>());
      return esum;
    }

    template < typename _Potential > inline real
    FixedPairListTypesInteractionTemplate < _Potential >::
    computeEnergyDeriv() {
      std::cout << "Warning! At the moment computeEnergyDeriv() in FixedPairListTypesInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline real
    FixedPairListTypesInteractionTemplate < _Potential >::
    computeEnergyAA() {
      std::cout << "Warning! At the moment computeEnergyAA() in FixedPairListTypesInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline real
    FixedPairListTypesInteractionTemplate < _Potential >::
    computeEnergyAA(int atomtype) {
      std::cout << "Warning! At the moment computeEnergyAA(int atomtype) in FixedPairListTypesInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline real
    FixedPairListTypesInteractionTemplate < _Potential >::
    computeEnergyCG() {
      std::cout << "Warning! At the moment computeEnergyCG() in FixedPairListTypesInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential > inline real
    FixedPairListTypesInteractionTemplate < _Potential >::
    computeEnergyCG(int atomtype) {
      std::cout << "Warning! At the moment computeEnergyCG(int atomtype) in FixedPairListTypesInteractionTemplate does not work." << std::endl;
      return 0.0;
    }

    template < typename _Potential >
    inline void
    FixedPairListTypesInteractionTemplate < _Potential >::
    computeVirialX(std::vector<real> &p_xx_total, int bins) {
        std::cout << "Warning! At the moment computeVirialX in FixedPairListTypesInteractionTemplate does not work." << std::endl << "Therefore, the corresponding interactions won't be included in calculation." << std::endl;
    }

    template < typename _Potential > inline real
    FixedPairListTypesInteractionTemplate < _Potential >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute the virial for the Fixed Pair List with types");
      std::cout << "Warning! computeVirial in FixedPairListTypesInteractionTemplate has not been tested." << std::endl;

      real w = 0.0;
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedPairList::PairList::Iterator it(*fixedpairList);
           it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.type();
        int type2 = p2.type();
        const Potential &potential = getPotential(type1, type2);
        // shared_ptr<Potential> potential = getPotential(type1, type2);

        Real3D force(0.0, 0.0, 0.0);
        Real3D r21;
        bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
        if(potential._computeForce(force, p1, p2, r21)) {
        // if(potential->_computeForce(force, p1, p2)) {
          //Real3D r21 = p1.position() - p2.position();
          w = w + r21 * force;
        }
      }

      // reduce over all CPUs
      real wsum;
      boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
      return wsum;
    }

    template < typename _Potential > inline void
    FixedPairListTypesInteractionTemplate < _Potential >::computeVirialTensor(Tensor& w){

      std::cout << "Warning! At the moment computeVirialTensor() in FixedPairListTypesInteractionTemplate does not work." << std::endl;

      /*LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");

      Tensor wlocal(0.0);
      //const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedPairList::PairList::Iterator it(*fixedpairList);
           it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        Real3D r21;
        //bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
        int type1 = p1.type();
        int type2 = p2.type();
        const Potential &potential = getPotential(type1, type2);
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
      w += wsum;*/
    }

    template < typename _Potential > inline void
    FixedPairListTypesInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor& w, real z){
      std::cout << "Warning! At the moment computeVirialTensor() in FixedPairListTypesInteractionTemplate does not work." << std::endl;
      /*LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");

      Tensor wlocal(0.0);
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedPairList::PairList::Iterator it(*fixedpairList);
           it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        Real3D p1pos = p1.position();
        Real3D p2pos = p2.position();
        int type1 = p1.type();
        int type2 = p2.type();
        const Potential &potential = getPotential(type1, type2);

        if(  (p1pos[2]>=z && p2pos[2]<=z) ||
             (p1pos[2]<=z && p2pos[2]>=z) ){
          Real3D r21;
          bc.getMinimumImageVectorBox(r21, p1pos, p2pos);
          Real3D force(0.0, 0.0, 0.0);
          if(potential._computeForce(force, p1, p2)) {
          // if(potential->_computeForce(force, p1, p2)) {
            Real3D r21 = p1.position() - p2.position();
            wlocal += Tensor(r21, force);
          }
        }
      }

      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, wlocal, wsum, std::plus<Tensor>());
      w += wsum;*/
    }

    template < typename _Potential > inline void
    FixedPairListTypesInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor *w, int n){
      std::cout << "Warning! At the moment computeVirialTensor() in FixedPairListTypesInteractionTemplate does not work." << std::endl;
      /*LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");

      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      Real3D Li = bc.getBoxL();
      Tensor *wlocal = new Tensor[n];
      for(int i=0; i<n; i++) wlocal[i] = Tensor(0.0);
      for (FixedPairList::PairList::Iterator it(*fixedpairList);
           it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        Real3D p1pos = p1.position();
        Real3D p2pos = p2.position();
        int type1 = p1.type();
        int type2 = p2.type();
        const Potential &potential = getPotential(type1, type2);

        int position1 = (int)( n * p1pos[2]/Li[2]);
        int position2 = (int)( n * p1pos[2]/Li[2]);

        int maxpos = std::max(position1, position2);
        int minpos = std::min(position1, position2);

        Real3D r21;
        bc.getMinimumImageVectorBox(r21, p1pos, p2pos);
        Real3D force(0.0, 0.0, 0.0);
        Tensor ww;
        if(potential._computeForce(force, p1, p2)) {
        // if(potential->_computeForce(force, p1, p2)) {
          Real3D r21 = p1.position() - p2.position();
          ww += Tensor(r21, force);
        }

        int i = minpos + 1;
        while(i<=maxpos){
          wlocal[i] += ww;
          i++;
        }
      }

      Tensor *wsum = new Tensor[n];
      boost::mpi::all_reduce(*mpiWorld, wlocal, n, wsum, std::plus<Tensor>());

      for(int j=0; j<n; j++){
        w[j] += wsum[j];
      }

      delete [] wsum;
      delete [] wlocal;*/
    }

    template < typename _Potential >
    inline real
    FixedPairListTypesInteractionTemplate< _Potential >::getMaxCutoff() {
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
