/*
  Copyright (C) 2017,2018
      Max Planck Institute for Polymer Research

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
#ifndef _INTERACTION_FIXEDLOCALTUPLECOMLISTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDLOCALTUPLECOMLISTINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "FixedLocalTupleList.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "storage/Storage.hpp"
#include "Interaction.hpp"
#include "types.hpp"
#include "esutil/Error.hpp"

namespace espressopp {
    namespace interaction {
	template < typename _Potential >
	class FixedLocalTupleComListInteractionTemplate: public Interaction, SystemAccess {

	protected:
	    typedef _Potential Potential;

	private:
	    long unsigned int num_of_part; // total number of particles
	    long unsigned int num_of_subchain; // total number of subchains
	    boost::unordered_map<long unsigned int, Real3D> com_origin;
	    boost::unordered_map<long unsigned int, Real3D> com_temporal;
	    boost::unordered_map<long unsigned int, real> total_mass;

	    int N_Constrain; // number of particles constraining the interaction
	    int nParticles;  // number of particles in this node

	public:
	    FixedLocalTupleComListInteractionTemplate
	    (shared_ptr < System > _system,
	     shared_ptr < FixedLocalTupleList > _fixedtupleList,
	     shared_ptr < Potential > _potential)
	      : SystemAccess(_system), fixedtupleList(_fixedtupleList), potential(_potential) {
		if (! potential) {
		    LOG4ESPP_ERROR(theLogger, "NULL potential");
		}
		System& system = getSystemRef();
		esutil::Error err(system.comm);

		long unsigned int localN = system.storage -> getNRealParticles();
		boost::mpi::all_reduce(*system.comm, localN, num_of_part, std::plus<int>());

		N_Constrain = 0;

		if (!fixedtupleList->empty()) {

		    FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList);
		    std::vector<Particle*> pList = it->second;

		    N_Constrain = pList.size() + 1;

		    // Check the length of individula particle list in FixedLocalTuple
		    for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
			std::vector<Particle*> pList = it->second;

			if (pList.size() + 1 != N_Constrain) {
			    std::stringstream msg;
			    msg << "ERROR: Tuple Length is not constant\n";
			    err.setException(msg.str());
			}
		    }
		}

		std::vector<int> nList;
		boost::mpi::all_gather(*system.comm, N_Constrain, nList);
		int check_N = 0;
		for (int i = 0; i < nList.size(); i++) {
		    if (nList[i] != 0) {
			check_N = nList[i];
			break;
		    }
		}
		for (int i = 0; i < nList.size(); i++) {
		    if (nList[i] != 0) {
			if (nList[i] != check_N) {
			    std::stringstream msg;
			    msg << "ERROR: Tuple Length is not constant\n";
			    err.setException(msg.str());
			}
		    }
		}
		N_Constrain = check_N;

		num_of_subchain = num_of_part/N_Constrain;

		int* pid_list;
		longint pidK = 0;
		pid_list = new int[num_of_subchain];

		// Check the particle id in FixedLocalTuple
		for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
		    Particle *p = it->first;
		    pid_list[pidK] = p->id();

		    for(long unsigned int i = 0; i < pidK; i++) {
			if (pid_list[i]/N_Constrain == pid_list[pidK]/N_Constrain) {
			    std::stringstream msg;
			    msg << "ERROR: Particle ID is redundant\n";
			    err.setException(msg.str());
			}
		    }
		    pidK++;
		}

		delete [] pid_list;
		pid_list = NULL;

		Real3D* totCOM = new Real3D[num_of_subchain];
		Real3D* COM = new Real3D[num_of_subchain];

		computeCOM();

		for(long unsigned int i = 0; i < num_of_subchain; i++) {
		    totCOM[i] = 0.;

		    if (com_temporal.count(i)) {
			COM[i] = com_temporal[i];
		    } else {
			COM[i] = 0.;
		    }
		}

		boost::mpi::all_reduce(*system.comm, COM, num_of_subchain, totCOM, std::plus<Real3D>() );

		for (long unsigned int i = 0; i < num_of_subchain; i++) {
		    com_origin[i] = totCOM[i];
		}

		delete [] COM;
		COM = NULL;
		delete [] totCOM;
		totCOM = NULL;
	    }

	    virtual ~FixedLocalTupleComListInteractionTemplate() {};

	    // set the center of mass of subchains
	    void setCom(longint id, const Real3D& pos) {
		com_origin[id] = pos;
	    }

	    void
	    setFixedLocalTupleList(shared_ptr < FixedLocalTupleList > _fixedtupleList) {
		fixedtupleList = _fixedtupleList;
	    }

	    shared_ptr < FixedLocalTupleList > getFixedLocalTupleList() {
		return fixedtupleList;
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
	    virtual real getMaxCutoff() { return 0.0; }
	    virtual int bondType() { return Nonbonded; }

	protected:
	    void computeCOM() {

		com_temporal.clear();
		total_mass.clear();

		const bc::BC& bc = *getSystemRef().bc;  // boundary conditions

		for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
		    Particle *p = it->first;
		    Real3D pos_i = p->position();
		    real mass = p->mass();
		    Real3D tmp_com = mass*pos_i;
		    real tmp_mass = mass;
		    std::vector<Particle*> pList = it->second;

		    for (long unsigned int j = 0; j < N_Constrain - 1; j++) {
			p = pList[j];
			Real3D pos_j = p->position();
			Real3D dist;
			bc.getMinimumImageVectorBox(dist, pos_j, pos_i);
			pos_j = pos_i + dist;
			mass = p->mass();
			tmp_com += mass*pos_j;
			tmp_mass += mass;
			pos_i = pos_j;
		    }
		    com_temporal[(p->id() - 1)/N_Constrain] = tmp_com/tmp_mass;
		    total_mass[(p->id() - 1)/N_Constrain] = tmp_mass;
		}
	    }

	    int ntypes;
	    shared_ptr < FixedLocalTupleList > fixedtupleList;
	    shared_ptr < Potential > potential;
	};

	//////////////////////////////////////////////////
	// INLINE IMPLEMENTATION
	//////////////////////////////////////////////////
	template < typename _Potential > inline void
	FixedLocalTupleComListInteractionTemplate < _Potential >::addForces() {
	    LOG4ESPP_INFO(_Potential::theLogger, "adding forces of FixedLocalTupleComList");

	    const bc::BC& bc = *getSystemRef().bc;  // boundary conditions

	    computeCOM();

	    for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
		Particle *p = it->first;
		std::vector<Particle*> pList = it->second;
		Real3D diff, force;
		long unsigned int pid = p->id() - 1;
		bc.getMinimumImageVectorBox(diff,
					    com_origin[pid/N_Constrain],
					    com_temporal[pid/N_Constrain]);
		force = potential->_computeForce(diff, total_mass[pid/N_Constrain]);
		p->force() += p->mass()*force;

		for (int i = 0; i < N_Constrain - 1; i++) {
		    p = pList[i];
		    p->force() += p->mass()*force;
		}
	    }
	}

	template < typename _Potential > inline real
	FixedLocalTupleComListInteractionTemplate < _Potential >::
	computeEnergy() {

	    LOG4ESPP_INFO(theLogger, "compute energy of the FixedLocalTupleComList pairs");

	    const bc::BC& bc = *getSystemRef().bc;  // boundary conditions

	    real e = 0.0;
	    computeCOM();

	    for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
		Particle *p = it->first;
		Real3D diff;
		long unsigned int pid = p->id() - 1;
		bc.getMinimumImageVectorBox(diff,
					    com_origin[pid/N_Constrain],
					    com_temporal[pid/N_Constrain]);

		e += potential->_computeEnergy(diff);
	    }

	    real esum;
	    boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
	    return esum;
	}

	template < typename _Potential > inline real
	FixedLocalTupleComListInteractionTemplate < _Potential >::
	computeEnergyDeriv() {
	    std::cout << "Warning! At the moment computeEnergyDeriv() in FixedLocalTupleComListInteractionTemplate does not work." << std::endl;
	    return 0.0;
	}

	template < typename _Potential > inline real
	FixedLocalTupleComListInteractionTemplate < _Potential >::
	computeEnergyAA() {
	    std::cout << "Warning! At the moment computeEnergyAA() in FixedLocalTupleComListInteractionTemplate does not work." << std::endl;
	    return 0.0;
	}

	template < typename _Potential > inline real
	FixedLocalTupleComListInteractionTemplate < _Potential >::
	computeEnergyAA(int atomtype) {
	    std::cout << "Warning! At the moment computeEnergyAA(int atomtype) in FixedLocalTupleComListInteractionTemplate does not work." << std::endl;
	    return 0.0;
	}

	template < typename _Potential > inline real
	FixedLocalTupleComListInteractionTemplate < _Potential >::
	computeEnergyCG() {
	    std::cout << "Warning! At the moment computeEnergyCG() in FixedLocalTupleComListInteractionTemplate does not work." << std::endl;
	    return 0.0;
	}

	template < typename _Potential > inline real
	FixedLocalTupleComListInteractionTemplate < _Potential >::
	computeEnergyCG(int atomtype) {
	    std::cout << "Warning! At the moment computeEnergyCG(int atomtype) in FixedLocalTupleComListInteractionTemplate does not work." << std::endl;
	    return 0.0;
	}

	template < typename _Potential >
	inline void
	FixedLocalTupleComListInteractionTemplate < _Potential >::
	computeVirialX(std::vector<real> &p_xx_total, int bins) {
	    std::cout << "Warning! At the moment computeVirialX in FixedLocalTupleComListInteractionTemplate does not work." << std::endl << "Therefore, the corresponding interactions won't be included in calculation." << std::endl;
	}


	template < typename _Potential > inline real
	FixedLocalTupleComListInteractionTemplate < _Potential >::
	computeVirial() {
	    LOG4ESPP_INFO(theLogger, "compute the virial of the FixedLocalTupleCom");

	    real w = 0.0;
	    const bc::BC& bc = *getSystemRef().bc;  // boundary conditions

	    computeCOM();

	    for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
		Particle *p = it->first;
		std::vector<Particle*> pList = it->second;
		Real3D diff, dist, force;
		long unsigned int pid = p->id() - 1;
		bc.getMinimumImageVectorBox(diff,
					    com_origin[pid/N_Constrain],
					    com_temporal[pid/N_Constrain]);
		bc.getMinimumImageVectorBox(dist,
					    p->position(),
					    com_origin[pid/N_Constrain]);
		force = potential->_computeForce(diff, total_mass[pid/N_Constrain]);

		// TODO: formulas are not correct yet?
		w += dist * p->mass() * force;

		for (int i = 0; i < N_Constrain - 1; i++) {
		    p = pList[i];
		    bc.getMinimumImageVectorBox(dist,
						p->position(),
						com_origin[pid/N_Constrain]);
		    w += dist * p->mass() * force;
		}
	    }

	    real wsum;
	    boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
	    return wsum;
	}

	template < typename _Potential > inline void
	FixedLocalTupleComListInteractionTemplate < _Potential >::computeVirialTensor(Tensor& w){
	    LOG4ESPP_INFO(theLogger, "compute the virial tensor of the FixedLocalTupleCom");

	    Tensor wlocal(0.0);
	    const bc::BC& bc = *getSystemRef().bc;  // boundary conditions

	    computeCOM();

	    for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
		Particle *p = it->first;
		std::vector<Particle*> pList = it->second;
		Real3D diff, dist, force;
		long unsigned int pid = p->id() - 1;
		bc.getMinimumImageVectorBox(diff,
					    com_origin[pid/N_Constrain],
					    com_temporal[pid/N_Constrain]);
		bc.getMinimumImageVectorBox(dist,
					    p->position(),
					    com_origin[pid/N_Constrain]);
		force = potential->_computeForce(diff, total_mass[pid/N_Constrain]);

		// TODO: formulas are not correct yet?
		wlocal += Tensor(dist, p->mass()*force);

		for (int i = 0; i < N_Constrain - 1; i++) {
		    p = pList[i];
		    bc.getMinimumImageVectorBox(dist,
						p->position(),
						com_origin[pid/N_Constrain]);
		    wlocal += Tensor(dist, p->mass()*force);
		}
	    }

	    // reduce over all CPUs
	    Tensor wsum(0.0);
	    boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
	    w += wsum;
	}

	template < typename _Potential > inline void
	FixedLocalTupleComListInteractionTemplate < _Potential >::
	computeVirialTensor(Tensor& w, real z){
	    std::cout<<"Warning! Calculating virial layerwise is not supported for "
		"in FixedLocalTupleComListInteractionTemplate."<<std::endl;
	}

	template < typename _Potential > inline void
	FixedLocalTupleComListInteractionTemplate < _Potential >::
	computeVirialTensor(Tensor *w, int n){
	    std::cout<<"Warning! Calculating virial layerwise is not supported for "
		"in FixedLocalTupleComListInteractionTemplate."<<std::endl;
	}
    }
}
#endif
