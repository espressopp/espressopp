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
#ifndef _INTERACTION_FIXEDLOCALTUPLERGLISTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDLOCALTUPLERGLISTINTERACTIONTEMPLATE_HPP

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
	class FixedLocalTupleRgListInteractionTemplate: public Interaction, SystemAccess {

	protected:
	    typedef _Potential Potential;

	private:
	    long unsigned int num_of_part; // total number of particles
	    long unsigned int num_of_subchain; // total number of subchains
	    boost::unordered_map<long unsigned int, real> rg_origin;

	    int N_Constrain; // number of particles constraining the interaction

	public:
	    FixedLocalTupleRgListInteractionTemplate
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

		boost::unordered_map<long unsigned int, Real3D> center = computeCenter();

		real* totRG = new real[num_of_subchain];
		real* RG = new real[num_of_subchain];

		boost::unordered_map<long unsigned int, real> rg = computeRG(center);

		for(long unsigned int i = 0; i < num_of_subchain; i++) {
		    totRG[i] = 0.;

		    if (rg.count(i)) {
			RG[i] = rg[i];
		    } else {
			RG[i] = 0.;
		    }
		}

		boost::mpi::all_reduce(*system.comm, RG, num_of_subchain, totRG, std::plus<real>() );

		for(long unsigned int i = 0; i < num_of_subchain; i++) {
		    rg_origin[i] = totRG[i];
		}

		delete [] RG;
		RG = NULL;
		delete [] totRG;
		totRG = NULL;
	    }

	    virtual ~FixedLocalTupleRgListInteractionTemplate() {};

	    // set the center of mass of subchains
	    void setRG(longint id, const real rg) {
		rg_origin[id] = rg;
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
	    // calcuclate the center of a subchain \mathbf{C} in the local cell list.
	    // \mathbf{C} = \frac{1}{N_Constrain} \sum^{N_Constrain}_{i=1} \mathbf{r}_i
	    // \mathbf{C} is not the center of mass of a subchain but the center of a subchain
	    boost::unordered_map<long unsigned int, Real3D> computeCenter() {

		boost::unordered_map<long unsigned int, Real3D> center;

		const bc::BC& bc = *getSystemRef().bc;  // boundary conditions

		for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
		    Particle *p = it->first;
		    Real3D pos_i = p->position();
		    Real3D tmp_center = pos_i;
		    std::vector<Particle*> pList = it->second;

		    for (long unsigned int j = 1; j < N_Constrain; j++) {
			p = pList[j - 1];
			Real3D pos_j = p->position();
			Real3D dist;
			bc.getMinimumImageVectorBox(dist, pos_j, pos_i);
			pos_j = pos_i + dist;
			tmp_center += pos_j;
			pos_i = pos_j;
		    }
		    center[(p->id() - 1)/N_Constrain] = tmp_center/N_Constrain;
		}

		return center;
	    }

	    // calcuclate the gyration radius of subchain in the local cell list.
	    boost::unordered_map<long unsigned int, real>
	    computeRG(boost::unordered_map<long unsigned int, Real3D> center) {
		boost::unordered_map<long unsigned int, real> rg;

		const bc::BC& bc = *getSystemRef().bc;  // boundary conditions

		for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
		    Particle *p = it->first;
		    Real3D pos_i = p->position();
		    std::vector<Particle*> pList = it->second;
		    long unsigned int id = (p->id() - 1)/N_Constrain;

		    Real3D dist;
		    bc.getMinimumImageVectorBox(dist, pos_i, center[id]);
		    real tmp_rg = dist.sqr();
		    for (long unsigned int j = 1; j < N_Constrain; j++) {
			p = pList[j - 1];
			Real3D pos_j = p->position();
			bc.getMinimumImageVectorBox(dist, pos_j, pos_i);
			pos_j = pos_i + dist;
			bc.getMinimumImageVectorBox(dist, pos_j, center[id]);
			tmp_rg += dist.sqr();
			pos_i = pos_j;
		    }
		    rg[id] = tmp_rg/N_Constrain;
		}

		return rg;
	    }

	    int ntypes;
	    shared_ptr < FixedLocalTupleList > fixedtupleList;
	    shared_ptr < Potential > potential;
	};

	//////////////////////////////////////////////////
	// INLINE IMPLEMENTATION
	//////////////////////////////////////////////////
	template < typename _Potential > inline void
	FixedLocalTupleRgListInteractionTemplate < _Potential >::addForces() {
	    LOG4ESPP_INFO(_Potential::theLogger, "adding forces of FixedLocalTupleRgList");

	    const bc::BC& bc = *getSystemRef().bc;  // boundary conditions

	    boost::unordered_map<long unsigned int, Real3D> center;
	    center = computeCenter();

	    boost::unordered_map<long unsigned int, real> rg_sq;
	    rg_sq = computeRG(center);

	    for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
		Particle *p = it->first;
		std::vector<Particle*> pList = it->second;
		Real3D diff;
		long unsigned int id = (p->id() - 1)/N_Constrain;

		bc.getMinimumImageVectorBox(diff,
					    p->position(),
					    center[id]);

		real diff_rg = rg_origin[id] - rg_sq[id];

		p->force() += p->mass()*potential->_computeForce(diff, diff_rg, N_Constrain);

		for (long unsigned int i = 0; i < N_Constrain - 1; i++) {
		    Particle* pt = pList[i];
		    bc.getMinimumImageVectorBox(diff,
						pt->position(),
						center[id]);
		    pt->force() +=
			p->mass()*potential->_computeForce(diff, diff_rg, N_Constrain);
		}

	    }
	}

	template < typename _Potential > inline real
	FixedLocalTupleRgListInteractionTemplate < _Potential >::
	computeEnergy() {

	    LOG4ESPP_INFO(theLogger, "compute energy of the FixedLocalTupleRgList pairs");

	    real e = 0.0;
	    boost::unordered_map<long unsigned int, Real3D> center;
	    center = computeCenter();

	    boost::unordered_map<long unsigned int, real> rg;
	    rg = computeRG(center);

	    for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
		Particle *p = it->first;

		long unsigned int id = (p->id() - 1)/N_Constrain;
		real diff = rg_origin[id] - rg[id];

		e += N_Constrain*potential->_computeEnergy(diff);
	    }

	    real esum;
	    boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
	    return esum;
	}

	template < typename _Potential > inline real
	FixedLocalTupleRgListInteractionTemplate < _Potential >::
	computeEnergyDeriv() {
	    std::cout << "Warning! At the moment computeEnergyDeriv() in FixedLocalTupleRgListInteractionTemplate does not work." << std::endl;
	    return 0.0;
	}

	template < typename _Potential > inline real
	FixedLocalTupleRgListInteractionTemplate < _Potential >::
	computeEnergyAA() {
	    std::cout << "Warning! At the moment computeEnergyAA() in FixedLocalTupleRgListInteractionTemplate does not work." << std::endl;
	    return 0.0;
	}

	template < typename _Potential > inline real
	FixedLocalTupleRgListInteractionTemplate < _Potential >::
	computeEnergyAA(int atomtype) {
	    std::cout << "Warning! At the moment computeEnergyAA(int atomtype) in FixedLocalTupleRgListInteractionTemplate does not work." << std::endl;
	    return 0.0;
	}

	template < typename _Potential > inline real
	FixedLocalTupleRgListInteractionTemplate < _Potential >::
	computeEnergyCG() {
	    std::cout << "Warning! At the moment computeEnergyCG() in FixedLocalTupleRgListInteractionTemplate does not work." << std::endl;
	    return 0.0;
	}

	template < typename _Potential > inline real
	FixedLocalTupleRgListInteractionTemplate < _Potential >::
	computeEnergyCG(int atomtype) {
	    std::cout << "Warning! At the moment computeEnergyCG(int atomtype) in FixedLocalTupleRgListInteractionTemplate does not work." << std::endl;
	    return 0.0;
	}

	template < typename _Potential >
	inline void
	FixedLocalTupleRgListInteractionTemplate < _Potential >::
	computeVirialX(std::vector<real> &p_xx_total, int bins) {
	    std::cout << "Warning! At the moment computeVirialX in FixedLocalTupleRgListInteractionTemplate does not work." << std::endl << "Therefore, the corresponding interactions won't be included in calculation." << std::endl;
	}


	template < typename _Potential > inline real
	FixedLocalTupleRgListInteractionTemplate < _Potential >::
	computeVirial() {
	    LOG4ESPP_INFO(theLogger, "compute the virial of the FixedLocalTupleRg");

	    real w = 0.0;
	    const bc::BC& bc = *getSystemRef().bc;  // boundary conditions

	    boost::unordered_map<long unsigned int, Real3D> center;
	    center = computeCenter();

	    boost::unordered_map<long unsigned int, real> rg_sq;
	    rg_sq = computeRG(center);

	    for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
		Particle *p = it->first;
		std::vector<Particle*> pList = it->second;
		Real3D diff;
		long unsigned int id = (p->id() - 1)/N_Constrain;

		bc.getMinimumImageVectorBox(diff,
					    p->position(),
					    center[id]);

		real diff_rg = rg_origin[id] - rg_sq[id];

		w += diff*p->mass()*potential->_computeForce(diff, diff_rg, N_Constrain);

		for (long unsigned int i = 0; i < N_Constrain - 1; i++) {
		    Particle* pt = pList[i];
		    bc.getMinimumImageVectorBox(diff,
						pt->position(),
						center[id]);
		    w += diff*p->mass()*potential->_computeForce(diff, diff_rg, N_Constrain);
		}

	    }

	    real wsum;
	    boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
	    return wsum;
	}

	template < typename _Potential > inline void
	FixedLocalTupleRgListInteractionTemplate < _Potential >::computeVirialTensor(Tensor& w){
	    LOG4ESPP_INFO(theLogger, "compute the virial tensor of the FixedLocalTupleRg");

	    Tensor wlocal(0.0);
	    const bc::BC& bc = *getSystemRef().bc;  // boundary conditions

	    boost::unordered_map<long unsigned int, Real3D> center;
	    center = computeCenter();

	    boost::unordered_map<long unsigned int, real> rg_sq;
	    rg_sq = computeRG(center);

	    for (FixedLocalTupleList::TupleList::Iterator it(*fixedtupleList); it.isValid(); ++it) {
		Particle *p = it->first;
		std::vector<Particle*> pList = it->second;
		Real3D diff;
		long unsigned int id = (p->id() - 1)/N_Constrain;

		bc.getMinimumImageVectorBox(diff,
					    p->position(),
					    center[id]);

		real diff_rg = rg_origin[id] - rg_sq[id];

		wlocal += Tensor(diff, p->mass()*potential->_computeForce(diff, diff_rg, N_Constrain));

		for (long unsigned int i = 0; i < N_Constrain - 1; i++) {
		    Particle* pt = pList[i];
		    bc.getMinimumImageVectorBox(diff,
						pt->position(),
						center[id]);
		    wlocal += Tensor(diff, p->mass()*potential->_computeForce(diff, diff_rg, N_Constrain));
		}

	    }

	    // reduce over all CPUs
	    Tensor wsum(0.0);
	    boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
	    w += wsum;
	}

	template < typename _Potential > inline void
	FixedLocalTupleRgListInteractionTemplate < _Potential >::
	computeVirialTensor(Tensor& w, real z){
	    std::cout<<"Warning! Calculating virial layerwise is not supported for "
		"in FixedLocalTupleRgListInteractionTemplate."<<std::endl;
	}

	template < typename _Potential > inline void
	FixedLocalTupleRgListInteractionTemplate < _Potential >::
	computeVirialTensor(Tensor *w, int n){
	    std::cout<<"Warning! Calculating virial layerwise is not supported for "
		"in FixedLocalTupleRgListInteractionTemplate."<<std::endl;
	}
    }
}
#endif
