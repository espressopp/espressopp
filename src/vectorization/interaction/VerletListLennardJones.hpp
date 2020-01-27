/*
  Copyright (C) 2019-2020
      Max Planck Institute for Polymer Research & JGU Mainz
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
#ifndef _VECTORIZATION_INTERACTION_VERLETLISTLENNARDJONES_HPP
#define _VECTORIZATION_INTERACTION_VERLETLISTLENNARDJONES_HPP

//#include <typeinfo>


#include "types.hpp"
#include "interaction/Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "vectorization/Vectorization.hpp"
#include "vectorization/VerletList.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"

#include "storage/Storage.hpp"
#include "LennardJones.hpp"

namespace espressopp { namespace vectorization {
  namespace interaction {

    using espressopp::vectorization::VerletList;
    using espressopp::interaction::Interaction;
    using espressopp::interaction::Nonbonded;

    class VerletListLennardJones: public Interaction {

    protected:
      typedef LennardJones _Potential;
      typedef _Potential Potential;

      struct LJCoefficients
      {
        LJCoefficients(real const& ff1, real const& ff2): ff1(ff1), ff2(ff2){}
        LJCoefficients(){}
        real ff1, ff2;
      };

    public:
      VerletListLennardJones(shared_ptr<VerletList> _verletList)
        : verletList(_verletList)
      {
        potentialArray    = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());
        ntypes = 0;
        np_types = 0;
        p_types = 0;
      }

      virtual ~VerletListLennardJones() {};

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
        LOG4ESPP_INFO(_Potential::theLogger, "added potential for type1=" << type1 << " type2=" << type2);
        if (type1 != type2) { // add potential in the other direction
           potentialArray.at(type2, type1) = potential;
           LOG4ESPP_INFO(_Potential::theLogger, "automatically added the same potential for type1=" << type2 << " type2=" << type1);
        }
        rebuildPotential();
      }

      // this is used in the innermost force-loop
      Potential &getPotential(int type1, int type2) {
        if(type1>=potentialArray.size_n() || type2>=potentialArray.size_m()) needRebuildPotential = true;
        return potentialArray.at(type1, type2);
      }

      // this is mainly used to access the potential from Python (e.g. to change parameters of the potential)
      shared_ptr<Potential> getPotentialPtr(int type1, int type2) {
        return  make_shared<Potential>(potentialArray.at(type1, type2));
      }

      void rebuildPotential()
      {
        np_types = potentialArray.size_n();
        p_types = potentialArray.size_m();
        ffs = AlignedVector<LJCoefficients>(np_types*p_types);
        cutoffSqr = AlignedVector<real>(np_types*p_types);
        AlignedVector<LJCoefficients>::iterator it1 = ffs.begin();
        AlignedVector<real>::iterator it3 = cutoffSqr.begin();
        for(auto& p: potentialArray){
          *(it1++) = LJCoefficients(p.getff1(),p.getff2());
          *(it3++) = p.getCutoffSqr();
        }
        needRebuildPotential = false;
      }
      template <bool ONETYPE, bool VEC_MODE_AOS> void addForces_impl();
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
      // not needed esutil::Array2D<shared_ptr<Potential>, esutil::enlarge> potentialArrayPtr;

      size_t np_types, p_types;
      AlignedVector<LJCoefficients> ffs;
      AlignedVector<real> cutoffSqr;
      bool needRebuildPotential = true;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    inline void
    VerletListLennardJones::
    addForces() {
      LOG4ESPP_DEBUG(_Potential::theLogger, "loop over verlet list pairs and add forces");

      // lookup table for LJ variables
      // ideal for low number of types vs number of particle pairs
      // trigger rebuild on setPotential and on size modification from getPotential
      int vlmaxtype = verletList->getMaxType();
      Potential max_pot = getPotential(vlmaxtype,vlmaxtype);
      if(needRebuildPotential) rebuildPotential();
      bool VEC_MODE_AOS = verletList->getParticleArray().mode_aos();
      if(np_types==1 && p_types==1)
        if(VEC_MODE_AOS)
          addForces_impl<true,true>();
        else
          addForces_impl<true,false>();
      else
        if(VEC_MODE_AOS)
          addForces_impl<false,true>();
        else
          addForces_impl<false,false>();
    }

    template <bool ONETYPE, bool VEC_MODE_AOS>
    inline void
    VerletListLennardJones::
    addForces_impl()
    {
      {
        real ff1_, ff2_, cutoffSqr_;
        if(ONETYPE){
          ff1_ = ffs[0].ff1;
          ff2_ = ffs[0].ff2;
          cutoffSqr_ = cutoffSqr[0];
        }

        auto& particleArray = verletList->getParticleArray();
        auto& neighborList = verletList->getNeighborList();

        const Real3DInt *pa_pos  = &(particleArray.position[0]);
        Real4D *pa_force         = &(particleArray.force[0]);

        const ulongint* __restrict pa_type = &(particleArray.type[0]);
        const real* __restrict pa_p_x = &(particleArray.p_x[0]);
        const real* __restrict pa_p_y = &(particleArray.p_y[0]);
        const real* __restrict pa_p_z = &(particleArray.p_z[0]);
        real* __restrict pa_f_x       = &(particleArray.f_x[0]);
        real* __restrict pa_f_y       = &(particleArray.f_y[0]);
        real* __restrict pa_f_z       = &(particleArray.f_z[0]);

        {
          const auto* __restrict plist  = &(neighborList.plist[0]);
          const auto* __restrict prange = &(neighborList.prange[0]);
          const auto* __restrict nplist = &(neighborList.nplist[0]);
          const int ip_max = neighborList.plist.size();

          int in_min=0;
          for(int ip=0; ip<ip_max; ip++)
          {
            int p = plist[ip];
            int p_lookup;
            real p_x, p_y, p_z;

            if(VEC_MODE_AOS)
            {
              p_x = pa_pos[p].x;
              p_y = pa_pos[p].y;
              p_z = pa_pos[p].z;
              if(!ONETYPE){
                p_lookup  = pa_pos[p].t*np_types;
              }
            }
            else
            {
              if(!ONETYPE){
                p_lookup  = pa_type[p]*np_types;
              }
              p_x = pa_p_x[p];
              p_y = pa_p_y[p];
              p_z = pa_p_z[p];
            }

            real f_x              = 0.0;
            real f_y              = 0.0;
            real f_z              = 0.0;

            const int in_max=prange[ip];

            #pragma vector always
            #pragma vector aligned
            #pragma ivdep
            for(int in=in_min; in<in_max; in++)
            {
              auto np_ii = nplist[in];
              {
                int np_lookup;
                real dist_x, dist_y, dist_z;
                if(VEC_MODE_AOS)
                {
                  dist_x   = p_x - pa_pos[np_ii].x;
                  dist_y   = p_y - pa_pos[np_ii].y;
                  dist_z   = p_z - pa_pos[np_ii].z;
                  if(!ONETYPE) np_lookup = pa_pos[np_ii].t+p_lookup;
                }
                else
                {
                  dist_x   = p_x - pa_p_x[np_ii];
                  dist_y   = p_y - pa_p_y[np_ii];
                  dist_z   = p_z - pa_p_z[np_ii];
                  if(!ONETYPE) np_lookup = pa_type[np_ii]+p_lookup;
                }

                real distSqr  = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;
                if(!ONETYPE){
                  cutoffSqr_ = cutoffSqr[np_lookup];
                }

                #if defined(ESPP_VECTOR_MASK)
                if(distSqr <= cutoffSqr_)
                #endif
                {
                  real frac2 = 1.0 / distSqr;
                  real frac6 = frac2 * frac2 * frac2;
                  real ffactor;

                  if(ONETYPE)
                    ffactor = ff1_ * frac6 - ff2_;
                  else
                    ffactor = ffs[np_lookup].ff1 * frac6 - ffs[np_lookup].ff2;

                  #if !defined(ESPP_VECTOR_MASK)
                  if(distSqr > cutoffSqr_) ffactor = 0.0;
                  #endif

                  ffactor = frac6 * ffactor * frac2;

                  f_x += dist_x * ffactor;
                  f_y += dist_y * ffactor;
                  f_z += dist_z * ffactor;

                  if(VEC_MODE_AOS)
                  {
                    auto& np_force = pa_force[np_ii];
                    np_force.x -= dist_x * ffactor;
                    np_force.y -= dist_y * ffactor;
                    np_force.z -= dist_z * ffactor;
                  }
                  else
                  {
                    pa_f_x[np_ii] -= dist_x * ffactor;
                    pa_f_y[np_ii] -= dist_y * ffactor;
                    pa_f_z[np_ii] -= dist_z * ffactor;
                  }

                }
              }
            }
            if(VEC_MODE_AOS)
            {
              auto& p_force = pa_force[p];
              p_force.x += f_x;
              p_force.y += f_y;
              p_force.z += f_z;
            }
            else
            {
              pa_f_x[p] += f_x;
              pa_f_y[p] += f_y;
              pa_f_z[p] += f_z;
            }

            in_min = in_max;
          }

        }
      }
    }

    inline real
    VerletListLennardJones::
    computeEnergy() {
      LOG4ESPP_DEBUG(_Potential::theLogger, "loop over verlet list pairs and sum up potential energies");

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
        LOG4ESPP_TRACE(_Potential::theLogger, "id1=" << p1.id() << " id2=" << p2.id() << " potential energy=" << e);
      }

      // reduce over all CPUs
      real esum;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, es, esum, std::plus<real>());
      return esum;
    }

    inline real
    VerletListLennardJones::
    computeEnergyDeriv() {
      LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeEnergyDeriv() is not yet implemented.");
      return 0.0;
    }

    inline real
    VerletListLennardJones::
    computeEnergyAA() {
      LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeEnergyAA() is not yet implemented.");
      return 0.0;
    }

    inline real
    VerletListLennardJones::
    computeEnergyAA(int atomtype) {
      LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeEnergyAA(int atomtype) is not yet implemented.");
      return 0.0;
    }

    inline real
    VerletListLennardJones::
    computeEnergyCG() {
      LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeEnergyCG() is not yet implemented.");
      return 0.0;
    }

    inline real
    VerletListLennardJones::
    computeEnergyCG(int atomtype) {
      LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeEnergyCG(int atomtype) is not yet implemented.");
      return 0.0;
    }

    inline void
    VerletListLennardJones::
    computeVirialX(std::vector<real> &p_xx_total, int bins) {
      LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeVirialX() is not yet implemented.");
    }

    inline real
    VerletListLennardJones::
    computeVirial() {
      LOG4ESPP_DEBUG(_Potential::theLogger, "loop over verlet list pairs and sum up virial");

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

    inline void
    VerletListLennardJones::
    computeVirialTensor(Tensor& w) {
      LOG4ESPP_DEBUG(_Potential::theLogger, "loop over verlet list pairs and sum up virial tensor");

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
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
      w += wsum;
    }

    // local pressure tensor for layer, plane is defined by z coordinate
    inline void
    VerletListLennardJones::
    computeVirialTensor(Tensor& w, real z) {
      LOG4ESPP_DEBUG(_Potential::theLogger, "loop over verlet list pairs and sum up virial tensor over one z-layer");

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
          if(potential._computeForce(force, p1, p2)) {
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
    inline void
    VerletListLennardJones::
    computeVirialTensor(Tensor *w, int n) {
      LOG4ESPP_DEBUG(_Potential::theLogger, "loop over verlet list pairs and sum up virial tensor in bins along z-direction");

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
        if(potential._computeForce(force, p1, p2)) {
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

    inline real
    VerletListLennardJones::
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
}}
#endif
