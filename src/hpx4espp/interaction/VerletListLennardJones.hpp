/*
  Copyright (C) 2019-2022
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
#ifndef _HPX4ESPP_INTERACTION_VERLETLISTLENNARDJONES_HPP
#define _HPX4ESPP_INTERACTION_VERLETLISTLENNARDJONES_HPP

//#include <typeinfo>

#include <hpx/config.hpp>

#include "types.hpp"
#include "interaction/Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "hpx4espp/VerletList.hpp"
#include "hpx4espp/utils/multithreading.hpp"
#include "hpx4espp/utils/assert.hpp"
#include "hpx4espp/utils/algorithms/for_loop.hpp"

#include "storage/Storage.hpp"
#include "LennardJones.hpp"

#include "vec/interaction/VerletListLennardJones.hpp"

namespace espressopp
{
namespace hpx4espp
{
namespace interaction
{
using espressopp::hpx4espp::VerletList;
using espressopp::interaction::Interaction;
using espressopp::interaction::Nonbonded;

class VerletListLennardJones : public Interaction
{
protected:
    typedef LennardJones _Potential;
    typedef _Potential Potential;

public:
    typedef vec::interaction::VerletListLennardJones::LJCoefficients LJCoefficients;

    VerletListLennardJones(shared_ptr<VerletList> _verletList) : verletList(_verletList)
    {
        potentialArray =
            espressopp::esutil::Array2D<Potential, espressopp::esutil::enlarge>(0, 0, Potential());
        ntypes = 0;
        np_types = 0;
        p_types = 0;
    }

    virtual ~VerletListLennardJones(){};

    void setVerletList(shared_ptr<VerletList> _verletList) { verletList = _verletList; }

    shared_ptr<VerletList> getVerletList() { return verletList; }

    void setPotential(int type1, int type2, const Potential& potential)
    {
        // typeX+1 because i<ntypes
        ntypes = std::max(ntypes, std::max(type1 + 1, type2 + 1));
        potentialArray.at(type1, type2) = potential;
        LOG4ESPP_INFO(_Potential::theLogger,
                      "added potential for type1=" << type1 << " type2=" << type2);
        if (type1 != type2)
        {  // add potential in the other direction
            potentialArray.at(type2, type1) = potential;
            LOG4ESPP_INFO(_Potential::theLogger, "automatically added the same potential for type1="
                                                     << type2 << " type2=" << type1);
        }
        rebuildPotential();
    }

    // this is used in the innermost force-loop
    Potential& getPotential(int type1, int type2)
    {
        if (type1 >= potentialArray.size_n() || type2 >= potentialArray.size_m())
            needRebuildPotential = true;
        return potentialArray.at(type1, type2);
    }

    // this is mainly used to access the potential from Python (e.g. to change parameters of the
    // potential)
    shared_ptr<Potential> getPotentialPtr(int type1, int type2)
    {
        return make_shared<Potential>(potentialArray.at(type1, type2));
    }

    void rebuildPotential()
    {
        np_types = potentialArray.size_n();
        p_types = potentialArray.size_m();
        ffs = vec::AlignedVector<LJCoefficients>(np_types * p_types);
        cutoffSqr = vec::AlignedVector<real>(np_types * p_types);
        vec::AlignedVector<LJCoefficients>::iterator it1 = ffs.begin();
        vec::AlignedVector<real>::iterator it3 = cutoffSqr.begin();
        for (auto& p : potentialArray)
        {
            *(it1++) = LJCoefficients(p.getff1(), p.getff2());
            *(it3++) = p.getCutoffSqr();
        }
        needRebuildPotential = false;
    }
    virtual void addForces();
    virtual real computeEnergy();
    virtual real computeEnergyDeriv();
    virtual real computeEnergyAA();
    virtual real computeEnergyCG();
    virtual real computeEnergyAA(int atomtype);
    virtual real computeEnergyCG(int atomtype);
    virtual void computeVirialX(std::vector<real>& p_xx_total, int bins);
    virtual real computeVirial();
    virtual void computeVirialTensor(Tensor& w);
    virtual void computeVirialTensor(Tensor& w, real z);
    virtual void computeVirialTensor(Tensor* w, int n);
    virtual real getMaxCutoff();
    virtual int bondType() { return Nonbonded; }

protected:
    template <bool ONETYPE>
    void addForces_impl();

    template <bool ONETYPE, bool N3L>
    void addForces_impl_vs(vec::ParticleArray& particleArray,
                           vec::ParticleArray& particleArrayNbr,
                           VerletList::NeighborList const& neighborList);

    int ntypes;
    shared_ptr<VerletList> verletList;
    espressopp::esutil::Array2D<Potential, espressopp::esutil::enlarge> potentialArray;
    // not needed espressopp::esutil::Array2D<shared_ptr<Potential>, espressopp::esutil::enlarge>
    // potentialArrayPtr;

    size_t np_types, p_types;
    vec::AlignedVector<LJCoefficients> ffs;
    vec::AlignedVector<real> cutoffSqr;
    bool needRebuildPotential = true;
};

//////////////////////////////////////////////////
// INLINE IMPLEMENTATION
//////////////////////////////////////////////////
inline void VerletListLennardJones::addForces()
{
    LOG4ESPP_DEBUG(_Potential::theLogger, "loop over verlet list pairs and add forces");

    // lookup table for LJ variables
    // ideal for low number of types vs number of particle pairs
    // trigger rebuild on setPotential and on size modification from getPotential
    int vlmaxtype = verletList->getMaxType();
    Potential max_pot = getPotential(vlmaxtype, vlmaxtype);
    if (needRebuildPotential) rebuildPotential();
    if (np_types == 1 && p_types == 1)
        addForces_impl<true>();
    else
        addForces_impl<false>();
}

template <bool ONETYPE>
inline void VerletListLennardJones::addForces_impl()
{
    auto& vss = verletList->getVirtualStorage();
    auto& nls = verletList->getNeighborLists();

    /// for loop over virtual subdomains
    const size_t nvs = vss.size();
    utils::parallelForLoop(size_t(0), nvs, [this, &vss, &nls](size_t inode) {
        auto& vs = vss[inode];
        const auto& nl = nls[inode];
        addForces_impl_vs<ONETYPE, 1>(vs.particles, vs.particles, nl.realNbrs);
        addForces_impl_vs<ONETYPE, 1>(vs.particles, vs.particles, nl.externalNbrs);

        for (const auto& nnode_cnl : nl.internalNbrs)
        {
            const auto nnode = nnode_cnl.first;
            const auto& cnl = nnode_cnl.second;
            addForces_impl_vs<ONETYPE, 0>(vs.particles, vss[nnode].particles, cnl);
        }
    });
}

template <bool ONETYPE, bool N3L>
inline void VerletListLennardJones::addForces_impl_vs(vec::ParticleArray& pa,
                                                      vec::ParticleArray& paNbr,
                                                      VerletList::NeighborList const& nl)
{
    vec::interaction::VerletListLennardJones::addForces_impl<ONETYPE, N3L>(pa, paNbr, nl, ffs,
                                                                           cutoffSqr, np_types);
}

inline real VerletListLennardJones::computeEnergy()
{
#if 0
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
#endif
    LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeEnergy() is not yet implemented.");
    return 0.0;
}

inline real VerletListLennardJones::computeEnergyDeriv()
{
    LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeEnergyDeriv() is not yet implemented.");
    return 0.0;
}

inline real VerletListLennardJones::computeEnergyAA()
{
    LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeEnergyAA() is not yet implemented.");
    return 0.0;
}

inline real VerletListLennardJones::computeEnergyAA(int atomtype)
{
    LOG4ESPP_WARN(_Potential::theLogger,
                  "Warning! computeEnergyAA(int atomtype) is not yet implemented.");
    return 0.0;
}

inline real VerletListLennardJones::computeEnergyCG()
{
    LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeEnergyCG() is not yet implemented.");
    return 0.0;
}

inline real VerletListLennardJones::computeEnergyCG(int atomtype)
{
    LOG4ESPP_WARN(_Potential::theLogger,
                  "Warning! computeEnergyCG(int atomtype) is not yet implemented.");
    return 0.0;
}

inline void VerletListLennardJones::computeVirialX(std::vector<real>& p_xx_total, int bins)
{
    LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeVirialX() is not yet implemented.");
}

inline real VerletListLennardJones::computeVirial()
{
#if 0
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
#endif
    LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeVirial() is not yet implemented.");
    return 0.0;
}

inline void VerletListLennardJones::computeVirialTensor(Tensor& w)
{
#if 0
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
#endif
    LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeVirialTensor() is not yet implemented.");
}

// local pressure tensor for layer, plane is defined by z coordinate
inline void VerletListLennardJones::computeVirialTensor(Tensor& w, real z)
{
#if 0
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
#endif
    LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeVirialTensor() is not yet implemented.");
}

// it will calculate the pressure in 'n' layers along Z axis
// the first layer has coordinate 0.0 the last - (Lz - Lz/n)
inline void VerletListLennardJones::computeVirialTensor(Tensor* w, int n)
{
#if 0
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
#endif
    LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeVirialTensor() is not yet implemented.");
}

inline real VerletListLennardJones::getMaxCutoff()
{
    real cutoff = 0.0;
    for (int i = 0; i < ntypes; i++)
    {
        for (int j = 0; j < ntypes; j++)
        {
            cutoff = std::max(cutoff, getPotential(i, j).getCutoff());
            // cutoff = std::max(cutoff, getPotential(i, j)->getCutoff());
        }
    }
    return cutoff;
}
}  // namespace interaction
}  // namespace hpx4espp
}  // namespace espressopp
#endif
