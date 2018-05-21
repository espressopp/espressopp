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
#ifndef _INTERACTION_VERLETLISTPIADRESSNODRIFTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_VERLETLISTPIADRESSNODRIFTINTERACTIONTEMPLATE_HPP

#include "System.hpp"
#include "bc/BC.hpp"
#include "types.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "VerletListAdress.hpp"
#include "FixedTupleListAdress.hpp"
#include "esutil/Array2D.hpp"
#include "SystemAccess.hpp"

namespace espressopp {
  namespace interaction {
    template < typename _Potential >
    class VerletListPIadressNoDriftInteractionTemplate: public Interaction {

      protected:
        typedef _Potential Potential;

      public:
        VerletListPIadressNoDriftInteractionTemplate
        (shared_ptr<VerletListAdress> _verletList, shared_ptr<FixedTupleListAdress> _fixedtupleList, int _ntrotter, bool _speedup)
          : verletList(_verletList), fixedtupleList(_fixedtupleList), ntrotter(_ntrotter), speedup(_speedup) {

          potentialArray = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());
          ntypes = 0;
        }

        void
        setVerletList(shared_ptr < VerletListAdress > _verletList) {
          verletList = _verletList;
        }

        shared_ptr<VerletListAdress> getVerletList() {
          return verletList;
        }

        void
        setFixedTupleList(shared_ptr<FixedTupleListAdress> _fixedtupleList) {
          fixedtupleList = _fixedtupleList;
        }

        shared_ptr < FixedTupleListAdress > getFixedTupleList() {
          return fixedtupleList;
        }

        void
        setNTrotter(int _ntrotter) {
          ntrotter = _ntrotter;
        }

        int getNTrotter() {
          return ntrotter;
        }

        void
        setSpeedup(bool _speedup) {
          speedup = _speedup;
        }

        bool getSpeedup() {
          return speedup;
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

        Potential &getPotential(int type1, int type2) {
          return potentialArray.at(type1, type2);
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
        virtual int bondType() {
          return Nonbonded;
        }

      protected:
        int ntypes;
        shared_ptr<VerletListAdress> verletList;
        shared_ptr<FixedTupleListAdress> fixedtupleList;
        esutil::Array2D<Potential, esutil::enlarge> potentialArray;

        int ntrotter;
        bool speedup; // if true approximate rings in classical region by single particles
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _Potential > inline void
    VerletListPIadressNoDriftInteractionTemplate < _Potential >::
    addForces() {
      // Pairs not inside the QM/Hybrid Zone (i.e. CL region)
      for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;

        // Calculate forces in CL region
        if (speedup == true) {
          const Potential &potential = getPotential(p1.type(), p2.type());
          Real3D forcecl(0.0, 0.0, 0.0);
          if(potential._computeForce(forcecl, p1, p2)) {
            p1.force() += forcecl;
            p2.force() -= forcecl;
          }
        }
        else {
          FixedTupleListAdress::iterator it3;
          FixedTupleListAdress::iterator it4;
          it3 = fixedtupleList->find(&p1);
          it4 = fixedtupleList->find(&p2);

          if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

            // Get the PI bead lists (i.e. the AdResS particles)
            std::vector<Particle*> atList1;
            std::vector<Particle*> atList2;
            atList1 = it3->second;
            atList2 = it4->second;

            // Iterate the two iterators in a parallel fashion
            std::vector<Particle*>::iterator itv2 = atList2.begin();
            for (std::vector<Particle*>::iterator itv = atList1.begin();
                 itv != atList1.end(); ++itv) {

              if (itv2 == atList2.end()) {
                std::stringstream ss;
                ss << "Tuplelists seem to have different lengths or are not started properly. Corresponding CG particles " << p1.id() << " and " << p2.id();
                throw std::runtime_error(ss.str());
              }

              // Get the individual PI beads
              Particle &p3 = **itv;
              Particle &p4 = **itv2;

              if (p3.pib() != p4.pib()) {
                std::stringstream ss;
                ss << "Path Integral bead numbers do not match in VerletListPIadressInteractionTemplate for particles " << p3.id() << " and " << p4.id();
                throw std::runtime_error(ss.str());
              }

              // Calculate forces
              const Potential &potential = getPotential(p3.type(), p4.type());
              Real3D forcecl(0.0, 0.0, 0.0);
              if(potential._computeForce(forcecl, p3, p4)) {
                forcecl *= 1.0/ntrotter;
                p3.force() += forcecl;
                p4.force() -= forcecl;
              }

              //Iterate the second iterator
              ++itv2;

            }

          }
          else {
            std::stringstream ss;
            ss << "One of the VP particles not found in tuples. Particle 1 - ID: " << p1.id() << " Ghost: " << p1.ghost() << " Position: " << p1.position() << " Particle 2 - ID: " << p2.id() << " Ghost: " << p2.ghost() << " Position: " << p2.position();
            throw std::runtime_error(ss.str());
          }

        }

      }

      // Pairs inside the QM/Hybrid Zone
      for (PairList::Iterator it(verletList->getAdrPairs()); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;

        FixedTupleListAdress::iterator it3;
        FixedTupleListAdress::iterator it4;
        it3 = fixedtupleList->find(&p1);
        it4 = fixedtupleList->find(&p2);

        if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

          // Get the PI bead lists (i.e. the AdResS particles)
          std::vector<Particle*> atList1;
          std::vector<Particle*> atList2;
          atList1 = it3->second;
          atList2 = it4->second;

          // Iterate the two iterators in a parallel fashion
          std::vector<Particle*>::iterator itv2 = atList2.begin();
          for (std::vector<Particle*>::iterator itv = atList1.begin();
               itv != atList1.end(); ++itv) {

            if (itv2 == atList2.end()) {
              std::stringstream ss;
              ss << "Tuplelists seem to have different lengths or are not started properly. Corresponding CG particles " << p1.id() << " and " << p2.id();
              throw std::runtime_error(ss.str());
            }

            // Get the individual PI beads
            Particle &p3 = **itv;
            Particle &p4 = **itv2;

            if (p3.pib() != p4.pib()) {
              std::stringstream ss;
              ss << "Path Integral bead numbers do not match in VerletListPIadressInteractionTemplate for particles " << p3.id() << " and " << p4.id();
              throw std::runtime_error(ss.str());
            }

            // Calculate forces
            const Potential &potential = getPotential(p3.type(), p4.type());
            Real3D forceqm(0.0, 0.0, 0.0);
            if(potential._computeForce(forceqm, p3, p4)) {
              forceqm *= 1.0/ntrotter;
              p3.force() += forceqm;
              p4.force() -= forceqm;
            }

            //Iterate the second iterator
            ++itv2;

          }

        }
        else {
          std::stringstream ss;
          ss << "One of the VP particles not found in tuples. Particle 1 - ID: " << p1.id() << " Ghost: " << p1.ghost() << " Position: " << p1.position() << " Particle 2 - ID: " << p2.id() << " Ghost: " << p2.ghost() << " Position: " << p2.position();
          throw std::runtime_error(ss.str());
        }

      }

    }

    // Energy calculation does currently only work if integrator.run( ) (also with 0) and decompose have been executed before. This is due to the initialization of the tuples.
    template < typename _Potential >
    inline real
    VerletListPIadressNoDriftInteractionTemplate < _Potential >::
    computeEnergy() {
      // Pairs not inside the QM/Hybrid Zone (i.e. CL region)
      real e = 0.0;
      for (PairList::Iterator it(verletList->getPairs());
           it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        int type1 = p1.type();
        int type2 = p2.type();

        if(speedup == true) {
          const Potential &potential = getPotential(type1, type2);
          e += potential._computeEnergy(p1, p2);
        }
        else {
          FixedTupleListAdress::iterator it3;
          FixedTupleListAdress::iterator it4;
          it3 = fixedtupleList->find(&p1);
          it4 = fixedtupleList->find(&p2);

          if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

            // Get the PI bead lists (i.e. the AdResS particles)
            std::vector<Particle*> atList1;
            std::vector<Particle*> atList2;
            atList1 = it3->second;
            atList2 = it4->second;

            // Iterate the two iterators in a parallel fashion
            std::vector<Particle*>::iterator itv2 = atList2.begin();
            for (std::vector<Particle*>::iterator itv = atList1.begin();
                 itv != atList1.end(); ++itv) {

              if (itv2 == atList2.end()) {
                std::stringstream ss;
                ss << "Tuplelists seem to have different lengths or are not started properly. Corresponding CG particles " << p1.id() << " and " << p2.id();
                throw std::runtime_error(ss.str());
              }

              // Get the individual PI beads
              Particle &p3 = **itv;
              Particle &p4 = **itv2;

              if (p3.pib()!= p4.pib()) {
                std::stringstream ss;
                ss << "Path Integral bead numbers do not match in VerletListPIadressInteractionTemplate for particles " << p3.id() << " and " << p4.id();
                throw std::runtime_error(ss.str());
              }

              // Calculate energy
              const Potential &potential = getPotential(p3.type(), p4.type());
              e += (1.0/ntrotter)*potential._computeEnergy(p3, p4);

              //Iterate the second iterator
              ++itv2;

            }

          }
          else {
            std::stringstream ss;
            ss << "One of the VP particles not found in tuples. Particle 1 - ID: " << p1.id() << " Ghost: " << p1.ghost() << " Position: " << p1.position() << " Particle 2 - ID: " << p2.id() << " Ghost: " << p2.ghost() << " Position: " << p2.position();
            throw std::runtime_error(ss.str());
          }
        }
      }

      // Pairs inside the QM/Hybrid Zone
      for (PairList::Iterator it(verletList->getAdrPairs());
           it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;

        FixedTupleListAdress::iterator it3;
        FixedTupleListAdress::iterator it4;
        it3 = fixedtupleList->find(&p1);
        it4 = fixedtupleList->find(&p2);

        if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

          // Get the PI bead lists (i.e. the AdResS particles)
          std::vector<Particle*> atList1;
          std::vector<Particle*> atList2;
          atList1 = it3->second;
          atList2 = it4->second;

          // Iterate the two iterators in a parallel fashion
          std::vector<Particle*>::iterator itv2 = atList2.begin();
          for (std::vector<Particle*>::iterator itv = atList1.begin();
               itv != atList1.end(); ++itv) {

            if (itv2 == atList2.end()) {
              std::stringstream ss;
              ss << "Tuplelists seem to have different lengths or are not started properly. Corresponding CG particles " << p1.id() << " and " << p2.id();
              throw std::runtime_error(ss.str());
            }

            // Get the individual PI beads
            Particle &p3 = **itv;
            Particle &p4 = **itv2;

            if (p3.pib()!= p4.pib()) {
              std::stringstream ss;
              ss << "Path Integral bead numbers do not match in VerletListPIadressInteractionTemplate for particles " << p3.id() << " and " << p4.id();
              throw std::runtime_error(ss.str());
            }

            // Calculate energy
            const Potential &potential = getPotential(p3.type(), p4.type());
            e += (1.0/ntrotter)*potential._computeEnergy(p3, p4);

            //Iterate the second iterator
            ++itv2;

          }

        }
        else {
          std::stringstream ss;
          ss << "One of the VP particles not found in tuples. Particle 1 - ID: " << p1.id() << " Ghost: " << p1.ghost() << " Position: " << p1.position() << " Particle 2 - ID: " << p2.id() << " Ghost: " << p2.ghost() << " Position: " << p2.position();
          throw std::runtime_error(ss.str());
        }

      }
      real esum;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, e, esum, std::plus<real>());
      return esum;
    }


    template < typename _Potential > inline real
    VerletListPIadressNoDriftInteractionTemplate < _Potential >::
    computeEnergyDeriv() {
      std::cout << "Warning! At the moment computeEnergyDeriv() in VerletListPIadressNoDriftInteractionTemplate does not work." << std::endl;
      return 0.0;
    }


    template < typename _Potential >
    inline real
    VerletListPIadressNoDriftInteractionTemplate < _Potential >::
    computeEnergyAA() {
      real e = 0.0;

      for (PairList::Iterator it(verletList->getAdrPairs());
           it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;

        FixedTupleListAdress::iterator it3;
        FixedTupleListAdress::iterator it4;
        it3 = fixedtupleList->find(&p1);
        it4 = fixedtupleList->find(&p2);

        if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

          // Get the PI bead lists (i.e. the AdResS particles)
          std::vector<Particle*> atList1;
          std::vector<Particle*> atList2;
          atList1 = it3->second;
          atList2 = it4->second;

          // Iterate the two iterators in a parallel fashion
          std::vector<Particle*>::iterator itv2 = atList2.begin();
          for (std::vector<Particle*>::iterator itv = atList1.begin();
               itv != atList1.end(); ++itv) {

            if (itv2 == atList2.end()) {
              std::stringstream ss;
              ss << "Tuplelists seem to have different lengths or are not started properly. Corresponding CG particles " << p1.id() << " and " << p2.id();
              throw std::runtime_error(ss.str());
            }

            // Get the individual PI beads
            Particle &p3 = **itv;
            Particle &p4 = **itv2;

            if (p3.pib() != p4.pib()) {
              std::stringstream ss;
              ss << "Path Integral bead numbers do not match in VerletListPIadressInteractionTemplate for particles " << p3.id() << " and " << p4.id();
              throw std::runtime_error(ss.str());
            }

            // Calculate energy
            const Potential &potential = getPotential(p3.type(), p4.type());
            e += (1.0/ntrotter)*potential._computeEnergy(p3, p4);

            //Iterate the second iterator
            ++itv2;

          }

        }
        else {
          std::stringstream ss;
          ss << "One of the VP particles not found in tuples. Particle 1 - ID: " << p1.id() << " Ghost: " << p1.ghost() << " Position: " << p1.position() << " Particle 2 - ID: " << p2.id() << " Ghost: " << p2.ghost() << " Position: " << p2.position();
          throw std::runtime_error(ss.str());
        }

      }

      for (PairList::Iterator it(verletList->getPairs());
           it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;

        FixedTupleListAdress::iterator it3;
        FixedTupleListAdress::iterator it4;
        it3 = fixedtupleList->find(&p1);
        it4 = fixedtupleList->find(&p2);

        if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

          // Get the PI bead lists (i.e. the AdResS particles)
          std::vector<Particle*> atList1;
          std::vector<Particle*> atList2;
          atList1 = it3->second;
          atList2 = it4->second;

          // Iterate the two iterators in a parallel fashion
          std::vector<Particle*>::iterator itv2 = atList2.begin();
          for (std::vector<Particle*>::iterator itv = atList1.begin();
               itv != atList1.end(); ++itv) {

            if (itv2 == atList2.end()) {
              std::stringstream ss;
              ss << "Tuplelists seem to have different lengths or are not started properly. Corresponding CG particles " << p1.id() << " and " << p2.id();
              throw std::runtime_error(ss.str());
            }

            // Get the individual PI beads
            Particle &p3 = **itv;
            Particle &p4 = **itv2;

            if (p3.pib() != p4.pib()) {
              std::stringstream ss;
              ss << "Path Integral bead numbers do not match in VerletListPIadressInteractionTemplate for particles " << p3.id() << " and " << p4.id();
              throw std::runtime_error(ss.str());
            }

            // Calculate energy
            const Potential &potential = getPotential(p3.type(), p4.type());
            e += (1.0/ntrotter)*potential._computeEnergy(p3, p4);

            //Iterate the second iterator
            ++itv2;

          }

        }
        else {
          std::stringstream ss;
          ss << "One of the VP particles not found in tuples. Particle 1 - ID: " << p1.id() << " Ghost: " << p1.ghost() << " Position: " << p1.position() << " Particle 2 - ID: " << p2.id() << " Ghost: " << p2.ghost() << " Position: " << p2.position();
          throw std::runtime_error(ss.str());
        }

      }

      real esum;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, e, esum, std::plus<real>());
      return esum;
    }


    template < typename _Potential >
    inline real
    VerletListPIadressNoDriftInteractionTemplate < _Potential >::
    computeEnergyAA(int atomtype) {
      real e = 0.0;

      for (PairList::Iterator it(verletList->getAdrPairs());
           it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;

        if((p1.type() == atomtype) || (p2.type() == atomtype)) {

          FixedTupleListAdress::iterator it3;
          FixedTupleListAdress::iterator it4;
          it3 = fixedtupleList->find(&p1);
          it4 = fixedtupleList->find(&p2);

          if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

            // Get the PI bead lists (i.e. the AdResS particles)
            std::vector<Particle*> atList1;
            std::vector<Particle*> atList2;
            atList1 = it3->second;
            atList2 = it4->second;

            // Iterate the two iterators in a parallel fashion
            std::vector<Particle*>::iterator itv2 = atList2.begin();
            for (std::vector<Particle*>::iterator itv = atList1.begin();
                 itv != atList1.end(); ++itv) {

              if (itv2 == atList2.end()) {
                std::stringstream ss;
                ss << "Tuplelists seem to have different lengths or are not started properly. Corresponding CG particles " << p1.id() << " and " << p2.id();
                throw std::runtime_error(ss.str());
              }

              // Get the individual PI beads
              Particle &p3 = **itv;
              Particle &p4 = **itv2;

              if (p3.pib() != p4.pib()) {
                std::stringstream ss;
                ss << "Path Integral bead numbers do not match in VerletListPIadressInteractionTemplate for particles " << p3.id() << " and " << p4.id();
                throw std::runtime_error(ss.str());
              }


              // Calculate energy
              const Potential &potential = getPotential(p3.type(), p4.type());

              if((p1.type() == atomtype) &&  (p2.type() == atomtype)) {
                e += (1.0/ntrotter)*potential._computeEnergy(p3, p4);
              }
              else {
                e += 0.5*(1.0/ntrotter)*potential._computeEnergy(p3, p4);
              }
              //Iterate the second iterator
              ++itv2;

            }

          }
          else {
            std::stringstream ss;
            ss << "One of the VP particles not found in tuples. Particle 1 - ID: " << p1.id() << " Ghost: " << p1.ghost() << " Position: " << p1.position() << " Particle 2 - ID: " << p2.id() << " Ghost: " << p2.ghost() << " Position: " << p2.position();
            throw std::runtime_error(ss.str());
          }

        }

      }

      for (PairList::Iterator it(verletList->getPairs());
           it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;

        if((p1.type() == atomtype) || (p2.type() == atomtype)) {

          FixedTupleListAdress::iterator it3;
          FixedTupleListAdress::iterator it4;
          it3 = fixedtupleList->find(&p1);
          it4 = fixedtupleList->find(&p2);

          if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

            // Get the PI bead lists (i.e. the AdResS particles)
            std::vector<Particle*> atList1;
            std::vector<Particle*> atList2;
            atList1 = it3->second;
            atList2 = it4->second;

            // Iterate the two iterators in a parallel fashion
            std::vector<Particle*>::iterator itv2 = atList2.begin();
            for (std::vector<Particle*>::iterator itv = atList1.begin();
                 itv != atList1.end(); ++itv) {

              if (itv2 == atList2.end()) {
                std::stringstream ss;
                ss << "Tuplelists seem to have different lengths or are not started properly. Corresponding CG particles " << p1.id() << " and " << p2.id();
                throw std::runtime_error(ss.str());
              }

              // Get the individual PI beads
              Particle &p3 = **itv;
              Particle &p4 = **itv2;

              if (p3.pib() != p4.pib()) {
                std::stringstream ss;
                ss << "Path Integral bead numbers do not match in VerletListPIadressInteractionTemplate for particles " << p3.id() << " and " << p4.id();
                throw std::runtime_error(ss.str());
              }

              // Calculate energy
              const Potential &potential = getPotential(p3.type(), p4.type());

              if((p1.type() == atomtype) &&  (p2.type() == atomtype)) {
                e += (1.0/ntrotter)*potential._computeEnergy(p3, p4);
              }
              else {
                e += 0.5*(1.0/ntrotter)*potential._computeEnergy(p3, p4);
              }

              //Iterate the second iterator
              ++itv2;

            }

          }
          else {
            std::stringstream ss;
            ss << "One of the VP particles not found in tuples. Particle 1 - ID: " << p1.id() << " Ghost: " << p1.ghost() << " Position: " << p1.position() << " Particle 2 - ID: " << p2.id() << " Ghost: " << p2.ghost() << " Position: " << p2.position();
            throw std::runtime_error(ss.str());
          }

        }

      }

      real esum;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, e, esum, std::plus<real>());
      return esum;
    }


    template < typename _Potential >
    inline real
    VerletListPIadressNoDriftInteractionTemplate < _Potential >::
    computeEnergyCG() {
      real e = 0.0;

      for (PairList::Iterator it(verletList->getAdrPairs());
           it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;

        FixedTupleListAdress::iterator it3;
        FixedTupleListAdress::iterator it4;
        it3 = fixedtupleList->find(&p1);
        it4 = fixedtupleList->find(&p2);

        if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

          // Get the PI bead lists (i.e. the AdResS particles)
          std::vector<Particle*> atList1;
          std::vector<Particle*> atList2;
          atList1 = it3->second;
          atList2 = it4->second;

          // Iterate the two iterators in a parallel fashion
          std::vector<Particle*>::iterator itv2 = atList2.begin();
          for (std::vector<Particle*>::iterator itv = atList1.begin();
               itv != atList1.end(); ++itv) {

            if (itv2 == atList2.end()) {
              std::stringstream ss;
              ss << "Tuplelists seem to have different lengths or are not started properly. Corresponding CG particles " << p1.id() << " and " << p2.id();
              throw std::runtime_error(ss.str());
            }

            // Get the individual PI beads
            Particle &p3 = **itv;
            Particle &p4 = **itv2;

            if (p3.pib() != p4.pib()) {
              std::stringstream ss;
              ss << "Path Integral bead numbers do not match in VerletListPIadressInteractionTemplate for particles " << p3.id() << " and " << p4.id();
              throw std::runtime_error(ss.str());
            }

            // Calculate energy
            const Potential &potential = getPotential(p3.type(), p4.type());
            e += (1.0/ntrotter)*potential._computeEnergy(p3, p4);

            //Iterate the second iterator
            ++itv2;

          }

        }
        else {
          std::stringstream ss;
          ss << "One of the VP particles not found in tuples. Particle 1 - ID: " << p1.id() << " Ghost: " << p1.ghost() << " Position: " << p1.position() << " Particle 2 - ID: " << p2.id() << " Ghost: " << p2.ghost() << " Position: " << p2.position();
          throw std::runtime_error(ss.str());
        }

      }

      for (PairList::Iterator it(verletList->getPairs());
           it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;

        FixedTupleListAdress::iterator it3;
        FixedTupleListAdress::iterator it4;
        it3 = fixedtupleList->find(&p1);
        it4 = fixedtupleList->find(&p2);

        if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

          // Get the PI bead lists (i.e. the AdResS particles)
          std::vector<Particle*> atList1;
          std::vector<Particle*> atList2;
          atList1 = it3->second;
          atList2 = it4->second;

          // Iterate the two iterators in a parallel fashion
          std::vector<Particle*>::iterator itv2 = atList2.begin();
          for (std::vector<Particle*>::iterator itv = atList1.begin();
               itv != atList1.end(); ++itv) {

            if (itv2 == atList2.end()) {
              std::stringstream ss;
              ss << "Tuplelists seem to have different lengths or are not started properly. Corresponding CG particles " << p1.id() << " and " << p2.id();
              throw std::runtime_error(ss.str());
            }

            // Get the individual PI beads
            Particle &p3 = **itv;
            Particle &p4 = **itv2;

            if (p3.pib() != p4.pib()) {
              std::stringstream ss;
              ss << "Path Integral bead numbers do not match in VerletListPIadressInteractionTemplate for particles " << p3.id() << " and " << p4.id();
              throw std::runtime_error(ss.str());
            }

            // Calculate energy
            const Potential &potential = getPotential(p3.type(), p4.type());
            e += (1.0/ntrotter)*potential._computeEnergy(p3, p4);

            //Iterate the second iterator
            ++itv2;

          }

        }
        else {
          std::stringstream ss;
          ss << "One of the VP particles not found in tuples. Particle 1 - ID: " << p1.id() << " Ghost: " << p1.ghost() << " Position: " << p1.position() << " Particle 2 - ID: " << p2.id() << " Ghost: " << p2.ghost() << " Position: " << p2.position();
          throw std::runtime_error(ss.str());
        }

      }

      real esum;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, e, esum, std::plus<real>());
      return esum;
    }


    template < typename _Potential >
    inline real
    VerletListPIadressNoDriftInteractionTemplate < _Potential >::
    computeEnergyCG(int atomtype) {
      real e = 0.0;

      for (PairList::Iterator it(verletList->getAdrPairs());
           it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;

        if((p1.type() == atomtype) || (p2.type() == atomtype)) {

          FixedTupleListAdress::iterator it3;
          FixedTupleListAdress::iterator it4;
          it3 = fixedtupleList->find(&p1);
          it4 = fixedtupleList->find(&p2);

          if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

            // Get the PI bead lists (i.e. the AdResS particles)
            std::vector<Particle*> atList1;
            std::vector<Particle*> atList2;
            atList1 = it3->second;
            atList2 = it4->second;

            // Iterate the two iterators in a parallel fashion
            std::vector<Particle*>::iterator itv2 = atList2.begin();
            for (std::vector<Particle*>::iterator itv = atList1.begin();
                 itv != atList1.end(); ++itv) {

              if (itv2 == atList2.end()) {
                std::stringstream ss;
                ss << "Tuplelists seem to have different lengths or are not started properly. Corresponding CG particles " << p1.id() << " and " << p2.id();
                throw std::runtime_error(ss.str());
              }

              // Get the individual PI beads
              Particle &p3 = **itv;
              Particle &p4 = **itv2;

              if (p3.pib() != p4.pib()) {
                std::stringstream ss;
                ss << "Path Integral bead numbers do not match in VerletListPIadressInteractionTemplate for particles " << p3.id() << " and " << p4.id();
                throw std::runtime_error(ss.str());
              }

              // Calculate energy
              const Potential &potential = getPotential(p3.type(), p4.type());

              if((p1.type() == atomtype) &&  (p2.type() == atomtype)) {
                e += (1.0/ntrotter)*potential._computeEnergy(p3, p4);
              }
              else {
                e += 0.5*(1.0/ntrotter)*potential._computeEnergy(p3, p4);
              }

              //Iterate the second iterator
              ++itv2;

            }

          }
          else {
            std::stringstream ss;
            ss << "One of the VP particles not found in tuples. Particle 1 - ID: " << p1.id() << " Ghost: " << p1.ghost() << " Position: " << p1.position() << " Particle 2 - ID: " << p2.id() << " Ghost: " << p2.ghost() << " Position: " << p2.position();
            throw std::runtime_error(ss.str());
          }

        }

      }

      for (PairList::Iterator it(verletList->getPairs());
           it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;

        if((p1.type() == atomtype) || (p2.type() == atomtype)) {

          FixedTupleListAdress::iterator it3;
          FixedTupleListAdress::iterator it4;
          it3 = fixedtupleList->find(&p1);
          it4 = fixedtupleList->find(&p2);

          if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

            // Get the PI bead lists (i.e. the AdResS particles)
            std::vector<Particle*> atList1;
            std::vector<Particle*> atList2;
            atList1 = it3->second;
            atList2 = it4->second;

            // Iterate the two iterators in a parallel fashion
            std::vector<Particle*>::iterator itv2 = atList2.begin();
            for (std::vector<Particle*>::iterator itv = atList1.begin();
                 itv != atList1.end(); ++itv) {

              if (itv2 == atList2.end()) {
                std::stringstream ss;
                ss << "Tuplelists seem to have different lengths or are not started properly. Corresponding CG particles " << p1.id() << " and " << p2.id();
                throw std::runtime_error(ss.str());
              }

              // Get the individual PI beads
              Particle &p3 = **itv;
              Particle &p4 = **itv2;

              if (p3.pib() != p4.pib()) {
                std::stringstream ss;
                ss << "Path Integral bead numbers do not match in VerletListPIadressInteractionTemplate for particles " << p3.id() << " and " << p4.id();
                throw std::runtime_error(ss.str());
              }

              // Calculate energy
              const Potential &potential = getPotential(p3.type(), p4.type());

              if((p1.type() == atomtype) &&  (p2.type() == atomtype)) {
                e += (1.0/ntrotter)*potential._computeEnergy(p3, p4);
              }
              else {
                e += 0.5*(1.0/ntrotter)*potential._computeEnergy(p3, p4);
              }

              //Iterate the second iterator
              ++itv2;

            }

          }
          else {
            std::stringstream ss;
            ss << "One of the VP particles not found in tuples. Particle 1 - ID: " << p1.id() << " Ghost: " << p1.ghost() << " Position: " << p1.position() << " Particle 2 - ID: " << p2.id() << " Ghost: " << p2.ghost() << " Position: " << p2.position();
            throw std::runtime_error(ss.str());
          }

        }

      }

      real esum;
      boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, e, esum, std::plus<real>());
      return esum;
    }


    template < typename _Potential > inline void
    VerletListPIadressNoDriftInteractionTemplate < _Potential >::
    computeVirialX(std::vector<real> &p_xx_total, int bins) {
      std::cout << "Warning! At the moment computeVirialX in VerletListPIadressNoDrift does not work"<<std::endl;
      return;
    }


    template < typename _Potential > inline real
    VerletListPIadressNoDriftInteractionTemplate < _Potential >::
    computeVirial() {
      real w = 0.0;

      // Pairs not inside the QM/Hybrid Zone (i.e. CL region)
      for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it) {
        // Get particles from pairlist
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;

        // Calculate forces in CL region
        if (speedup == true) {
          const Potential &potential = getPotential(p1.type(), p2.type());
          Real3D forcecl(0.0, 0.0, 0.0);
          if(potential._computeForce(forcecl, p1, p2)) {
            Real3D dist(0.0, 0.0, 0.0);
            verletList->getSystem()->bc->getMinimumImageVector(dist, p1.position(), p2.position());
            w += dist * forcecl;
          }
        }
        else {
          FixedTupleListAdress::iterator it3;
          FixedTupleListAdress::iterator it4;
          it3 = fixedtupleList->find(&p1);
          it4 = fixedtupleList->find(&p2);

          if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

            // Get the PI bead lists (i.e. the AdResS particles)
            std::vector<Particle*> atList1;
            std::vector<Particle*> atList2;
            atList1 = it3->second;
            atList2 = it4->second;

            // Iterate the two iterators in a parallel fashion
            std::vector<Particle*>::iterator itv2 = atList2.begin();
            for (std::vector<Particle*>::iterator itv = atList1.begin();
                 itv != atList1.end(); ++itv) {

              if (itv2 == atList2.end()) {
                std::stringstream ss;
                ss << "Tuplelists seem to have different lengths or are not started properly. Corresponding CG particles " << p1.id() << " and " << p2.id();
                throw std::runtime_error(ss.str());
              }

              // Get the individual PI beads
              Particle &p3 = **itv;
              Particle &p4 = **itv2;

              if (p3.pib() != p4.pib()) {
                std::stringstream ss;
                ss << "Path Integral bead numbers do not match in VerletListPIadressInteractionTemplate for particles " << p3.id() << " and " << p4.id();
                throw std::runtime_error(ss.str());
              }

              // Calculate CL forces
              const Potential &potential = getPotential(p3.type(), p4.type());
              Real3D forcecl(0.0, 0.0, 0.0);
              if(potential._computeForce(forcecl, p3, p4)) {
                forcecl *= 1.0/ntrotter;
                Real3D dist(0.0, 0.0, 0.0);
                verletList->getSystem()->bc->getMinimumImageVector(dist, p3.position(), p4.position());
                w += dist * forcecl;
              }

              //Iterate the second iterator
              ++itv2;

            }

          }
          else {
            std::stringstream ss;
            ss << "One of the VP particles not found in tuples. Particle 1 - ID: " << p1.id() << " Ghost: " << p1.ghost() << " Position: " << p1.position() << " Particle 2 - ID: " << p2.id() << " Ghost: " << p2.ghost() << " Position: " << p2.position();
            throw std::runtime_error(ss.str());
          }

        }

      }

      // Pairs inside the QM/Hybrid Zone
      for (PairList::Iterator it(verletList->getAdrPairs()); it.isValid(); ++it) {
        real w1, w2;
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;

        FixedTupleListAdress::iterator it3;
        FixedTupleListAdress::iterator it4;
        it3 = fixedtupleList->find(&p1);
        it4 = fixedtupleList->find(&p2);

        if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

          // Get the PI bead lists (i.e. the AdResS particles)
          std::vector<Particle*> atList1;
          std::vector<Particle*> atList2;
          atList1 = it3->second;
          atList2 = it4->second;

          // Iterate the two iterators in a parallel fashion
          std::vector<Particle*>::iterator itv2 = atList2.begin();
          for (std::vector<Particle*>::iterator itv = atList1.begin();
               itv != atList1.end(); ++itv) {

            if (itv2 == atList2.end()) {
              std::stringstream ss;
              ss << "Tuplelists seem to have different lengths or are not started properly. Corresponding CG particles " << p1.id() << " and " << p2.id();
              throw std::runtime_error(ss.str());
            }

            // Get the individual PI beads
            Particle &p3 = **itv;
            Particle &p4 = **itv2;

            if (p3.pib() != p4.pib()) {
              std::stringstream ss;
              ss << "Path Integral bead numbers do not match in VerletListPIadressInteractionTemplate for particles " << p3.id() << " and " << p4.id();
              throw std::runtime_error(ss.str());
            }

            const Potential &potential = getPotential(p3.type(), p4.type());
            Real3D forceqm(0.0, 0.0, 0.0);
            if(potential._computeForce(forceqm, p3, p4)) {
              forceqm *= 1.0/ntrotter;
              Real3D dist(0.0, 0.0, 0.0);
              verletList->getSystem()->bc->getMinimumImageVector(dist, p3.position(), p4.position());
              w += dist * forceqm;
            }

            //Iterate the second iterator
            ++itv2;

          }

        }
        else {
          std::stringstream ss;
          ss << "One of the VP particles not found in tuples. Particle 1 - ID: " << p1.id() << " Ghost: " << p1.ghost() << " Position: " << p1.position() << " Particle 2 - ID: " << p2.id() << " Ghost: " << p2.ghost() << " Position: " << p2.position();
          throw std::runtime_error(ss.str());
        }

      }

      real wsum;
      wsum = 0.0;
      boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
      return wsum;
    }

    template < typename _Potential > inline void
    VerletListPIadressNoDriftInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor& w) {
      std::cout << "Warning! At the moment IK computeVirialTensor in VerletListPIadressNoDrift does not work"<<std::endl;
      return;
    }

    template < typename _Potential > inline void
    VerletListPIadressNoDriftInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor& w, real z) {
      std::cout << "Warning! At the moment IK computeVirialTensor in VerletListPIadressNoDrift does not work"<<std::endl;
      return;
    }

    template < typename _Potential > inline void
    VerletListPIadressNoDriftInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor *w, int n) {
      std::cout << "Warning! At the moment IK computeVirialTensor in VerletListPIadressNoDrift does not work"<<std::endl;
      return;
    }

    template < typename _Potential >
    inline real
    VerletListPIadressNoDriftInteractionTemplate< _Potential >::
    getMaxCutoff() {
      real cutoff = 0.0;
      for (int i = 0; i < ntypes; i++) {
        for (int j = 0; j < ntypes; j++) {
          cutoff = std::max(cutoff, getPotential(i, j).getCutoff());
        }
      }
      return cutoff;
    }
  }
}
#endif
