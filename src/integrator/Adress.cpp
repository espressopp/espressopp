/*
  Copyright (C) 2012-2018
      Max Planck Institute for Polymer Research
  Copyright (C) 2008-2011
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

#include "python.hpp"
#include "mpi.hpp"
#include "Adress.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "FixedTupleListAdress.hpp"
#include "iterator/CellListAllPairsIterator.hpp"
#include "iterator/CellListIterator.hpp"
#include <iomanip>

namespace espressopp {

  namespace integrator {

    using namespace espressopp::iterator;

    Adress::Adress(shared_ptr<System> _system, shared_ptr<VerletListAdress> _verletList, shared_ptr<FixedTupleListAdress> _fixedtupleList, bool _KTI, int _regionupdates, int _multistep)
        : Extension(_system), verletList(_verletList), fixedtupleList(_fixedtupleList), KTI(_KTI), regionupdates(_regionupdates), multistep(_multistep){
        LOG4ESPP_INFO(theLogger, "construct Adress");
        type = Extension::Adress;

        // AdResS stuff
        dhy = verletList->getHy();
        pidhy2 = M_PI/(dhy * 2.0);
        dex = verletList->getEx();
        dex2 = dex * dex;
        dexdhy = dex + verletList->getHy();
        dexdhy2 = dexdhy * dexdhy;
        updatecount = 0;
        communicateAdrPositions();

    }


    Adress::~Adress() {
      LOG4ESPP_INFO(theLogger, "~Adress");
      disconnect();
    }

    void Adress::disconnect(){
        _SetPosVel.disconnect();
        _initForces.disconnect();
        _integrate1.disconnect();
        _integrate2.disconnect();
        // _inIntP.disconnect();
        //_aftCalcF.disconnect();
        _recalc2.disconnect();
        _befIntV.disconnect();
        _integrateSlow.disconnect();
        _aftCalcSlow.disconnect();
    }

    void Adress::connect() {

        // connection to after runInit()
        _SetPosVel = integrator->runInit.connect(
                boost::bind(&Adress::SetPosVel, this), boost::signals2::at_front);

        // connection to after initForces()
        _initForces = integrator->aftInitF.connect(
                boost::bind(&Adress::initForces, this), boost::signals2::at_front);

        // connection to inside of integrate1()
        _integrate1 = integrator->inIntP.connect(
                boost::bind(&Adress::integrate1, this, _1), boost::signals2::at_front);

        // // connection to inside of integrate1()
        // _inIntP = integrator->inIntP.connect(
        //         boost::bind(&Adress::communicateAdrPositions, this), boost::signals2::at_front);

        // connection to after integrate2()
        _integrate2 = integrator->aftIntV.connect(
                boost::bind(&Adress::integrate2, this), boost::signals2::at_front);

        // connection to after integrate2()
        _integrateSlow = integrator->aftIntSlow.connect(
                boost::bind(&Adress::integrateSlow, this), boost::signals2::at_front);

        // Note: Both this extension as well as Langevin Thermostat access singal aftCalcF. This might lead to undefined behavior.
        // Therefore, we use other signals here, to make sure the Thermostat would be always called first, before force distributions take place.
        // connection to after _aftCalcF()
        //_aftCalcF = integrator->aftCalcF.connect(
        //        boost::bind(&Adress::aftCalcF, this));

        // connection to after aftCalcSlow
        _aftCalcSlow = integrator->aftCalcSlow.connect(
                boost::bind(&Adress::aftCalcF, this), boost::signals2::at_back);

        // connection to after _recalc2()
        _recalc2 = integrator->recalc2.connect(
                boost::bind(&Adress::aftCalcF, this), boost::signals2::at_front);

        // connection to after _befIntV()
        _befIntV = integrator->befIntV.connect(
                boost::bind(&Adress::aftCalcF, this), boost::signals2::at_front);
    }


    void Adress::SetPosVel(){

        System& system = getSystemRef();

        // Set the positions and velocity of CG particles & update weights.
        CellList localCells = system.storage->getLocalCells();
        for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {

              Particle &vp = *cit;

              FixedTupleListAdress::iterator it3;
              it3 = fixedtupleList->find(&vp);

              if (it3 != fixedtupleList->end()) {

                  std::vector<Particle*> atList;
                  atList = it3->second;

                  // Compute center of mass
                  Real3D cmp(0.0, 0.0, 0.0); // center of mass position
                  Real3D cmv(0.0, 0.0, 0.0); // center of mass velocity
                  for (std::vector<Particle*>::iterator it2 = atList.begin();
                                       it2 != atList.end(); ++it2) {
                      Particle &at = **it2;
                      cmp += at.mass() * at.position();
                      cmv += at.mass() * at.velocity();
                  }
                  cmp /= vp.getMass();
                  cmv /= vp.getMass();

                  // update (overwrite) the position and velocity of the VP
                  vp.position() = cmp;
                  vp.velocity() = cmv;

                  if (KTI == false) {
                      // calculate distance to nearest adress particle or center
                      std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
                      Real3D pa = **it2; // position of adress particle
                      Real3D d1(0.0, 0.0, 0.0);
                      verletList->getSystem()->bc->getMinimumImageVector(d1, vp.position(), pa);
                      real min1sq;
                      if (verletList->getAdrRegionType()) { // spherical adress region
                        min1sq = d1.sqr(); // set min1sq before loop
                        ++it2;
                        for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                             pa = **it2;
                             verletList->getSystem()->bc->getMinimumImageVector(d1, vp.position(), pa);
                             real distsq1 = d1.sqr();
                             if (distsq1 < min1sq) min1sq = distsq1;
                        }
                      }
                      else { //slab-type adress region
                        min1sq = d1[0]*d1[0];   // set min1sq before loop
                        ++it2;
                        for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                             pa = **it2;
                             verletList->getSystem()->bc->getMinimumImageVector(d1, vp.position(), pa);
                             real distsq1 = d1[0]*d1[0];
                             if (distsq1 < min1sq) min1sq = distsq1;

                        }
                      }

                      real w = weight(min1sq);
                      vp.lambda() = w;

                      real wDeriv = weightderivative(min1sq);
                      vp.lambdaDeriv() = wDeriv;

                  }

              }
              else { // this should not happen
                  std::cout << " VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
                  std::cout << " (" << vp.position() << ")\n";
                  exit(1);
                  return;
              }


        }

    }


    void Adress::initForces(){

        System& system = getSystemRef();

        // AT reals
        ParticleList& adrATparticles = system.storage->getAdrATParticles();
        for (std::vector<Particle>::iterator it = adrATparticles.begin();
                it != adrATparticles.end(); ++it) {
            it->force() = 0.0;
        }

        // AT ghosts
        typedef std::list<ParticleList> ParticleListAdr;
        ParticleListAdr& adrATparticlesG = system.storage->getAdrATParticlesG();
        for (ParticleListAdr::iterator it = adrATparticlesG.begin();
                it != adrATparticlesG.end(); ++it) {

            for (ParticleList::iterator it2 = it->begin();
                    it2 != it->end(); ++it2) {
                it2->force() = 0.0;
                it2->drift() = 0.0;
            }
        }
    }



    void Adress::integrate1(real& maxSqDist){

        System& system = getSystemRef();
        real dt = integrator->getTimeStep();

        ParticleList& adrATparticles = system.storage->getAdrATParticles();
        for (std::vector<Particle>::iterator it = adrATparticles.begin();
                it != adrATparticles.end(); it++) {

            real sqDist = 0.0;
            real dtfm = 0.5 * dt / it->mass();

            // Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t)
            it->velocity() += dtfm * it->force();

            // Propagate positions (only NVT): p(t + dt) = p(t) + dt * v(t+0.5*dt)
            Real3D deltaP = dt * it->velocity();
            it->position() += deltaP;
            sqDist += deltaP * deltaP;

            maxSqDist = std::max(maxSqDist, sqDist);
        }

        // Set the positions and velocity of CG particles
        CellList localCells = system.storage->getLocalCells();
        for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {


              Particle &vp = *cit;

              FixedTupleListAdress::iterator it3;
              it3 = fixedtupleList->find(&vp);

              if (it3 != fixedtupleList->end()) {

                  std::vector<Particle*> atList;
                  atList = it3->second;

                  // Compute center of mass
                  Real3D cmp(0.0, 0.0, 0.0); // center of mass position
                  Real3D cmv(0.0, 0.0, 0.0); // center of mass velocity
                  for (std::vector<Particle*>::iterator it2 = atList.begin();
                                       it2 != atList.end(); ++it2) {
                      Particle &at = **it2;
                      cmp += at.mass() * at.position();
                      cmv += at.mass() * at.velocity();
                  }
                  cmp /= vp.getMass();
                  cmv /= vp.getMass();

                  // update (overwrite) the position and velocity of the VP
                  vp.position() = cmp;
                  vp.velocity() = cmv;

              }
              else { // this should not happen
                  std::cout << " VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
                  std::cout << " (" << vp.position() << ")\n";
                  exit(1);
                  return;
              }
        }

        // Communicate new position of region defining particles
        communicateAdrPositions();

        // Update resolution values if KTI == false
        if (KTI == false) {

          CellList localCells2 = system.storage->getLocalCells();
          for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {


                Particle &vp = *cit;

                FixedTupleListAdress::iterator it3;
                it3 = fixedtupleList->find(&vp);

                if (it3 != fixedtupleList->end()) {

                        // calculate distance to nearest adress particle or center
                        std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
                        Real3D pa = **it2; // position of adress particle
                        Real3D d1(0.0, 0.0, 0.0);
                        real min1sq;
                        verletList->getSystem()->bc->getMinimumImageVector(d1, vp.position(), pa);
                        if (verletList->getAdrRegionType()) { // spherical adress region
                          min1sq = d1.sqr(); // set min1sq before loop
                          ++it2;
                          for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                               pa = **it2;
                               verletList->getSystem()->bc->getMinimumImageVector(d1, vp.position(), pa);
                               real distsq1 = d1.sqr();
                               if (distsq1 < min1sq) min1sq = distsq1;
                          }
                        }
                        else { //slab-type adress region
                          min1sq = d1[0]*d1[0];   // set min1sq before loop
                          ++it2;
                          for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                               pa = **it2;
                               verletList->getSystem()->bc->getMinimumImageVector(d1, vp.position(), pa);
                               real distsq1 = d1[0]*d1[0];
                               if (distsq1 < min1sq) min1sq = distsq1;
                          }
                        }


                        real w = weight(min1sq);
                        vp.lambda() = w;

                        real wDeriv = weightderivative(min1sq);
                        vp.lambdaDeriv() = wDeriv;

                        // This loop is required when applying routines which use atomistic lambdas.
                        /*for (std::vector<Particle*>::iterator it2 = atList.begin();
                                         it2 != atList.end(); ++it2) {
                            Particle &at = **it2;
                            at.lambda() = vp.lambda();
                            at.lambdaDeriv() = vp.lambdaDeriv();
                        }*/

                }
                else { // this should not happen
                    std::cout << " VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
                    std::cout << " (" << vp.position() << ")\n";
                    exit(1);
                    return;
                }

          }

        }

    }


    void Adress::integrate2() {

        System& system = getSystemRef();
        real dt = integrator->getTimeStep();

        // propagete real AT particles
        ParticleList& adrATparticles = system.storage->getAdrATParticles();
        for (std::vector<Particle>::iterator it = adrATparticles.begin();
                it != adrATparticles.end(); ++it) {

            real dtfm = 0.5 * dt / it->mass();

            // Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t)
            it->velocity() += dtfm * it->force();
        }

        //Update CG velocities
        CellList localCells = system.storage->getLocalCells();
        for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {

              Particle &vp = *cit;

              FixedTupleListAdress::iterator it3;
              it3 = fixedtupleList->find(&vp);

              if (it3 != fixedtupleList->end()) {

                  std::vector<Particle*> atList;
                  atList = it3->second;

                  Real3D cmv(0.0, 0.0, 0.0); // center of mass velocity
                  for (std::vector<Particle*>::iterator it2 = atList.begin();
                                       it2 != atList.end(); ++it2) {
                      Particle &at = **it2;
                      cmv += at.mass() * at.velocity();
                  }
                  cmv /= vp.getMass();
                  vp.velocity() = cmv;

              }
              else { // this should not happen
                  std::cout << " VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
                  std::cout << " (" << vp.position() << ")\n";
                  exit(1);
                  return;
              }


        }



    }


    void Adress::integrateSlow() {

        System& system = getSystemRef();
        real dt = integrator->getTimeStep();

        // propagete real AT particles
        ParticleList& adrATparticles = system.storage->getAdrATParticles();
        for (std::vector<Particle>::iterator it = adrATparticles.begin();
                it != adrATparticles.end(); ++it) {

            real dtfm = 0.5 * multistep * dt / it->mass();

            // Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t)
            it->velocity() += dtfm * it->force();
        }

        //Update CG velocities
        CellList localCells = system.storage->getLocalCells();
        for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {

              Particle &vp = *cit;

              FixedTupleListAdress::iterator it3;
              it3 = fixedtupleList->find(&vp);

              if (it3 != fixedtupleList->end()) {

                  std::vector<Particle*> atList;
                  atList = it3->second;

                  Real3D cmv(0.0, 0.0, 0.0); // center of mass velocity
                  for (std::vector<Particle*>::iterator it2 = atList.begin();
                                       it2 != atList.end(); ++it2) {
                      Particle &at = **it2;
                      cmv += at.mass() * at.velocity();
                  }
                  cmv /= vp.getMass();
                  vp.velocity() = cmv;

              }
              else { // this should not happen
                  std::cout << " VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
                  std::cout << " (" << vp.position() << ")\n";
                  exit(1);
                  return;
              }


        }



    }

    void Adress::communicateAdrPositions(){
       //if adrCenter is not set, the center of adress zone moves along with some particles
       //the coordinates of the center(s) (adrPositions) must be communicated to all nodes

      //As this is a bit complicated: When upating only every other step, it's important to use and communicate not the pointer to the actual region defining particle's position, but to a copy of it. Otherwise, while the other CPUs would create new copies that don't move in between the communication, the pointer on the original CPU would still point to the original particle, which keeps moving. Hence, we create a totally new position vector to hold copies of the positions. It has to keep existing outside of this function, hence static. Furthermore, we want to reserve the memory once at the beginning, as soon as we know the what to reserve. Hence, the firstcall construct. The reservation is necessary, as the vector elements are not supposed to move in memory during push_back.

      if (!(verletList->getAdrCenterSet())) {

        static bool firstcall = true;
        static std::vector<Real3D> adrposlist;

        if((firstcall) && (verletList->getAdrList().size() > 0)){
          adrposlist.reserve(verletList->getAdrList().size());
          firstcall = false;
        }

        if(updatecount == 0){

          // Clear the vector of copies
          adrposlist.clear();

          if(verletList->getAdrList().size() == 1){

            // Old Version (works only for single moving region but is faster than below)
            CellList realcells = getSystem()->storage->getRealCells(); //should be realcells
            int root,mayberoot,lroot;
            lroot=0;
            verletList->adrPositions.clear();
                for (CellListIterator it(realcells); it.isValid(); ++it) {
                    if (verletList->getAdrList().count(it->id()) == 1) {

                        // Update the copy and append address to adrPositions
                        adrposlist.push_back(it->position());
                        verletList->adrPositions.push_back(&(adrposlist.back()));

                        //verletList->adrPositions.push_back(&(it->position())); // without the additional copy, old version
                        lroot = getSystem()->comm->rank();
                    }
            }
            mayberoot = (lroot ? lroot : 0);
            boost::mpi::all_reduce(*getSystem()->comm,mayberoot,root,boost::mpi::maximum<int>());
            mpi::broadcast(*getSystem()->comm,verletList->adrPositions,root);

          }
          else if(verletList->getAdrList().size() > 1){

            // New version (works also for many moving regions but is slower)
            System& system = getSystemRef();
            std::vector<Real3D*> procAdrPositions;
            CellList realcells = getSystem()->storage->getRealCells();

            verletList->adrPositions.clear();
            for (CellListIterator it(realcells); it.isValid(); ++it) {

                if (verletList->getAdrList().count(it->id()) == 1) {

                   // Update the copy and append address to adrPositions
                   adrposlist.push_back(it->position());
                   procAdrPositions.push_back(&(adrposlist.back()));

		        //procAdrPositions.push_back(&(it->position())); // without the additional copy, old version
                }

            }

            std::vector<std::vector<Real3D*> > procAdrPositionsAll;
            boost::mpi::all_gather(*system.comm, procAdrPositions, procAdrPositionsAll);

            for (std::vector<std::vector<Real3D*> >::iterator itr=procAdrPositionsAll.begin(); itr != procAdrPositionsAll.end(); ++itr) {
                for (std::vector<Real3D*>::iterator itr2=(*itr).begin(); itr2 != (*itr).end(); ++itr2){
                  verletList->adrPositions.push_back(*itr2);
                }
            }

          }
          else{
              std::cout << "When using moving AdResS regions there can be only either one or several of them. It seems like the list of PIDs of region defining particles is empty for some reason.\n";
              exit(1);
              return;
          }

        }

        updatecount += 1;
        if(updatecount == regionupdates){updatecount = 0;}

      }

    }



    // AdResS Weighting function
    real Adress::weight(real distanceSqr){
        if (dex2 > distanceSqr) return 1.0;
        else if (dexdhy2 < distanceSqr) return 0.0;
        else {
            real argument = sqrt(distanceSqr) - dex;
            //return 1.0-(30.0/(pow(dhy, 5.0)))*(1.0/5.0*pow(argument, 5.0)-dhy/2.0*pow(argument, 4.0)+1.0/3.0*pow(argument, 3.0)*dhy*dhy);
            return pow(cos(pidhy2 * argument),2.0); // for cosine squared weighting function
        }
    }
    real Adress::weightderivative(real distanceSqr){
        if (dex2 > distanceSqr) return 0.0;
        else if (dexdhy2 < distanceSqr) return 0.0;
        else{
          real argument = sqrt(distanceSqr) - dex;
          //return -(30.0/(pow(dhy, 5.0)))*(pow(argument, 4.0)-2.0*dhy*pow(argument, 3.0)+argument*argument*dhy*dhy);
          return -pidhy2 * 2.0 * cos(pidhy2*argument) * sin(pidhy2*argument); // for cosine squared weighting function
        }
    }



    void Adress::aftCalcF(){
        System& system = getSystemRef();
        CellList localCells = system.storage->getLocalCells();
        for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {

        Particle &vp = *cit;

        FixedTupleListAdress::iterator it3;
        it3 = fixedtupleList->find(&vp);

        if (it3 != fixedtupleList->end()) {

            std::vector<Particle*> atList;
            atList = it3->second;

            // update force of AT particles belonging to a VP
            Real3D vpfm = vp.force() / vp.getMass();
            for (std::vector<Particle*>::iterator it2 = atList.begin();
                                 it2 != atList.end(); ++it2) {
                Particle &at = **it2;

                at.force() += at.mass() * vpfm;
            }
        }
        else { // this should not happen
            std::cout << " particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
            std::cout << " (" << vp.position() << ")\n";
            exit(1);
            return;
        }
      }

    }



    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void Adress::registerPython() {
      using namespace espressopp::python;

      class_<Adress, shared_ptr<Adress>, bases<Extension> >
        ("integrator_Adress", init<shared_ptr<System>, shared_ptr<VerletListAdress>, shared_ptr<FixedTupleListAdress>, bool, int, int >())
        .def("connect", &Adress::connect)
        .def("disconnect", &Adress::disconnect)
        ;
    }

  }

}
