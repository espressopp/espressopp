/*
  Copyright (C) 2016
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

#include "python.hpp"

#include "Rattle.hpp"

#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "System.hpp"
#include "bc/BC.hpp"

//Refs: 
//Andersen, H. C. ``Rattle: A ``velocity'' version of the Shake algorithm for molecular dynamics calculations'', J. Comp. Physics, 52, 24-34 (1983)
//Allen & Tildesley, ``Computer Simulation of Liquids'', OUP, 1987
//Tip: to see how to add the calculation of kinetic energy and constraint virial, see Allen&Tildesley file f09-rattle.f77

namespace espressopp {
  using namespace iterator;
  namespace integrator {

    LOG4ESPP_LOGGER(Rattle::theLogger, "Rattle");

    Rattle::Rattle(shared_ptr<System> _system, 
        real _maxit, real _tol, real _rptol)
    : Extension(_system),
        maxit(_maxit), tol(_tol), rptol(_rptol) {

        LOG4ESPP_INFO(theLogger, "construct Rattle");

    }

    Rattle::~Rattle() {
      LOG4ESPP_INFO(theLogger, "~Rattle");
    }

    void Rattle::disconnect(){
      _befIntP.disconnect();
      _aftIntP.disconnect();
      _aftIntV.disconnect();
    }

    void Rattle::connect(){
      _befIntP  = integrator->befIntP.connect( boost::bind(&Rattle::saveOldPos, this));
      _aftIntP  = integrator->aftIntP.connect( boost::bind(&Rattle::applyPositionConstraints, this));
      _aftIntV  = integrator->aftIntV.connect( boost::bind(&Rattle::applyVelocityConstraints, this));
    }

    void Rattle::addBond(int pid1, int pid2, real constraintDist, real mass1, real mass2) {
      if (mass2 > mass1) {
        std::ostringstream msg;
        msg << "In Rattle, the heavy atom should be listed before the hydrogen in each constrained bond" << std::endl;
        throw std::runtime_error( msg.str() );
      }
      ConstrainedBond newbond;
      newbond.pidHeavy = pid1;
      newbond.pidHyd = pid2;
      newbond.constraintDist2 = constraintDist*constraintDist;
      newbond.invmassHeavy = 1.0/mass1;
      newbond.invmassHyd = 1.0/mass2;
      constrainedBonds.insert(std::make_pair(pid2,newbond));
      constrainedBondsKeys.push_back(pid2);
    }

    void Rattle::saveOldPos() {
      oldPos.clear();
      System& system = getSystemRef();
      boost::unordered_map<longint, ConstrainedBond>::iterator it;
      //collect coordinates of bonds on this CPU
      //loop over bonds
      for (it=constrainedBonds.begin(); it != constrainedBonds.end(); it++ ) {
        longint heavyID = (it->second).pidHeavy;
        longint lightID = (it->second).pidHyd;
        Particle* hp = system.storage->lookupAdrATParticle(heavyID);
        if (hp) { //heavy atom is on this node
          Particle* lp = system.storage->lookupAdrATParticle(lightID);
          if (!lp) {
            std::ostringstream msg;
            msg << "In Rattle, cannot find light particle " << lightID << ", all light and heavy particles in a group of rigid bonds must be on the same node" << std::endl;
            throw std::runtime_error( msg.str() );
          }
          oldPos.insert(std::make_pair(lightID,
                        std::pair<Real3D,Real3D>(hp->getPos(),lp->getPos())));
        }
      }
    }

    void Rattle::applyPositionConstraints() {

      int iteration = 0;
      bool done = false;
      real dt = integrator->getTimeStep();
      System& system = getSystemRef();
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions

      boost::unordered_map<longint, Real3D> currPosition;
      boost::unordered_map<longint, Real3D> currVelocity;
      boost::unordered_map<longint, bool> movedLastTime; //was this particle moved last time?
      boost::unordered_map<longint, bool> movingThisTime; //is the particle being moved this time?

      if (oldPos.size() == 0) {return;} //no rigid bonds on this node

      //loop over bonds on this node
      OldPos::iterator it;
      for (it=oldPos.begin(); it != oldPos.end(); it++ ) { 
        longint pidHyd = it->first;
        longint pidHeavy = constrainedBonds[pidHyd].pidHeavy;
        //In the maps below, the values for each key pidHeavy may be initialised more than, if pidHeavy is involved in more than one constrained bond. This is not a problem
        movedLastTime[pidHyd] = true;
        movedLastTime[pidHeavy] = true;
        movingThisTime[pidHyd] = false;
        movingThisTime[pidHeavy] = false;
        currPosition[pidHyd] = system.storage->lookupAdrATParticle(pidHyd)->getPos();
        currPosition[pidHeavy] = system.storage->lookupAdrATParticle(pidHeavy)->getPos(); 
        currVelocity[pidHyd] = system.storage->lookupAdrATParticle(pidHyd)->getV();
        currVelocity[pidHeavy] = system.storage->lookupAdrATParticle(pidHeavy)->getV(); 
      }

      //constraint interations
      while (!done && iteration < maxit) {
        done = true;
        //loop over constrained bonds on this cpu
        for (it=oldPos.begin(); it != oldPos.end(); it++) {
          //get pids of the two particles in this bond
          longint a = it->first; //pidHyd
          longint b = constrainedBonds[a].pidHeavy;
          if (movedLastTime[a] || movedLastTime[b]) { 
            //compare current distance to desired constraint distance
            Real3D pab;
            bc.getMinimumImageVectorBox(pab,currPosition[a],currPosition[b]); //a-b, current positions which change during iterations
            real pabsq = pab.sqr();
            real constraint_absq = constrainedBonds[a].constraintDist2;
            real diffsq = pabsq - constraint_absq; 
            if (fabs(diffsq) > (constraint_absq*tol) ) {
              //get ab vector before unconstrained position update
              Real3D rab;
              bc.getMinimumImageVectorBox(rab,oldPos[a].second,oldPos[a].first); //pos at time t (end of last timestep), a-b;
              real rab_dot_pab = rab * pab; //r_ab(t) * r_ab,curr(t+dt)
              if (rab_dot_pab < (constraint_absq*rptol)) { //i.e. if angle is too large
                std::ostringstream msg;
                msg << "Constraint failure in RATTLE" << std::endl;
                throw std::runtime_error( msg.str() );
              }
              real rma = constrainedBonds[a].invmassHyd;
              real rmb = constrainedBonds[a].invmassHeavy;
              real gab = diffsq / (2.0 * (rma + rmb) * rab_dot_pab);
              //direct constraint along bond vector at end of previous timestep
              Real3D displ = gab * rab; 
              currPosition[a] -= rma * displ;
              currPosition[b] += rmb * displ;

              displ /= dt;
              currVelocity[a] -= rma * displ;
              currVelocity[b] += rmb * displ;

              movingThisTime[a] = true;
              movingThisTime[b] = true;
              done = false;
            }
          }
        }
        boost::unordered_map<longint, bool>::iterator it3;
        for (it3 = movedLastTime.begin(); it3 != movedLastTime.end(); it3++) {
          longint id = it3->first;
          movedLastTime[id] = movingThisTime[id];
          movingThisTime[id] = false;
        }

        iteration += 1;
      }

      if (!done) {
        std::ostringstream msg;
        msg << "Too many position constraint iterations in Rattle" << std::endl;
        throw std::runtime_error( msg.str() );
      }

      //store new values for positions
      boost::unordered_map<longint, Real3D>::iterator it2;
      for (it2 = currPosition.begin(); it2 != currPosition.end(); it2++) {
        longint id = it2->first;
        Particle* p = system.storage->lookupAdrATParticle(id);
        p->position() = currPosition[id];
        p->velocity() = currVelocity[id];
      }
    }

    void Rattle::applyVelocityConstraints() {

      int iteration = 0;
      bool done = false;
      real dt = integrator->getTimeStep();
      System& system = getSystemRef();
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions

      boost::unordered_map<longint, Real3D> currPosition;
      boost::unordered_map<longint, Real3D> currVelocity;
      boost::unordered_map<longint, bool> changedLastTime; //was this particle velocity changed last time?
      boost::unordered_map<longint, bool> changingThisTime; //is the particle velocity being changed this time?

      //get all constrained bonds (light particles) on this CPU
      //can't use contents of oldPos to identify which constrained bonds are on this CPU because particles may have changed CPU since applyPositionConstraints()
      lightPart.clear();
      boost::unordered_map<longint, ConstrainedBond>::iterator it3;
      for (it3=constrainedBonds.begin(); it3 != constrainedBonds.end(); it3++ ) {
        longint lightID = (it3->second).pidHyd;
        Particle* lp = system.storage->lookupAdrATParticle(lightID);
        if (lp) { //light atom is on this node
          lightPart.push_back(lightID);
        }
      }

      //loop over constrained bonds on this node
      std::vector<int>::iterator it;
      for (it=lightPart.begin(); it != lightPart.end(); it++ ) {
        longint pidHyd = *it;
        longint pidHeavy = constrainedBonds[pidHyd].pidHeavy;
        //In the maps below, the values for each key pidHeavy may be initialised more than, if pidHeavy is involved in more than one constrained bond. This is not a problem
        changedLastTime[pidHyd] = true;
        changedLastTime[pidHeavy] = true; 
        changingThisTime[pidHyd] = false;
        changingThisTime[pidHeavy] = false; 
        currPosition[pidHyd] = system.storage->lookupAdrATParticle(pidHyd)->getPos();
        currPosition[pidHeavy] = system.storage->lookupAdrATParticle(pidHeavy)->getPos(); 
        currVelocity[pidHyd] = system.storage->lookupAdrATParticle(pidHyd)->getV();
        currVelocity[pidHeavy] = system.storage->lookupAdrATParticle(pidHeavy)->getV(); 
      }

      //constraint interations
      while (!done && iteration < maxit) {
        done = true;
        //loop over constrained bonds on this cpu
        for (it=lightPart.begin(); it != lightPart.end(); it++) {
          //get pids of the two particles in this bond
          longint a = *it; //pidHyd
          longint b = constrainedBonds[a].pidHeavy;
          if (changedLastTime[a] || changedLastTime[b]) { 
            Real3D vab = currVelocity[a] - currVelocity[b];
            Real3D rab;
            bc.getMinimumImageVectorBox(rab,currPosition[a],currPosition[b]);
            real rab_dot_vab = rab * vab;
            real rma = constrainedBonds[a].invmassHyd;
            real rmb = constrainedBonds[a].invmassHeavy;
            real constraint_absq = constrainedBonds[a].constraintDist2;
            real gab = -1.0 * rab_dot_vab / ( (rma + rmb) * constraint_absq);
            if (fabs(gab) > tol) {
              Real3D deltav = gab * rab;
              currVelocity[a] += rma * deltav;
              currVelocity[b] -= rmb * deltav;

              changingThisTime[a] = true;
              changingThisTime[b] = true;
              done = false;
            }
          }
        }

        boost::unordered_map<longint, bool>::iterator it3;
        for (it3 = changedLastTime.begin(); it3 != changedLastTime.end(); it3++) {
          longint id = it3->first;
          changedLastTime[id] = changingThisTime[id];
          changingThisTime[id] = false;
        }

        iteration += 1;
      }

      if (!done) {
        std::ostringstream msg;
        msg << "Too many velocity constraint iterations in Rattle" << std::endl;
        throw std::runtime_error( msg.str() );
      }

      //store new values for velocities
      boost::unordered_map<longint, Real3D>::iterator it2;
      for (it2 = currVelocity.begin(); it2 != currVelocity.end(); it2++) {
        longint id = it2->first;
        Particle* p = system.storage->lookupAdrATParticle(id);
        p->velocity() = currVelocity[id];
      }
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void Rattle::registerPython() {

      using namespace espressopp::python;

      class_<Rattle, shared_ptr<Rattle>, bases<Extension> >
        ("integrator_Rattle", init<shared_ptr<System>, real, real, real>())
         .def("addBond", &Rattle::addBond)
        ;
    }
  }
}
