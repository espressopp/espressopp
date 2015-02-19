/*
  Copyright (C) 2012,2013
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

namespace espresso {

  namespace integrator {

    using namespace espresso::iterator;

    Adress::Adress(shared_ptr<System> _system, shared_ptr<VerletListAdress> _verletList, shared_ptr<FixedTupleListAdress> _fixedtupleList, bool _KTI /*= false*/)
        : Extension(_system), verletList(_verletList), fixedtupleList(_fixedtupleList), KTI(_KTI){
        LOG4ESPP_INFO(theLogger, "construct Adress");
        type = Extension::Adress;
        
        // AdResS stuff
        dhy = verletList->getHy();
        pidhy2 = M_PI/(dhy * 2.0);
        dex = verletList->getEx();
        dex2 = dex * dex;
        dexdhy = dex + verletList->getHy();
        dexdhy2 = dexdhy * dexdhy;

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
        _inIntP.disconnect();
        //_aftCalcF.disconnect();
        _recalc2.disconnect();
        _befIntV.disconnect();
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

        // connection to inside of integrate1()
        _inIntP = integrator->inIntP.connect(
                boost::bind(&Adress::communicateAdrPositions, this), boost::signals2::at_front);

        // connection to after integrate2()
        _integrate2 = integrator->aftIntV.connect(
                boost::bind(&Adress::integrate2, this), boost::signals2::at_front);
        
        // Note: Both this extension as well as Langevin Thermostat access singal aftCalcF. This might lead to undefined behavior.
        // Therefore, we use other signals here, to make sure the Thermostat would be always called first, before force distributions take place.
        // connection to after _aftCalcF()
        //_aftCalcF = integrator->aftCalcF.connect(
        //        boost::bind(&Adress::aftCalcF, this));        
        
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
                  //real M = vp.getMass(); // sum of mass of AT particles
                  for (std::vector<Particle*>::iterator it2 = atList.begin();
                                       it2 != atList.end(); ++it2) {
                      Particle &at = **it2;
                      //Real3D d1 = at.position() - vp.position();
                      //Real3D d1;
                      //verletList->getSystem()->bc->getMinimumImageVectorBox(d1, at.position(), vp.position());
                      //cmp += at.mass() * d1;

                      cmp += at.mass() * at.position();
                      cmv += at.mass() * at.velocity();
                  }
                  cmp /= vp.getMass();
                  cmv /= vp.getMass();
                  //cmp += vp.position(); // cmp is a relative position
                  //std::cout << " cmp M: "  << M << "\n\n";
                  //std::cout << "  moving VP to " << cmp << ", velocitiy is " << cmv << "\n";

                  // update (overwrite) the position and velocity of the VP
                  vp.position() = cmp;
                  vp.velocity() = cmv;

                  if (KTI == false) {
                      // calculate distance to nearest adress particle or center
                      std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
                      Real3D pa = **it2; // position of adress particle
                      Real3D d1(0.0, 0.0, 0.0);
                      //Real3D d1 = vp.position() - pa;                                                      // X SPLIT VS SPHERE CHANGE
                      verletList->getSystem()->bc->getMinimumImageVector(d1, vp.position(), pa);
                      real min1sq;
                      //real d1 = vp.position()[0] - pa[0];                                                // X SPLIT VS SPHERE CHANGE
                      //real min1sq = d1.sqr();  // set min1sq before loop                                   // X SPLIT VS SPHERE CHANGE
                      //real min1sq = d1*d1;   // set min1sq before loop                                   // X SPLIT VS SPHERE CHANGE
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
                        min1sq = d1[0]*d1[0];   // set min1sq before loop                                   // X SPLIT VS SPHERE CHANGE
                        ++it2;
                        for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                             pa = **it2;
                             //d1 = vp.position() - pa;                                                          // X SPLIT VS SPHERE CHANGE
                             //d1 = vp.position()[0] - pa[0];
                             verletList->getSystem()->bc->getMinimumImageVector(d1, vp.position(), pa);        // X SPLIT VS SPHERE CHANGE
                             //real distsq1 = d1.sqr();                                                          // X SPLIT VS SPHERE CHANGE
                             real distsq1 = d1[0]*d1[0];                                                           // X SPLIT VS SPHERE CHANGE
                             //std::cout << pa << " " << sqrt(distsq1) << "\n";
                             if (distsq1 < min1sq) min1sq = distsq1;
                      
                        }
                      }

                      real w = weight(min1sq);                  
                      vp.lambda() = w;                  
                      //weights.insert(std::make_pair(&vp, w));

                      real wDeriv = weightderivative(sqrt(min1sq));
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
            
            //if(it->id()==2135){
            //   std::cout << "Force of atomistic particle (AdResS. sim.) with id " << it->id() << " is: " << std::setprecision(15) << it->force() << "\n";  // FOR DEBUGGING
            //}
            
            //std::cout << "Force of atomistic particle (AdResS. sim.) with id " << it->id() << " is: " << std::setprecision(15) << it->force() << "\n";  // FOR DEBUGGING
            //std::cout << "Position of atomistic particle (AdResS. sim.) with id " << it->id() << " is: " << std::setprecision(15) << it->position() << "\n";
            
            real sqDist = 0.0;
            real dtfm = 0.5 * dt / it->mass();

            // Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t)
            it->velocity() += dtfm * it->force();

            // Propagate positions (only NVT): p(t + dt) = p(t) + dt * v(t+0.5*dt)
            Real3D deltaP = dt * it->velocity();
            //std::cout << it->id() << ": from (" << it->position() << ")";
            it->position() += deltaP;
            sqDist += deltaP * deltaP;
            //std::cout << " to (" << it->position() << ") " << sqrt(sqDist) << "\n";

            maxSqDist = std::max(maxSqDist, sqDist);
        }               
        
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
                  //real M = vp.getMass(); // sum of mass of AT particles
                  for (std::vector<Particle*>::iterator it2 = atList.begin();
                                       it2 != atList.end(); ++it2) {
                      Particle &at = **it2;
                      //Real3D d1 = at.position() - vp.position();
                      //Real3D d1;
                      //verletList->getSystem()->bc->getMinimumImageVectorBox(d1, at.position(), vp.position());
                      //cmp += at.mass() * d1;

                      cmp += at.mass() * at.position();
                      cmv += at.mass() * at.velocity();
                  }
                  cmp /= vp.getMass();
                  cmv /= vp.getMass();
                  //cmp += vp.position(); // cmp is a relative position
                  //std::cout << " cmp M: "  << M << "\n\n";
                  //std::cout << "  moving VP to " << cmp << ", velocitiy is " << cmv << "\n";

                  // update (overwrite) the position and velocity of the VP
                  vp.position() = cmp;
                  vp.velocity() = cmv;

                  if (KTI == false) {
                  
                      // calculate distance to nearest adress particle or center
                      std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
                      Real3D pa = **it2; // position of adress particle
                      Real3D d1(0.0, 0.0, 0.0);
                      real min1sq;
                      //Real3D d1 = vp.position() - pa;                                                      // X SPLIT VS SPHERE CHANGE
                      verletList->getSystem()->bc->getMinimumImageVector(d1, vp.position(), pa);
                      //real d1 = vp.position()[0] - pa[0];                                                // X SPLIT VS SPHERE CHANGE
                      //real min1sq = d1.sqr();  // set min1sq before loop                                   // X SPLIT VS SPHERE CHANGE
                      //real min1sq = d1*d1;   // set min1sq before loop                                   // X SPLIT VS SPHERE CHANGE
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
                        min1sq = d1[0]*d1[0];   // set min1sq before loop                                   // X SPLIT VS SPHERE CHANGE
                        ++it2;
                        for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                             pa = **it2;
                             //d1 = vp.position() - pa;                                                          // X SPLIT VS SPHERE CHANGE
                             //d1 = vp.position()[0] - pa[0];
                             verletList->getSystem()->bc->getMinimumImageVector(d1, vp.position(), pa);        // X SPLIT VS SPHERE CHANGE
                             //real distsq1 = d1.sqr();                                                          // X SPLIT VS SPHERE CHANGE
                             real distsq1 = d1[0]*d1[0];                                                           // X SPLIT VS SPHERE CHANGE
                             //std::cout << pa << " " << sqrt(distsq1) << "\n";
                             if (distsq1 < min1sq) min1sq = distsq1;
                        }
                      }


                      real w = weight(min1sq);                  
                      vp.lambda() = w;                  
                      //weights.insert(std::make_pair(&vp, w));

                      real wDeriv = weightderivative(sqrt(min1sq));
                      vp.lambdaDeriv() = wDeriv;
                      
                      // This loop is required when applying routines which use atomistic lambdas.
                      /*for (std::vector<Particle*>::iterator it2 = atList.begin();
                                       it2 != atList.end(); ++it2) {
                          Particle &at = **it2;
                          at.lambda() = vp.lambda();
                          at.lambdaDeriv() = vp.lambdaDeriv();
                      }*/
                  
                  }
                  
              }
              else { // this should not happen
                  std::cout << " VP particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
                  std::cout << " (" << vp.position() << ")\n";
                  exit(1);
                  return;
              }
            
            
        } 

        //std::cout << " " << maxSqDist << "\n";
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

                  // Compute center of mass
                  //Real3D cmp(0.0, 0.0, 0.0); // center of mass position
                  Real3D cmv(0.0, 0.0, 0.0); // center of mass velocity
                  //real M = vp.getMass(); // sum of mass of AT particles
                  for (std::vector<Particle*>::iterator it2 = atList.begin();
                                       it2 != atList.end(); ++it2) {
                      Particle &at = **it2;
                      //Real3D d1 = at.position() - vp        shared_ptr<VerletListAdress> verletList;
                      //Real3D d1;
                      //verletList->getSystem()->bc->getMinimumImageVectorBox(d1, at.position(), vp.position());
                      //cmp += at.mass() * d1;

                      //cmp += at.mass() * at.position();
                      cmv += at.mass() * at.velocity();
                  }
                  //cmp /= vp.getMass();
                  cmv /= vp.getMass();
                  //cmp += vp.position(); // cmp is a relative position
                  //std::cout << " cmp M: "  << M << "\n\n";
                  //std::cout << "  moving VP to " << cmp << ", velocitiy is " << cmv << "\n";

                  // update (overwrite) the position and velocity of the VP
                  //vp.position() = cmp;
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

       // get local cells
       CellList realcells = getSystem()->storage->getRealCells(); //should be realcells
       int root,mayberoot,lroot;
       lroot=0;
       if (!(verletList->getAdrCenterSet())) {
          verletList->adrPositions.clear(); // clear position pointers
          for (CellListIterator it(realcells); it.isValid(); ++it) {
              if (verletList->getAdrList().count(it->id()) == 1) {
                  verletList->adrPositions.push_back(&(it->position()));
                  lroot = getSystem()->comm->rank(); //for the moment only works when there's only one particle in adrPositions
              }
              //TODO if length(adrPositions) > 1 print warning
          }
          mayberoot = (lroot ? lroot : 0);
          boost::mpi::all_reduce(*getSystem()->comm,mayberoot,root,boost::mpi::maximum<int>());
          mpi::broadcast(*getSystem()->comm,verletList->adrPositions,root); // only necessary for moving adrCenter
       }
    }
    
    
    
    // AdResS Weighting function
    real Adress::weight(real distanceSqr){
        if (dex2 > distanceSqr) return 1.0;
        else if (dexdhy2 < distanceSqr) return 0.0;
        else {
            real argument = sqrt(distanceSqr) - dex;
            return 1.0-(30.0/(pow(dhy, 5.0)))*(1.0/5.0*pow(argument, 5.0)-dhy/2.0*pow(argument, 4.0)+1.0/3.0*pow(argument, 3.0)*dhy*dhy);
            //return pow(cos(pidhy2 * argument),2.0); // for cosine squared weighting function
        }
    }
    real Adress::weightderivative(real distance){
        real argument = distance - dex;
        return -(30.0/(pow(dhy, 5.0)))*(pow(argument, 4.0)-2.0*dhy*pow(argument, 3.0)+argument*argument*dhy*dhy);
        //return -pidhy2 * 2.0 * cos(pidhy2*argument) * sin(pidhy2*argument); // for cosine squared weighting function
    }

    
    
    void Adress::aftCalcF(){        
        System& system = getSystemRef();
        CellList localCells = system.storage->getLocalCells();
        for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
        /*for (std::set<Particle*>::iterator it=adrZone.begin();
                it != adrZone.end(); ++it) {*/

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

                //vp.force() +=  (vp.getMass() * at.force()) / (3.0 * at.mass());
                
                at.force() += at.mass() * vpfm;
                //std::cout << "Force of atomistic particle (AdResS sim.) with id " << at.id() << " is: " << at.force() << "\n";
            }
        }
        else { // this should not happen
            std::cout << " particle " << vp.id() << "-" << vp.ghost() << " not found in tuples ";
            std::cout << " (" << vp.position() << ")\n";
            exit(1);
            return;
        }
      }
      
      /*for (std::set<Particle*>::iterator it=cgZone.begin();
                    it != cgZone.end(); ++it) {

            Particle &vp = **it;

            FixedTupleListAdress::iterator it3;
            it3 = fixedtupleList->find(&vp);

            if (it3 != fixedtupleList->end()) {

                std::vector<Particle*> atList1;
                atList1 = it3->second;

                Real3D vpfm = vp.force() / vp.getMass();
                for (std::vector<Particle*>::iterator itv = atList1.begin();
                        itv != atList1.end(); ++itv) {
                    Particle &at = **itv;
                    
                    //vp.force() +=  (vp.getMass() * at.force()) / (3.0 * at.mass());
                    
                    // at.velocity() = vp.velocity(); // overwrite velocity
                    at.force() += at.mass() * vpfm;
                    //std::cout << "f" << at.mass() * vpfm << " m " << at.mass() << " M "<<  vp.getMass() << " id " << at.id() << std::endl;
                }

            }
            else { // this should not happen
                std::cout << " VP particle not found in tuples: " << vp.id() << "-" << vp.ghost();
                exit(1);
                return;
            }
      }*/
    }
    
    
    
    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void Adress::registerPython() {
      using namespace espresso::python;

      class_<Adress, shared_ptr<Adress>, bases<Extension> >
        ("integrator_Adress", init<shared_ptr<System>, shared_ptr<VerletListAdress>, shared_ptr<FixedTupleListAdress>, bool >())
        .def("connect", &Adress::connect)
        .def("disconnect", &Adress::disconnect)
        ;
    }

  }

}
