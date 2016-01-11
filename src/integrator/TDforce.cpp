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
#include "TDforce.hpp"
#include "bc/BC.hpp"
#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "interaction/InterpolationLinear.hpp"
#include "interaction/InterpolationAkima.hpp"
#include "interaction/InterpolationCubic.hpp"


namespace espressopp {
  using namespace iterator;

  namespace integrator {

    LOG4ESPP_LOGGER(TDforce::theLogger, "TDforce");

    TDforce::TDforce(shared_ptr<System> system, shared_ptr<VerletListAdress> _verletList, real _startdist, real _enddist, int _edgeweightmultiplier)
    :Extension(system), verletList(_verletList), startdist(_startdist), enddist(_enddist), edgeweightmultiplier(_edgeweightmultiplier) {

        // startdist & enddist are the distances between which the TD force actually acts.
        // This is usually more or less the thickness of the hybrid region. However, the TD force sometimes is applied in a slighty wider area.

        type = Extension::FreeEnergyCompensation;

        if (verletList->getAdrCenterSet()) { //adress region centre is fixed in space
          center=**(verletList->getAdrPositions().begin());
        }
        else
        {
          center = (0.0,0.0,0.0);
        }

        LOG4ESPP_INFO(theLogger, "TDforce constructed");

        sphereAdr = verletList->getAdrRegionType(); //true=spherical, false=xslab
    }

    TDforce::~TDforce() {}



    void TDforce::connect(){
        _applyForce = integrator->aftCalcF.connect(
            boost::bind(&TDforce::applyForce, this));
    }

    void TDforce::disconnect(){
        _applyForce.disconnect();
    }

    //void TDforce::setAdrRegionType(bool _sphereAdr){
    //    sphereAdr = _sphereAdr;
    //}

    //bool TDforce::getAdrRegionType(){
    //    return sphereAdr;
    //}


    void TDforce::addForce(int itype, const char* _filename, int type) {
        boost::mpi::communicator world;
        filename = _filename;
        Table table;

        if (itype == 1) { // create a new InterpolationLinear
            table = make_shared <interaction::InterpolationLinear> ();
            table->read(world, _filename);
        }

        else if (itype == 2) { // create a new InterpolationAkima
            table = make_shared <interaction::InterpolationAkima> ();
            table->read(world, _filename);
        }

        else if (itype == 3) { // create a new InterpolationCubic
            table = make_shared <interaction::InterpolationCubic> ();
            table->read(world, _filename);
        }

        forces.insert(std::make_pair(type,table));
    }


    void TDforce::applyForce() {
          LOG4ESPP_DEBUG(theLogger, "apply TD force");

          System& system = getSystemRef();
          const bc::BC& bc = *getSystemRef().bc;  // boundary conditions

          // iterate over CG particles
          CellList cells = system.storage->getRealCells();
          for(CellListIterator cit(cells); !cit.isDone(); ++cit) {

              Table table;
              if (forces.find(cit->getType())!=forces.end()) { //because there may be CG particles to which TD force is not applied
                table = forces.find(cit->getType())->second; //TODO shouldn't do find twice
              }

              if (table) {




//                   if (!(verletList->getAdrCenterSet())) { //moving adress region
// //TODO if length(adrPositions) > 1 print warning (for the moment only works when there's only one particle in adrPositions)
//                     center=**(verletList->adrPositions.begin());
//                   }

//                   if (sphereAdr){ // spherical adress region
//                      // calculate distance from reference point
//                      Real3D dist3D;
//                      bc.getMinimumImageVectorBox(dist3D,cit->getPos(),center); //pos - center
//                      real dist = sqrt(dist3D.sqr());

//                      if (dist>0.0) {
//                        // read fforce from table
//                        real fforce = table->getForce(dist);
//                        fforce /= dist;

//                        // substract td force
//                        cit->force() += (dist3D * fforce);
//                      }
//                   } else {
//                      // use this if you need 1-dir force only!
//                      real d1 = cit->getPos()[0] - center[0];
//                      real d1abs = fabs(d1);
//                      real force = table->getForce(d1abs);
// //std::cout<<cit->id()<<" "<<cit->force()<<" "<<cit->getPos()[0]<<" "<<center[0]<<" "<<d1<<" "<<d1abs<<" "<<force<<std::endl;
//                      if (d1>0.0) {
//                        cit->force()[0] += force;
//                      } else {
//                        cit->force()[0] -= force;
//                      }
// //std::cout<<cit->id()<<" "<<cit->force()<<std::endl;
//                   }



                  if (!(verletList->getAdrCenterSet())) { // adress regions defined based on particles

                    if (sphereAdr){

                      // ONLY DEUGGING
                      //int count = 1;
                      // ONLY DEUGGING

                      real buffer = 0.000001;
                      real width = fabs(enddist - startdist + 2.0*buffer);  // width of region where TD force acts
                      std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();   // get positions
                      Real3D pa = **it2;
                      Real3D dist3D;
                      bc.getMinimumImageVectorBox(dist3D,cit->getPos(), pa); // calculate vector between particle and first center
                      //std::cout << "pa1: " << pa << "\n";
                      Real3D mindist3D = dist3D; // set mindist3D before loop

                      real dist3Dabs = sqrt(dist3D.sqr()); // calculate absolute distance

                      real weight = 0.0;  // initialize weight

                      // weighting scheme: only particles that are in a hybrid region contribute
                      if ((dist3Dabs < enddist + buffer) && (dist3Dabs > startdist - buffer)){
                        weight = 1.0 - (dist3Dabs - startdist - buffer) / width; // 0 at edge to CG region, 1 at edge to atomistic region

                        // Direction modification at the edges
                        weight = pow(weight, edgeweightmultiplier);

                        // std::cout << "dist3Dabs1: " << dist3Dabs << "\n";
                        // std::cout << "width1: " << width << "\n";
                        // std::cout << "buffer1: " << buffer << "\n";
                        // std::cout << "startdist1: " << startdist << "\n";
                        // std::cout << "enddist1: " << enddist << "\n";
                        // std::cout << "(dist3Dabs - startdist - buffer)1: " << (dist3Dabs - startdist - buffer) << "\n\n";
                      }
                      else{
                        weight = 0.0; // particle outside of hybrid region
                      }

                      // if ((cit->position()[0] < 1.50216) && (cit->position()[0] > 0.80216) && (cit->position()[1] > 2.5) && (cit->position()[1] < 3.5) && (cit->position()[2] > 2.5) && (cit->position()[2] < 3.5)){
                      //         std::cout << "dist3Dabs1: " << dist3Dabs << "\n";
                      //         std::cout << "weight1: " << weight << "\n";
                      // }
                      //real norm = weight; // normalization constant
                      Real3D direction = dist3D*weight/dist3Dabs; // normalized but weighted direction vector

                      // Loop over all other centers
                      ++it2;
                      for (; it2 != verletList->getAdrPositions().end(); ++it2) {

                            // ONLY DEUGGING
                            //count += 1;
                            // ONLY DEUGGING

                            pa = **it2;
                            //std::cout << "pa2: " << pa << "\n";
                            verletList->getSystem()->bc->getMinimumImageVector(dist3D, cit->getPos(), pa); // calculate vector between particle and other centers
                            if (dist3D.sqr() < mindist3D.sqr()) mindist3D = dist3D; // make it the minimum if shortest length so far

                            dist3Dabs = sqrt(dist3D.sqr()); // calculate absolute distance

                            // calculate weights again, as above
                            if ((dist3Dabs < enddist + buffer) && (dist3Dabs > startdist - buffer)){
                              weight = 1.0 - (dist3Dabs - startdist - buffer) / width;

                              // std::cout << "dist3Dabs2: " << dist3Dabs << "\n";
                              // std::cout << "width2: " << width << "\n";
                              // std::cout << "buffer2: " << buffer << "\n";
                              // std::cout << "startdist2: " << startdist << "\n";
                              // std::cout << "enddist2: " << enddist << "\n";
                              // std::cout << "(dist3Dabs - startdist - buffer)2: " << (dist3Dabs - startdist - buffer) << "\n\n";

                              // Direction modification at the edges (TEST)
                              weight = pow(weight, edgeweightmultiplier);

                            }
                            else{
                              weight = 0.0;
                            }

                            // if ((cit->position()[0] < 1.50216) && (cit->position()[0] > 0.80216) && (cit->position()[1] > 2.5) && (cit->position()[1] < 3.5) && (cit->position()[2] > 2.5) && (cit->position()[2] < 3.5)){
                            //     std::cout << "dist3Dabs2: " << dist3Dabs << "\n";
                            //     std::cout << "weight2: " << weight << "\n";
                            // }

                            //norm += weight; // add to normalization constant
                            direction += dist3D*weight/dist3Dabs; // add to direction vector
                      }

                      // ONLY DEUGGING
                      // if (count != 2){
                      //   std::cout << "FAIL, some center missing in TDforce.cpp - probably error in communication.\n";
                      //   exit(1);
                      //   return;
                      // }
                      // ONLY DEUGGING


                      real mindist3Dabs = sqrt(mindist3D.sqr()); // calculate overall smallest absolute distance

                      // if ((direction.sqr() > 0.0) && (cit->position()[0] < 1.50216) && (cit->position()[0] > 0.80216) && (cit->position()[1] > 2.5) && (cit->position()[1] < 3.5) && (cit->position()[2] > 2.5) && (cit->position()[2] < 3.5)){
                      //   std::cout << "cit->position(): " << cit->position() << "\n";
                      //   std::cout << "mindist3Dabs: " << mindist3Dabs << "\n";
                      //   std::cout << "sqrt(direction.sqr()): " << sqrt(direction.sqr()) << "\n";
                      //   std::cout << "direction: " << direction << "\n\n";
                      // }

                      //if (mindist3Dabs>0.0) {
                         // read fforce from table for overall smallest distance (should give zero, if particle is inside a full atomistic area!)
                         //real fforce = table->getForce(mindist3Dabs);

                         // if (norm > 0.0){
                         //   fforce /= norm; // normalize with norm
                         // }

                         // substract td force
                       if ( (mindist3Dabs > startdist) && (mindist3Dabs < enddist) && (direction.sqr() > 0.0) ) {
                          real fforce = table->getForce(mindist3Dabs);
                          fforce /= sqrt(direction.sqr());
                          cit->force() += (direction * fforce); // and apply into the weighted direction
                       }
                       //}

                    }
                    else{
                      std::cout << "In TDforce: Trying to use moving AdResS region with adaptive region slab geometry. This is not implemented and doesn't make much sense anyhow.\n";
                      exit(1);
                      return;
                    }

                  }
                  else{

                    if (sphereAdr){ // spherical adress region
                       // calculate distance from reference point
                       Real3D dist3D;
                       bc.getMinimumImageVectorBox(dist3D,cit->getPos(),center); //pos - center
                       real dist = sqrt(dist3D.sqr());

                       if (dist>0.0) {
                         // read fforce from table
                         real fforce = table->getForce(dist);
                         fforce /= dist;

                         // substract td force
                         cit->force() += (dist3D * fforce);
                       }
                    } else {
                       // use this if you need 1-dir force only!
                       real d1 = cit->getPos()[0] - center[0];
                       real d1abs = fabs(d1);
                       real force = table->getForce(d1abs);
                      //std::cout<<cit->id()<<" "<<cit->force()<<" "<<cit->getPos()[0]<<" "<<center[0]<<" "<<d1<<" "<<d1abs<<" "<<force<<std::endl;
                       if (d1>0.0) {
                         cit->force()[0] += force;
                       } else {
                         cit->force()[0] -= force;
                       }
                      //std::cout<<cit->id()<<" "<<cit->force()<<std::endl;
                    }

                  }








              }
          }
    }


    //void TDforce::setCenter(real x, real y, real z){
    //        center = Real3D(x, y, z);
    //}

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void TDforce::registerPython() {

      using namespace espressopp::python;

      //void (TDforce::*pySetCenter)(real x, real y, real z)
      //                  = &TDforce::setCenter;

      void (TDforce::*pyAddForce)(int itype, const char* filename, int type)
                        = &TDforce::addForce;

      //void (TDforce::*pySetAdrRegionType)(bool _sphereAdr)
      //      = &TDforce::setAdrRegionType;


      class_<TDforce, shared_ptr<TDforce>, bases<Extension> >
        ("integrator_TDforce", init< shared_ptr<System>, shared_ptr<VerletListAdress>, real, real, int >())
        .add_property("filename", &TDforce::getFilename)
        .def("connect", &TDforce::connect)
        .def("disconnect", &TDforce::disconnect)
        //.def("setCenter", pySetCenter)
        //.def("setAdrRegionType", pySetAdrRegionType)
        .def("addForce", pyAddForce)
        ;
    }

  }
}
