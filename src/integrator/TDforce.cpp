/*
  Copyright (C) 2017,2018
      Jakub Krajniak (jkrajniak at gmail.com), Max Planck Institute for Polymer Research
  Copyright (C) 2012,2013,2014,2015,2016
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

    TDforce::TDforce(shared_ptr<System> system, shared_ptr<VerletListAdress> _verletList, real _startdist, real _enddist, int _edgeweightmultiplier, bool _slow)
    :Extension(system), verletList(_verletList), startdist(_startdist), enddist(_enddist), edgeweightmultiplier(_edgeweightmultiplier), slow(_slow) {

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
      if(slow){
        _applyForce = integrator->aftCalcSlow.connect(
            boost::bind(&TDforce::applyForce, this), boost::signals2::at_front);
      }
      else{
        _applyForce = integrator->aftCalcF.connect(
            boost::bind(&TDforce::applyForce, this));
      }
    }

    void TDforce::disconnect(){
        _applyForce.disconnect();
    }

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

        if (forces.count(type) == 0) {
          forces.insert(std::make_pair(type,table));
        } else {
          forces[type].reset();  // decrease use_counter by one
          forces[type] = table;
        }
    }


    real TDforce::getForce(longint type_id, real dist) {
      std::unordered_map<int, Table>::iterator tableIt;
      tableIt = forces.find(type_id);
      Table table;
      if (tableIt != forces.end())
        table = tableIt->second;
      if (table)
        return table->getForce(dist);
      else
        throw std::runtime_error("Table not found");
    }


    void TDforce::applyForce() {
          LOG4ESPP_DEBUG(theLogger, "apply TD force");

          System& system = getSystemRef();
          const bc::BC& bc = *getSystemRef().bc;  // boundary conditions

          std::unordered_map<int, Table>::iterator tableIt;
          Table table;
          // iterate over CG particles
          CellList cells = system.storage->getRealCells();
          for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
              tableIt = forces.find(cit->getType());
              if (tableIt != forces.end()) { //because there may be CG particles to which TD force is not applied
                table = tableIt->second;
              }

              if (table) {

                  if (!(verletList->getAdrCenterSet())) { // adress regions defined based on particles

                    if (sphereAdr){

                      if(verletList->getAdrList().size() > 1){

                        if ((startdist == 0.0) && (enddist == 0.0)){
                          std::cout << "In TDforce: Trying to apply TD force with moving, particle-based spherical AdResS region. However, both startdist and enddist set to 0. This does not make sense.\n";
                          exit(1);
                          return;
                        }

                        real buffer = 0.000001;
                        real width = fabs(enddist - startdist + 2.0*buffer);  // width of region where TD force acts
                        std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();   // get positions
                        Real3D pa = **it2;
                        Real3D dist3D;
                        bc.getMinimumImageVectorBox(dist3D,cit->getPos(), pa); // calculate vector between particle and first center

                        Real3D mindist3D = dist3D; // set mindist3D before loop

                        real dist3Dabs = sqrt(dist3D.sqr()); // calculate absolute distance

                        real weight = 0.0;  // initialize weight

                        // weighting scheme: only particles that are in a hybrid region contribute
                        if ((dist3Dabs < enddist + buffer) && (dist3Dabs > startdist - buffer)){
                          weight = 1.0 - (dist3Dabs - startdist - buffer) / width; // 0 at edge to CG region, 1 at edge to atomistic region

                          // Direction modification at the edges
                          weight = pow(weight, edgeweightmultiplier);
                        }
                        else{
                          weight = 0.0; // particle outside of hybrid region
                        }

                        Real3D direction = dist3D*weight/dist3Dabs; // normalized but weighted direction vector

                        // Loop over all other centers
                        ++it2;
                        for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                              pa = **it2;
                              verletList->getSystem()->bc->getMinimumImageVector(dist3D, cit->getPos(), pa); // calculate vector between particle and other centers
                              if (dist3D.sqr() < mindist3D.sqr()) mindist3D = dist3D; // make it the minimum if shortest length so far
                              dist3Dabs = sqrt(dist3D.sqr()); // calculate absolute distance

                              // calculate weights again, as above
                              if ((dist3Dabs < enddist + buffer) && (dist3Dabs > startdist - buffer)){
                                weight = 1.0 - (dist3Dabs - startdist - buffer) / width;

                                // Direction modification at the edges
                                weight = pow(weight, edgeweightmultiplier);
                              }
                              else{
                                weight = 0.0;
                              }

                              //norm += weight; // add to normalization constant
                              direction += dist3D*weight/dist3Dabs; // add to direction vector
                        }

                        real mindist3Dabs = sqrt(mindist3D.sqr()); // calculate overall smallest absolute distance

                        // apply td force
                         if ( (mindist3Dabs > startdist) && (mindist3Dabs < enddist) && (direction.sqr() > 0.0) ) {
                            real fforce = table->getForce(mindist3Dabs);
                            fforce /= sqrt(direction.sqr());
                            cit->force() += (direction * fforce);
                         }

                      }
                      else{
                         // calculate distance from reference particle
                         center=**(verletList->adrPositions.begin());
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
                      }

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
                       if (d1>0.0) {
                         cit->force()[0] += force;
                       } else {
                         cit->force()[0] -= force;
                       }
                    }

                  }
              }
          }
    }

    real TDforce::computeTDEnergy() {

          real TDEnergy = 0.0;
          real TDEnergySum = 0.0;
          System& system = getSystemRef();
          const bc::BC& bc = *getSystemRef().bc;

          CellList cells = system.storage->getRealCells();
          for(CellListIterator cit(cells); !cit.isDone(); ++cit) {

              Table table = forces.find(cit->getType())->second;
              if (table) {

                  if (!(verletList->getAdrCenterSet())) {

                    if (sphereAdr){

                      if(verletList->getAdrList().size() > 1){

                        std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();
                        Real3D pa = **it2;
                        Real3D dist3D;
                        bc.getMinimumImageVectorBox(dist3D,cit->getPos(), pa);
                        Real3D mindist3D = dist3D;

                        // Loop over all other centers
                        ++it2;
                        for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                              pa = **it2;
                              verletList->getSystem()->bc->getMinimumImageVector(dist3D, cit->getPos(), pa);
                              if (dist3D.sqr() < mindist3D.sqr()) mindist3D = dist3D;
                        }

                        real mindist3Dabs = sqrt(mindist3D.sqr());
                         if ( mindist3Dabs >0.0 ) {
                            TDEnergy += table->getEnergy(mindist3Dabs);
                         }

                      }
                      else{
                         center=**(verletList->adrPositions.begin());
                         Real3D dist3D;
                         bc.getMinimumImageVectorBox(dist3D,cit->getPos(),center);
                         real dist = sqrt(dist3D.sqr());

                         if (dist>0.0) {
                           TDEnergy += table->getEnergy(dist);
                         }
                      }

                    }
                    else{
                      std::cout << "In TDforce: Trying to use moving AdResS region with adaptive region slab geometry. This is not implemented and doesn't make much sense anyhow.\n";
                      exit(1);
                      return 0.0;
                    }

                  }
                  else{

                    if (sphereAdr){
                       Real3D dist3D;
                       bc.getMinimumImageVectorBox(dist3D,cit->getPos(),center);
                       real dist = sqrt(dist3D.sqr());

                       if (dist>0.0) {
                         TDEnergy += table->getEnergy(dist);
                       }
                    } else {
                       real d1 = cit->getPos()[0] - center[0];
                       real d1abs = fabs(d1);
                       TDEnergy += table->getEnergy(d1abs);
                    }

                  }
              }

          }
          mpi::all_reduce(*getSystem()->comm, TDEnergy, TDEnergySum, std::plus<real>());
          return TDEnergySum;
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void TDforce::registerPython() {

      using namespace espressopp::python;

      void (TDforce::*pyAddForce)(int itype, const char* filename, int type)
                        = &TDforce::addForce;

      class_<TDforce, shared_ptr<TDforce>, bases<Extension> >
        ("integrator_TDforce", init< shared_ptr<System>, shared_ptr<VerletListAdress>, real, real, int, bool >())
        .add_property("filename", &TDforce::getFilename)
        .def("connect", &TDforce::connect)
        .def("disconnect", &TDforce::disconnect)
        .def("addForce", pyAddForce)
        .def("computeTDEnergy", &TDforce::computeTDEnergy)
        .def("getForce", &TDforce::getForce)
        ;
    }

  }
}
