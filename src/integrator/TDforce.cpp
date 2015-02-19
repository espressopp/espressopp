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


namespace espresso {
  using namespace iterator;

  namespace integrator {

    LOG4ESPP_LOGGER(TDforce::theLogger, "TDforce");

    TDforce::TDforce(shared_ptr<System> system, shared_ptr<VerletListAdress> _verletList)
    :Extension(system),verletList(_verletList) {

        type = Extension::FreeEnergyCompensation;

        if (verletList->getAdrCenterSet()) { //adress region centre is fixed in space
          //center=**(verletList->adrPositions.begin()); 
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


                  if (!(verletList->getAdrCenterSet())) { //moving adress region
//TODO if length(adrPositions) > 1 print warning (for the moment only works when there's only one particle in adrPositions)
                    center=**(verletList->adrPositions.begin());
                  }

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


    //void TDforce::setCenter(real x, real y, real z){
    //        center = Real3D(x, y, z);
    //}

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void TDforce::registerPython() {

      using namespace espresso::python;

      //void (TDforce::*pySetCenter)(real x, real y, real z)
      //                  = &TDforce::setCenter;

      void (TDforce::*pyAddForce)(int itype, const char* filename, int type)
                        = &TDforce::addForce;

      //void (TDforce::*pySetAdrRegionType)(bool _sphereAdr)
      //      = &TDforce::setAdrRegionType;


      class_<TDforce, shared_ptr<TDforce>, bases<Extension> >
        ("integrator_TDforce", init< shared_ptr<System>, shared_ptr<VerletListAdress> >())
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
