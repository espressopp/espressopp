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
#include <boost/python.hpp>
#include "Configurations.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "bc/BC.hpp"
#include "mpi.h"
#include <cmath>

using namespace espressopp;

#define DEFAULT_TAG 71

namespace espressopp {
  namespace analysis {

    using namespace iterator;

    LOG4ESPP_LOGGER(Configurations::logger, "Configurations");

    void Configurations::setCapacity(int max) 
    {
      if (max < 0) {
         LOG4ESPP_ERROR(logger, "number for maximal configurations must be positive");
         return;
      }

      maxConfigs = max;

      // if capacity has been reduced, delete configurations

      int nconfigs = configurations.size();

      if (maxConfigs < nconfigs) {

         int diff = nconfigs - maxConfigs;

         LOG4ESPP_INFO(logger, "delete " << diff << 
              " configurations due to restricted capacity");

         configurations.erase(configurations.begin(),
                              configurations.begin() + diff);
      }

    }

    int Configurations::getCapacity()
    {
      return maxConfigs;
    }

    int Configurations::getSize()
    {
      return configurations.size();
    }

   
    ConfigurationList Configurations::all()
    {
      return configurations;
    }

    /** Get a configuration from stack */

    ConfigurationPtr Configurations::get(int stackpos)
    {
      int nconfigs = configurations.size();
      if (0 <= stackpos and stackpos < nconfigs) {
        return configurations[nconfigs - 1 - stackpos];
      } else {
        LOG4ESPP_ERROR(logger, "Configurations::get <out-of-range>");
        return shared_ptr<Configuration>();
      }
    }

    ConfigurationPtr Configurations::back()
    {
      return configurations.back();
    }

    void Configurations::pushConfig(ConfigurationPtr config)
    {
      int nconfs = configurations.size();

      if (maxConfigs && nconfs >= maxConfigs) {

        LOG4ESPP_DEBUG(logger, "delete first configuration");

        // remove the first configuration

        configurations.erase(configurations.begin());
      }

      configurations.push_back(config);
    }

    void Configurations::gather() {

      System& system = getSystemRef();
  
      // determine number of local particles and total particles

      int myN = system.storage->getNRealParticles();
      int maxN;   // maximal number of particles one processor has
      int totalN; // totlal number of particles all processors have

      boost::mpi::all_reduce(*system.comm, myN, maxN, boost::mpi::maximum<int>());
      boost::mpi::all_reduce(*system.comm, myN, totalN, std::plus<int>());

      LOG4ESPP_INFO(logger, "#Partices: me = " << myN << ", max = " << maxN
                            << ", totalN = " << totalN);

      int*  ids         = new int [maxN];  // buffer for gather
      Real3D* coordinates;
      Real3D* velocities;
      Real3D* forces;
      real* radii;

      if (gatherPos)    coordinates = new Real3D [maxN];  // buffer for gather
      if (gatherVel)    velocities  = new Real3D [maxN];  // buffer for gather
      if (gatherForce)  forces      = new Real3D [maxN];  // buffer for gather
      if (gatherRadius) radii       = new real [maxN];  // buffer for gather

      // fill the buffer with my values

      CellList realCells = system.storage->getRealCells();

      int i = 0; 
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        ids[i] = cit->id();
        if (gatherPos) {
          Real3D pos = cit->position();
          Int3D img = cit->image();
          if (folded)
        	system.bc->foldPosition(pos, img);
          else
        	system.bc->unfoldPosition(pos, img);
          coordinates[i] = pos;
        }
        if (gatherVel)    velocities[i]  = cit->velocity();
        if (gatherForce)  forces[i]      = cit->force();
        if (gatherRadius) radii[i]       = cit->radius();
        i++;
      }

      if (i != myN) {
        LOG4ESPP_ERROR(logger, "mismatch for number of local particles");
      }

      // each proc sends its data to master process

      if (system.comm->rank() == 0) {

         ConfigurationPtr config = make_shared<Configuration> (gatherPos, gatherVel, gatherForce, gatherRadius); //totalN

         // root process collects data from all processors and sets it

         int nproc = system.comm->size();

         for (int iproc = 0; iproc < nproc; iproc++) {
           int nother;
           if (iproc) {
              int nIds, nVals;   // number of received values
              int tmp;
              boost::mpi::request req;
              boost::mpi::status  stat;
              LOG4ESPP_DEBUG(logger, "receive tags from " << iproc);
              req = system.comm->irecv<int>(iproc, DEFAULT_TAG, ids, maxN);
              system.comm->send(iproc, DEFAULT_TAG, 0);
              stat = req.wait();
              nIds = *stat.count<int>();
              if (gatherPos) {
                req = system.comm->irecv<Real3D>(iproc, DEFAULT_TAG, coordinates, maxN);
                system.comm->send(iproc, DEFAULT_TAG, 0);
                stat = req.wait();
                nVals = *stat.count<Real3D>();
                if (nVals != nIds) {
                  LOG4ESPP_ERROR(logger, "serious error collecting data, got " <<
                                nIds << " ids, but " << nVals << " coordinates");
                }
              }
              if (gatherVel) {
                req = system.comm->irecv<Real3D>(iproc, DEFAULT_TAG, velocities, maxN);
                system.comm->send(iproc, DEFAULT_TAG, 0);
                stat = req.wait();
                nVals = *stat.count<Real3D>();
                if (nVals != nIds) {
                  LOG4ESPP_ERROR(logger, "serious error collecting data, got " <<
                                nIds << " ids, but " << nVals << " velocities");
                }
              }
              if (gatherForce) {
                req = system.comm->irecv<Real3D>(iproc, DEFAULT_TAG, forces, maxN);
                system.comm->send(iproc, DEFAULT_TAG, 0);
                stat = req.wait();
                nVals = *stat.count<Real3D>();
                if (nVals != nIds) {
                  LOG4ESPP_ERROR(logger, "serious error collecting data, got " <<
                                nIds << " ids, but " << nVals << " forces");
                }
              }
              if (gatherRadius) {
                req = system.comm->irecv<real>(iproc, DEFAULT_TAG, radii, maxN);
                system.comm->send(iproc, DEFAULT_TAG, 0);
                stat = req.wait();
                nVals = *stat.count<real>();
                if (nVals != nIds) {
                  LOG4ESPP_ERROR(logger, "serious error collecting data, got " <<
                                nIds << " ids, but " << nVals << " radii");
                }
              }
              nother = nIds;
           } else {
             nother = myN;
           }
   
           LOG4ESPP_INFO(logger, "add " << nother << " coordinates of proc " << iproc);

           for (int i = 0; i < nother; i++) {
             //LOG4ESPP_DEBUG(logger, "set coordinates of particle with id = " << index <<
             //                       ": " << coordinates[3*i] << " " <<  coordinates[3*i+1] << " " << coordinates[3*i+2]);
             int index = ids[i];
             if (gatherPos)    config->setCoordinates(index, coordinates[i]);
             if (gatherVel)    config->setVelocities(index, velocities[i]);
             if (gatherForce)  config->setForces(index, forces[i]);
             if (gatherRadius) config->setRadius(index, radii[i]);
           }
        }

        LOG4ESPP_INFO(logger, "save the latest configuration");

        pushConfig(config);

      } else {

       LOG4ESPP_INFO(logger, "proc " << system.comm->rank() << " sends data " 
                      << " of " << myN << " particles");

       // not master process, send data to master process

       int tmp;

       boost::mpi::status stat;

       // wait for a signal (empty message) before sending

       system.comm->irecv<int>(0, DEFAULT_TAG, tmp);
       system.comm->send<int>(0, DEFAULT_TAG, ids, myN);
       if (gatherPos) {
         system.comm->irecv<int>(0, DEFAULT_TAG, tmp);
         system.comm->send<Real3D>(0, DEFAULT_TAG, coordinates, myN);
       }
       if (gatherVel) {
         system.comm->irecv<int>(0, DEFAULT_TAG, tmp);
         system.comm->send<Real3D>(0, DEFAULT_TAG, velocities, myN);
       }
       if (gatherForce) {
         system.comm->irecv<int>(0, DEFAULT_TAG, tmp);
         system.comm->send<Real3D>(0, DEFAULT_TAG, forces, myN);
       }
       if (gatherRadius) {
         system.comm->irecv<int>(0, DEFAULT_TAG, tmp);
         system.comm->send<real>(0, DEFAULT_TAG, radii, myN);
       }
      }

      // ToDo: remove first configuration if capacity is exhausted

      // master process saves the configuration

      if (gatherRadius) delete [] radii;
      if (gatherForce)  delete [] forces;
      if (gatherVel)    delete [] velocities;
      if (gatherPos)    delete [] coordinates;
      delete [] ids;
    }

    // Python wrapping

    void Configurations::registerPython() {

      using namespace espressopp::python;

      class_<ConfigurationList> ("_ConfigurationList", no_init)
      .def("__iter__", boost::python::iterator<ConfigurationList>())
      ;

      class_<Configurations>
        ("analysis_Configurations", init< shared_ptr< System > >())
      .def(init<shared_ptr< System >,bool,bool,bool,bool,bool>())
      .add_property("size", &Configurations::getSize)
      .add_property("capacity", &Configurations::getCapacity, 
                                &Configurations::setCapacity)
      .def("gather", &Configurations::gather)
      .def("__getitem__", &Configurations::get)
      .def("back", &Configurations::back)
      .def("all", &Configurations::all)
      .def("clear", &Configurations::clear)
      ;
    }
  }
}
