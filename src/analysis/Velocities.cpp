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
#include "Velocities.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "mpi.h"
#include <cmath>

using namespace espressopp;

#define DEFAULT_TAG 71

namespace espressopp {
  namespace analysis {

    using namespace iterator;

    LOG4ESPP_LOGGER(Velocities::logger, "Velocities");

    void Velocities::setCapacity(int max) 
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

    int Velocities::getCapacity()
    {
      return maxConfigs;
    }

    int Velocities::getSize()
    {
      return configurations.size();
    }

   
    ConfigurationList Velocities::all()
    {
      return configurations;
    }

    /** Get a configuration from stack */

    ConfigurationPtr Velocities::get(int stackpos)
    {
      int nconfigs = configurations.size();
      if (0 <= stackpos and stackpos < nconfigs) {
        return configurations[nconfigs - 1 - stackpos];
      } else {
        LOG4ESPP_ERROR(logger, "Velocities::get <out-of-range>");
        return shared_ptr<Configuration>();
      }
    }

    ConfigurationPtr Velocities::back()
    {
      return configurations.back();
    }

    void Velocities::pushConfig(ConfigurationPtr config)
    {
      int nconfs = configurations.size();

      if (maxConfigs && nconfs >= maxConfigs) {

        LOG4ESPP_DEBUG(logger, "delete first configuration");

        // remove the first configuration

        configurations.erase(configurations.begin());
      }

      configurations.push_back(config);
    }

    void Velocities::gather() {

      System& system = getSystemRef();
  
      // determine number of local particles and total particles

      int myN = system.storage->getNRealParticles();
      int maxN;   // maximal number of particles one processor has
      int totalN; // totlal number of particles all processors have

      boost::mpi::all_reduce(*system.comm, myN, maxN, boost::mpi::maximum<int>());
      boost::mpi::all_reduce(*system.comm, myN, totalN, std::plus<int>());

      LOG4ESPP_INFO(logger, "#Partices: me = " << myN << ", max = " << maxN
                            << ", totalN = " << totalN);

      real* coordinates = new real [3 * maxN];  // buffer for gather
      int*  ids         = new int [maxN];  // buffer for gather

      // fill the buffer with my values

      CellList realCells = system.storage->getRealCells();

      int i = 0; 

      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {

        ids[i] = cit->id();
        
        // replaced position with velocity
        Real3D& pos = cit->velocity();

        coordinates[3*i]   = pos[0];
        coordinates[3*i+1] = pos[1];
        coordinates[3*i+2] = pos[2];

        i++;
      }

      if (i != myN) {
        LOG4ESPP_ERROR(logger, "mismatch for number of local particles");
      }

      // each proc sends its data to master process

      if (system.comm->rank() == 0) {

         ConfigurationPtr config = make_shared<Configuration> (); //totalN

         // root process collects data from all processors and sets it

         int nproc = system.comm->size();

         for (int iproc = 0; iproc < nproc; iproc++) {

           int nother;

           if (iproc) {
   
              int nIds, nCoords;   // number of received values
              int tmp;

              boost::mpi::request req;
              boost::mpi::status  stat;

              LOG4ESPP_DEBUG(logger, "receive tags from " << iproc);

              req = system.comm->irecv<int>(iproc, DEFAULT_TAG, ids, maxN);
              system.comm->send(iproc, DEFAULT_TAG, 0);
              stat = req.wait();
              nIds = *stat.count<int>();

              req = system.comm->irecv<real>(iproc, DEFAULT_TAG, coordinates, 3*maxN);
              system.comm->send(iproc, DEFAULT_TAG, 0);
              stat = req.wait();
              nCoords = *stat.count<real>();
  
              // make sure to have 3 coordinate values for each id

              if (nCoords != 3 * nIds) {
                LOG4ESPP_ERROR(logger, "serious error collecting data, got " << 
                              nIds << " ids, but " << nCoords << " coordinates");
              }

              nother = nIds;

           } else {
             nother = myN;
           }
   
           LOG4ESPP_INFO(logger, "add " << nother << " coordinates of proc " << iproc);

           for (int i = 0; i < nother; i++) {

             int index = ids[i];

             LOG4ESPP_INFO(logger, "set coordianates of particle with id = " << index <<
                                    ": " << coordinates[3*i] << " " <<  coordinates[3*i+1] << " " << coordinates[3*i+2]);

             config->set(index, coordinates[3*i], coordinates[3*i+1], coordinates[3*i+2]);
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
       system.comm->irecv<int>(0, DEFAULT_TAG, tmp);
       system.comm->send<real>(0, DEFAULT_TAG, coordinates, 3*myN);
      }

      // ToDo: remove first configuration if capacity is exhausted

      // master process saves the configuration

      delete [] coordinates;
      delete [] ids;
    }

    // Python wrapping

    void Velocities::registerPython() {

      using namespace espressopp::python;

      //class_<ConfigurationList> ("_ConfigurationList", no_init)
      //.def("__iter__", boost::python::iterator<ConfigurationList>())
      //;

      class_<Velocities>
        ("analysis_Velocities", init< shared_ptr< System > >())
      .add_property("size", &Velocities::getSize)
      .add_property("capacity", &Velocities::getCapacity, 
                                &Velocities::setCapacity)
      .def("gather", &Velocities::gather)
      .def("__getitem__", &Velocities::get)
      .def("back", &Velocities::back)
      .def("all", &Velocities::all)
      .def("clear", &Velocities::clear)
      ;
    }
  }
}
