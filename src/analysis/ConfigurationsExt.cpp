/*
  Copyright (C) 2016
      Max Planck Institute for Polymer Research & JGU Mainz
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
#include "ConfigurationsExt.hpp"
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

    LOG4ESPP_LOGGER(ConfigurationsExt::logger, "ConfigurationsExt");

    void ConfigurationsExt::setCapacity(int max)
    {
      if (max < 0) {
         LOG4ESPP_ERROR(logger, "number for maximal configurations must be positive");
         return;
      }

      maxConfigs = max;

      // if capacity has been reduced, delete configurations

      int nconfigs = configurationsExt.size();

      if (maxConfigs < nconfigs) {

         int diff = nconfigs - maxConfigs;

         LOG4ESPP_INFO(logger, "delete " << diff <<
              " configurations due to restricted capacity");

         configurationsExt.erase(configurationsExt.begin(),
                              configurationsExt.begin() + diff);
      }

      unfolded = true; // by default it writes unfolded coordinates
    }

    int ConfigurationsExt::getCapacity()
    {
      return maxConfigs;
    }

    int ConfigurationsExt::getSize()
    {
      return configurationsExt.size();
    }


    ConfigurationExtList ConfigurationsExt::all()
    {
      return configurationsExt;
    }

    /** Get a configuration from stack */

    ConfigurationExtPtr ConfigurationsExt::get(int stackpos)
    {
      int nconfigs = configurationsExt.size();
      if (0 <= stackpos and stackpos < nconfigs) {
        return configurationsExt[nconfigs - 1 - stackpos];
      } else {
        LOG4ESPP_ERROR(logger, "ConfigurationsExt::get <out-of-range>");
        return shared_ptr<ConfigurationExt>();
      }
    }

    ConfigurationExtPtr ConfigurationsExt::back()
    {
      return configurationsExt.back();
    }

    void ConfigurationsExt::pushConfig(ConfigurationExtPtr config)
    {
      int nconfs = configurationsExt.size();

      if (maxConfigs && nconfs >= maxConfigs) {

        LOG4ESPP_DEBUG(logger, "delete first configuration");

        // remove the first configuration

        configurationsExt.erase(configurationsExt.begin());
      }

      configurationsExt.push_back(config);
    }

    void ConfigurationsExt::gather() {

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
      real* velocities  = new real [3 * maxN];  // buffer for gather
      int*  ids         = new int [maxN];  // buffer for gather

      // fill the buffer with my values

      CellList realCells = system.storage->getRealCells();

      int i = 0;

      if( unfolded ){
        for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {

          ids[i] = cit->id();

          Real3D& pos = cit->position();
          Real3D& vel = cit->velocity();
          Int3D& img = cit->image();
          Real3D L = system.bc->getBoxL();

          coordinates[3*i]   = pos[0] + img[0] * L[0];
          coordinates[3*i+1] = pos[1] + img[1] * L[1];
          coordinates[3*i+2] = pos[2] + img[2] * L[2];

          velocities[3*i]   = vel[0];
          velocities[3*i+1] = vel[1];
          velocities[3*i+2] = vel[2];

          i++;
        }
      }
      else{
        for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {

          ids[i] = cit->id();

          Real3D& pos = cit->position();
          Real3D& vel = cit->velocity();

          coordinates[3*i]   = pos[0];
          coordinates[3*i+1] = pos[1];
          coordinates[3*i+2] = pos[2];

          velocities[3*i]   = vel[0];
          velocities[3*i+1] = vel[1];
          velocities[3*i+2] = vel[2];

          i++;
        }
      }

      if (i != myN) {
        LOG4ESPP_ERROR(logger, "mismatch for number of local particles");
      }

      // each proc sends its data to master process

      if (system.comm->rank() == 0) {

         ConfigurationExtPtr config = make_shared<ConfigurationExt> (); //totalN

         // root process collects data from all processors and sets it

         int nproc = system.comm->size();

         for (int iproc = 0; iproc < nproc; iproc++) {

           int nother;

           if (iproc) {

              int nIds, nCoords, nVelocs; // number of received values
              //int tmp;

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


              req = system.comm->irecv<real>(iproc, DEFAULT_TAG, velocities, 3*maxN);
              system.comm->send(iproc, DEFAULT_TAG, 0);
              stat = req.wait();
              nVelocs = *stat.count<real>();

              // make sure to have 3 velocities values for each id

              if (nVelocs != 3 * nIds) {
                LOG4ESPP_ERROR(logger, "serious error collecting data, got " <<
                              nIds << " ids, but " << nVelocs << " velocities");
              }

              nother = nIds;

           } else {
             nother = myN;
           }


           LOG4ESPP_INFO(logger, "add " << nother << " coordinates of proc " << iproc);

           for (int i = 0; i < nother; i++) {

             int index = ids[i];

             LOG4ESPP_INFO(logger, "set coordianates of particle with id = " << index <<
                                   ": " << coordinates[3*i] << " " <<  coordinates[3*i+1] << " " << coordinates[3*i+2] <<
                                   "and velocities: " <<  velocities[3*i] << " " <<  velocities[3*i+1] << " " << velocities[3*i+2]);

             RealND _vec(6);  // p[0].p[1].p[2].v[0].v[1].v[2]

             for (int k=0; k<3; k++)
               _vec.setItem(k, coordinates[3*i + k]);

             for (int k=0; k<3; k++)
               _vec.setItem(k+3, velocities[3*i + k]);

             config->set(index, _vec);
           }
        }

        LOG4ESPP_INFO(logger, "save the latest configuration");

        pushConfig(config);

      } else {

       LOG4ESPP_INFO(logger, "proc " << system.comm->rank() << " sends data "
                      << " of " << myN << " particles");

       // not master process, send data to master process

       int tmp;

       boost::mpi::request requests[2];
       boost::mpi::status stat;

       // wait for a signal (empty message) before sending

       requests[0] = system.comm->irecv<int>(0, DEFAULT_TAG, tmp);
       system.comm->send<int>(0, DEFAULT_TAG, ids, myN);
       stat = requests[0].wait();
       requests[1] = system.comm->irecv<int>(0, DEFAULT_TAG, tmp);
       stat = requests[1].wait();
       system.comm->send<real>(0, DEFAULT_TAG, coordinates, 3*myN);
       system.comm->send<real>(0, DEFAULT_TAG, velocities, 3*myN);
      }

      // ToDo: remove first configuration if capacity is exhausted

      // master process saves the configuration

      delete [] coordinates;
      delete [] velocities;
      delete [] ids;
    }

    // Python wrapping

    void ConfigurationsExt::registerPython() {

      using namespace espressopp::python;

      class_<ConfigurationExtList> ("_ConfigurationExtList", no_init)
      .def("__iter__", boost::python::iterator<ConfigurationExtList>())
      ;

      class_<ConfigurationsExt>
        ("analysis_ConfigurationsExt", init< shared_ptr< System > >())
      .add_property("size", &ConfigurationsExt::getSize)
      .add_property("capacity", &ConfigurationsExt::getCapacity,
                                &ConfigurationsExt::setCapacity)
      .add_property("unfolded", &ConfigurationsExt::getUnfolded,
                                &ConfigurationsExt::setUnfolded)
      .def("gather", &ConfigurationsExt::gather)
      .def("__getitem__", &ConfigurationsExt::get)
      .def("back", &ConfigurationsExt::back)
      .def("all", &ConfigurationsExt::all)
      .def("clear", &ConfigurationsExt::clear)
      ;
    }
  }
}
