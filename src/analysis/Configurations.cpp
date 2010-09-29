#include "python.hpp"
#include "Configurations.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "mpi.h"
#include <cmath>

using namespace espresso;

#define DEFAULT_TAG 71

namespace espresso {
  namespace analysis {

    using namespace iterator;

    LOG4ESPP_LOGGER(Configurations::logger, "Configurations");

    Configurations::Configuration::Configuration(int nParticles) 
    {
      coordinates = new real[3*nParticles];
      this->nParticles = nParticles;
    }

    Configurations::Configuration::~Configuration()
    {
      delete coordinates;
    }

    void Configurations::Configuration::set(int index, real x, real y, real z)
    {
      if (index < 0 || index >= nParticles) {
         LOG4ESPP_ERROR(logger, "index = " << index << " out of range" <<
                         ", must be >= 0 and < " << nParticles);
      } else {
        coordinates[3*index]   = x;
        coordinates[3*index+1] = y;
        coordinates[3*index+2] = z;
      }
    }

    void Configurations::setCapacity(int max) 
    {
      if (max < 0) {
         LOG4ESPP_ERROR(logger, "number for maximal configurations must be positive");
      } else {
         maxConfigs = max;
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

    /** Get number of particles of a snapshop on the stack. */

    int Configurations::getNParticles(int stackpos)
    {
      int nconfigs = configurations.size();

      if (0 <= stackpos and stackpos < nconfigs) {
        return configurations[nconfigs - 1 - stackpos]->nParticles;
      }
    }

    Real3D Configurations::getCoordinates(int index, int stackpos) 
    {
      int nconfigs = configurations.size();

      if (0 <= stackpos and stackpos < nconfigs) {
        real* coords = configurations[nconfigs - 1 - stackpos]->coordinates;
        coords += 3*index;
        return Real3D(coords[0], coords[1], coords[2]);
      }
    }

    void Configurations::push() {

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

        ids[i] = cit->p.id;

        coordinates[3*i]   = cit->r.p[0];
        coordinates[3*i+1] = cit->r.p[1];
        coordinates[3*i+2] = cit->r.p[2];

        i++;
      }

      if (i != myN) {
        LOG4ESPP_ERROR(logger, "mismatch for number of local particles");
      }

      // each proc sends its data to master process

      if (system.comm->rank() == 0) {

         ConfigurationPtr config = make_shared<Configuration> (totalN);

         // root process collects data from all processors and sets it

         int nproc = system.comm->size();

         for (int iproc = 0; iproc < nproc; iproc++) {

           int nother;

           if (iproc) {
   
              MPI_Status status;
              MPI_Request request;

              int nIds, nDoubles;   // number of received values
              int tmp;

              LOG4ESPP_DEBUG(logger, "receive tags from " << iproc);

              MPI_Irecv(ids, maxN, MPI_INT, iproc, DEFAULT_TAG, *system.comm, &request);
              // send source processor an empty message that I am read to receive
              MPI_Send(&tmp, 0, MPI_INT, iproc, DEFAULT_TAG, *system.comm);
              // wait for the message
              MPI_Wait(&request, &status);
              MPI_Get_count(&status, MPI_INT, &nIds);

              MPI_Irecv(coordinates, 3*maxN, MPI_DOUBLE, iproc, DEFAULT_TAG, *system.comm, &request);
              MPI_Send(&tmp, 0, MPI_INT, iproc, DEFAULT_TAG, *system.comm);
              MPI_Wait(&request, &status);
              MPI_Get_count(&status, MPI_DOUBLE, &nDoubles);
  
              if (nDoubles != 3 * nIds) {
                LOG4ESPP_ERROR(logger, "serious error collecting data, got " << 
                              nIds << " ids, but " << nDoubles << " coordinates");
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

        configurations.push_back(config);

      } else {

       LOG4ESPP_INFO(logger, "proc " << system.comm->rank() << " sends data " 
                      << " of " << myN << " particles");

       // not master process, send data to master process

       int tmp;
       MPI_Status status;
       // wait for a signal (empty message) before sending
       MPI_Recv(&tmp, 0, MPI_INT, 0, DEFAULT_TAG, *system.comm, &status);
       MPI_Rsend(ids, myN, MPI_INT, 0, DEFAULT_TAG, *system.comm);
       MPI_Recv(&tmp, 0, MPI_INT, 0, DEFAULT_TAG, *system.comm, &status);
       MPI_Rsend(coordinates, 3*myN, MPI_DOUBLE, 0, DEFAULT_TAG, *system.comm);
       for (int i = 0; i < 3*myN; i++) printf ("coordinate[%d] = %g\n", i, coordinates[i]);
      }

      // ToDo: remove first configuration if capacity is exhausted

      // master process saves the configuration

      delete [] coordinates;
      delete [] ids;
    }

    real Configurations::compute() const {
      return -1.0;
    }

    void Configurations::registerPython() {
      using namespace espresso::python;
      class_<Configurations, bases< Observable > >
        ("analysis_Configurations", init< shared_ptr< System > >())
      .add_property("size", &Configurations::getSize)
      .add_property("capacity", &Configurations::getCapacity, &Configurations::setCapacity)
      .def("push", &Configurations::push)
      .def("getNParticles", &Configurations::getNParticles)
      .def("getCoordinates", &Configurations::getCoordinates)
      ;
    }
  }
}
