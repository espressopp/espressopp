#include "python.hpp"
#include "Temperature.hpp"

#define LOG4ESPP_LEVEL_DEBUG

namespace espresso {
  namespace analysis {
    real Temperature::compute() const {
      /* The idea is to loop over all particles (through the storage
         which probably means through the system) and compute the
         sum of v squared. Then call mpi_reduce on the controller.*/
      real T = 0.0;
      return T;
    }
  }
}
