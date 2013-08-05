#include <boost/python.hpp>
#include "mpi.h"
#include "python.hpp"
#include "DumpXYZ.hpp"

//#include "ParticleAccess.hpp"
//#include "SystemAccess.hpp"

using namespace espresso;

namespace espresso {
  namespace io {
      
    /*
    void dump(){
      
      //System& system1 = getSystemRef();
      /*
      System& system = getSystemRef();
      //int step = system.integrator->getStep();
      int myN = system.storage->getNRealParticles();
      int myrank = system.comm->rank();
      
      std::cout << "N:  " << myN << "  cpu: "<< myrank << std::endl;
      
      std::cout << " Check "<< std::endl;
    }
    */
      
    // Python wrapping
    void DumpXYZ::registerPython() {

      using namespace espresso::python;

      class_<DumpXYZ, bases<ParticleAccess>, boost::noncopyable >
      ("io_DumpXYZ", init< shared_ptr< System > >())
        //.def("dump", &DumpXYZ::dump)
      ;
    }
  }
}
