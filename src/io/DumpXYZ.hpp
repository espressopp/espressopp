// ESPP_CLASS
#ifndef _IO_DUMPXYZ_HPP
#define _IO_DUMPXYZ_HPP


#include "mpi.h"
#include "types.hpp"
#include "python.hpp"
#include "ParticleAccess.hpp"
#include "SystemAccess.hpp"

#include "storage/Storage.hpp"

// ******************************************************************************************

namespace espresso {
  namespace io{

    //typedef std::vector<ConfigurationExtPtr> ConfigurationExtList;

    class DumpXYZ : public ParticleAccess {

    public:

      /** Constructor, allow for unlimited snapshots. */

      DumpXYZ(shared_ptr<System> system) : ParticleAccess(system){ 
        //maxConfigs = 0; 
        
        //System& system1 = getSystemRef();
        //int step = system.integrator->getStep();
        //int myN = system1.storage->getNRealParticles();
        //int myrank = system1.comm->rank();

        //std::cout << "CONSTRUCTOR   N:  " << myN << "  cpu: "<< myrank << std::endl;

      }
      ~DumpXYZ() {}

      void perform_action(){
        //dump();
      }
      
      //void dump();
      
      /*

      bool getUnfolded(){return unfolded;}
      void setUnfolded(bool v){unfolded = v;}

      void gather();
      ConfigurationExtPtr get(int stackpos);
      ConfigurationExtPtr back();
      ConfigurationExtList all();
      void clear() { configurationsExt.clear(); }
      */

      static void registerPython();
    
    protected:

      //static LOG4ESPP_DECL_LOGGER(logger);

    private:
/*
      void pushConfig(ConfigurationExtPtr config);
      ConfigurationExtList configurationsExt;
      int maxConfigs;
      bool unfolded;  // one can choose folded or unfolded coordinates, by default it is unfolded
 */
    };
  }
}

#endif
