// ESPP_CLASS
#ifndef _IO_DUMPXYZ_HPP
#define _IO_DUMPXYZ_HPP

#include "ParticleAccess.hpp"
#include "integrator/MDIntegrator.hpp"

/////#include "../System.hpp"
//#include "analysis/ConfigurationsExt.hpp"

//using namespace espresso::analysis;

namespace espresso {
  namespace io{

    class DumpXYZ : public ParticleAccess {

    public:

      DumpXYZ(shared_ptr<System> system, 
              shared_ptr<integrator::MDIntegrator> _integrator,
              std::string _file_name,
              bool _unfolded) : 
                        ParticleAccess(system), 
                        integrator(_integrator),
                        file_name( _file_name ),
                        unfolded(_unfolded){ 
                          
        //ConfigurationsExt _conf( system );
        //conf = _conf;
      }
      ~DumpXYZ() {}

      void perform_action(){
        dump();
      }
      
      void dump();
      
      std::string getFilename(){return file_name;}
      void setFilename(std::string v){file_name = v;}
      bool getUnfolded(){return unfolded;}
      void setUnfolded(bool v){unfolded = v;}

      static void registerPython();
    
    protected:

      //static LOG4ESPP_DECL_LOGGER(logger);

    private:
      
      //ConfigurationsExt conf( System() );
      // integrator we need to know an integration step
      shared_ptr<integrator::MDIntegrator> integrator;
      
      std::string file_name;
      
      bool unfolded;  // one can choose folded or unfolded coordinates, by default it is folded
    };
  }
}

#endif
