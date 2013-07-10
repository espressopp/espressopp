// ESPP_CLASS
#ifndef _ANALYSIS_VELOCITYAUTOCORRELATION_HPP
#define _ANALYSIS_VELOCITYAUTOCORRELATION_HPP

#include "ConfigsParticleDecomp.hpp"

namespace espresso {
  namespace analysis {

    /*
     * Class derives from ConfigsParticleDecomp
    */
    class VelocityAutocorrelation : public ConfigsParticleDecomp {

    public:
      
      VelocityAutocorrelation(shared_ptr<System> system): ConfigsParticleDecomp (system){
        // by default calculation progress is printed
        setPrint_progress(true);
        
        key = "velocity";
      }
      ~VelocityAutocorrelation() {}
      
      virtual python::list compute() const;

      void setPrint_progress(bool _print_progress){
        print_progress = _print_progress;
      }
      bool getPrint_progress(){return print_progress;}
      
      static void registerPython();
    private:
      bool print_progress;
    };
  }
}

#endif
