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
      }
      ~VelocityAutocorrelation() {}
      
      virtual python::list compute() const;

      static void registerPython();
    };
  }
}

#endif
