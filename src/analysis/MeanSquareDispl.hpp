// ESPP_CLASS
#ifndef _ANALYSIS_MEANSQUAREDISPL_HPP
#define _ANALYSIS_MEANSQUAREDISPL_HPP

#include "ConfigsParticleDecomp.hpp"

namespace espresso {
  namespace analysis {

    /*
     * Class derives from ConfigsParticleDecomp
    */

    class MeanSquareDispl : public ConfigsParticleDecomp {

    public:
      
      MeanSquareDispl(shared_ptr<System> system): ConfigsParticleDecomp (system){
        key = "unfolded";
      }
      ~MeanSquareDispl() {}
      
      virtual python::list compute() const;

      static void registerPython();
    };
  }
}

#endif
