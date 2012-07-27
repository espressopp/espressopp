// ESPP_CLASS
#ifndef _ANALYSIS_MEANSQUAREDISPL_HPP
#define _ANALYSIS_MEANSQUAREDISPL_HPP

#include "ConfigsParticleDecomp.hpp"

namespace espresso {
  namespace analysis {

    /*
     * Class derives from ConfigsParticleDecomp.
     * 
     * This implementation of mean square displacement calculation does not take into
     * account particle masses. It is correct if all the particles has equal masses only.
     * Otherwise it should be modified.
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
