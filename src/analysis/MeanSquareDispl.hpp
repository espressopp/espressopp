// ESPP_CLASS
#ifndef _ANALYSIS_MEANSQUAREDISPL_HPP
#define _ANALYSIS_MEANSQUAREDISPL_HPP

#include "ConfigsParticleDecomp.hpp"

namespace espresso {
  namespace analysis {

    /*
     * Class derived from ConfigsParticleDecomp.
     * 
     * This implementation of mean square displacement calculation does not take into
     * account particle masses. It is correct if all the particles have equal masses only.
     * Otherwise it should be modified.
    */

    class MeanSquareDispl : public ConfigsParticleDecomp {

    public:
      
      MeanSquareDispl(shared_ptr<System> system): ConfigsParticleDecomp (system){
        // by default 
        setPrint_progress(true);
        
        key = "unfolded";
      }
      ~MeanSquareDispl() {}
      
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
