// ESPP_CLASS
#ifndef _ANALYSIS_VISCOSITY_HPP
#define _ANALYSIS_VISCOSITY_HPP

#include "types.hpp"
#include "Autocorrelation.hpp"
#include "python.hpp"

using namespace std;

namespace espresso {
  namespace analysis {

    /*
     * Class stores some single value in time for later calculation (using parallel
     * computation) of autocorrelation function. It is useful for example for viscosity
     * calculations.
     * 
     * !Important! It should be the same time period between snapshots.
    */
    
    // now the single value is Real3D
    // TODO probably template realization

    class Viscosity : public Autocorrelation {

    public:
      // Constructor, allow for unlimited snapshots.
      Viscosity(shared_ptr<System> system): Autocorrelation (system){
      }
      ~Viscosity(){
      }

      // Take a snapshot (save the current value of nonlinar component of pressure tensor)
      void gather();
      
      python::list compute(real t0, real dt, real T);

      static void registerPython();
      
    };
  }
}

#endif
