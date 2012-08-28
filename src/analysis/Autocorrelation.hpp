// ESPP_CLASS
#ifndef _ANALYSIS_AUTOCORRELATION_HPP
#define _ANALYSIS_AUTOCORRELATION_HPP

#include "types.hpp"
#include "SystemAccess.hpp"

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

    class Autocorrelation : public SystemAccess {

    public:
      // Constructor, allow for unlimited snapshots.
      Autocorrelation(shared_ptr<System> system): SystemAccess (system){
      }
      ~Autocorrelation(){
        valueList.clear();
      }

      // get number of available snapshots. Returns the size of ValueList
      unsigned int getListSize() const;

      // Take a snapshot (save the current value)
      void gather(Real3D);
      
      // Get a configuration from ConfigurationList
      Real3D getValue(unsigned int position) const;

      // it returns all values
      vector<Real3D> all() const;

      // it erases all the configurations from ConfigurationList
      void clear(){
        valueList.clear();
      }

      python::list compute();

      static void registerPython();
    
    private:

      void pushValue(Real3D);
 
      // the list of snapshots
      vector<Real3D> valueList;
    };
  }
}

#endif
