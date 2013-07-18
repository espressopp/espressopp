// ESPP_CLASS
#ifndef _ANALYSIS_CONFIGURATIONS_HPP
#define _ANALYSIS_CONFIGURATIONS_HPP

#include "types.hpp"
#include "SystemAccess.hpp"
#include "Configuration.hpp"

namespace espresso {
  namespace analysis {

    /** Class that stores particle positions for later analysis.

        Important: this class can also be used if the number of
        particles changes between different snapshots.
    */

    typedef std::vector<ConfigurationPtr> ConfigurationList;

    class Configurations : public SystemAccess {

    public:

      /** Constructor, allow for unlimited snapshots. */
      Configurations(shared_ptr<System> system) : SystemAccess (system) {
    	  gatherPos = true;
    	  gatherVel = false;
    	  gatherForce = false;
    	  gatherRadius = false;
    	  maxConfigs = 0;
      }
      Configurations(shared_ptr<System> system, bool _pos, bool _vel, bool _force, bool _radius, bool _folded)
                    : SystemAccess (system), gatherPos(_pos), gatherVel(_vel), gatherForce(_force), gatherRadius(_radius), folded(_folded)
      { maxConfigs = 0; }

      /** set number of maximal snapshots. */

      void setCapacity(int max);

      /** get number of maximal snapshots. */

      int getCapacity();

      /** get number of available snapshots. */

      int getSize();

      ~Configurations() {}

      /** Gake a snapshot of all current particle positions. */

      void gather();

      ConfigurationPtr get(int stackpos);

      ConfigurationPtr back();

      ConfigurationList all();

      void clear() { configurations.clear(); }

      static void registerPython();
    
    protected:

      static LOG4ESPP_DECL_LOGGER(logger);

    private:

      void pushConfig(ConfigurationPtr config);
 
      ConfigurationList configurations;

      int maxConfigs;

      bool gatherPos, gatherVel, gatherForce, gatherRadius, folded;
    };
  }
}

#endif
