// ESPP_CLASS
#ifndef _ANALYSIS_VELOCITIES_HPP
#define _ANALYSIS_VELOCITIES_HPP

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

    class Velocities : public SystemAccess {

    public:

      /** Constructor, allow for unlimited snapshots. */
      Velocities(shared_ptr<System> system) : SystemAccess (system)
      { maxConfigs = 0; }

      /** set number of maximal snapshots. */
      void setCapacity(int max);

      /** get number of maximal snapshots. */
      int getCapacity();

      /** get number of available snapshots. */
      int getSize();

      ~Velocities() {}

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
    };
  }
}

#endif
