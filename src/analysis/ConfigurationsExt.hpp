// ESPP_CLASS
#ifndef _ANALYSIS_ConfigurationsExt_HPP
#define _ANALYSIS_ConfigurationsExt_HPP

#include "types.hpp"
#include "SystemAccess.hpp"
#include "ConfigurationExt.hpp"

namespace espresso {
  namespace analysis {

    /** Class that stores particle positions for later analysis.

        Important: this class can also be used if the number of
        particles changes between different snapshots.
    */

    typedef std::vector<ConfigurationExtPtr> ConfigurationExtList;

    class ConfigurationsExt : public SystemAccess {

    public:

      /** Constructor, allow for unlimited snapshots. */

      ConfigurationsExt(shared_ptr<System> system) : SystemAccess (system)
      { maxConfigs = 0; }

      /** set number of maximal snapshots. */

      void setCapacity(int max);

      /** get number of maximal snapshots. */

      int getCapacity();

      /** get number of available snapshots. */

      int getSize();

      ~ConfigurationsExt() {}

      /** Gake a snapshot of all current particle positions. */

      void gather();

      ConfigurationExtPtr get(int stackpos);

      ConfigurationExtPtr back();

      ConfigurationExtList all();

      void clear() { configurationsExt.clear(); }

      static void registerPython();
    
    protected:

      static LOG4ESPP_DECL_LOGGER(logger);

    private:

      void pushConfig(ConfigurationExtPtr config);
 
      ConfigurationExtList configurationsExt;

      int maxConfigs;
    };
  }
}

#endif
