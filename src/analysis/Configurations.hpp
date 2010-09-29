// ESPP_CLASS
#ifndef _ANALYSIS_CONFIGURATION_HPP
#define _ANALYSIS_CONFIGURATION_HPP

#include "types.hpp"
#include "Observable.hpp"

namespace espresso {
  namespace analysis {

    /** Class that stores particle positions for later analysis.

        Important: this class can also be used if the number of
        particles changes between different snapshots.
    */

    class Configurations : public Observable {

    public:

      /** Constructor, allow for unlimited snapshots. */

      Configurations(shared_ptr<System> system) : Observable(system) { maxConfigs = 0; }

      /** set number of maximal snapshots. */

      void setCapacity(int max);

      /** get number of maximal snapshots. */

      int getCapacity();

      /** get number of available snapshots. */

      int getSize();

      ~Configurations() {}

      /** Gake a snapshot of all current particle positions. */

      void push();

      /** Get number of particles of a snapshop on the stack. */

      int getNParticles(int stackpos);

      Real3D getCoordinates(int index, int stackpos);

      static void registerPython();
    
      real compute() const;

    protected:

      static LOG4ESPP_DECL_LOGGER(logger);

    private:

      struct Configuration {
        Configuration(int nParticles);
        ~Configuration();
        void set(int index, real x, real y, real z);
        int nParticles;     // number of particles of the configuration
        real* coordinates;   // size will be 3 * nParticles, contains sorted positions
      };

      typedef shared_ptr<Configuration> ConfigurationPtr;

      std::vector<ConfigurationPtr> configurations;

      int maxConfigs;
    };
  }
}

#endif
