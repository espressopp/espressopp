// ESPP_CLASS
#ifndef _ANALYSIS_CONFIGURATION_HPP
#define _ANALYSIS_CONFIGURATION_HPP

#include "types.hpp"
#include "SystemAccess.hpp"
#include <map>

namespace espresso {
  namespace analysis {

    typedef shared_ptr<class Configuration> ConfigurationPtr;

    /** Iterator class for configuration to be used in Python */

    class ConfigurationIterator {

     public:

      ConfigurationIterator(std::map<size_t, Real3D>& coordinates);

      /** Get next particle id for which coordinates are available */

      int nextId();
      Real3D nextCoordinates();

     private:

      std::map<size_t, Real3D>::iterator it;
      std::map<size_t, Real3D>::iterator end;
    };

    /** Class that stores particle positions for later analysis. */

    class Configuration {

     public:

      Configuration(); //int nParticles

      ~Configuration();

      Real3D getCoordinates(size_t id);

      size_t getSize();

      void set(size_t id, real x, real y, real z);

      //int nParticles;     // number of particles of the configuration

      static void registerPython();

      class ConfigurationIterator getIterator();

     private:

      std::map<size_t, Real3D> coordinates;
    };

  }
}

#endif
