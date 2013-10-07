// ESPP_CLASS
#ifndef _ANALYSIS_CONFIGURATIONEXT_HPP
#define _ANALYSIS_CONFIGURATIONEXT_HPP

#include "types.hpp"
#include "SystemAccess.hpp"
#include "RealND.hpp"
#include <map>

namespace espresso {
  namespace analysis {

    typedef shared_ptr<class ConfigurationExt> ConfigurationExtPtr;

    /** Iterator class for configuration to be used in Python */

    class ConfigurationExtIterator {

     public:

      ConfigurationExtIterator(std::map<size_t, RealND >& coordinates);
      //ConfigurationExtIterator(std::map<size_t, Real3D>& coordinates);

      /** Get next particle id for which properties are available */

      int nextId();
      const RealND nextProperties();
      //Real3D nextCoordinates();
      //Real3D nextVelocities();

     private:

      std::map<size_t, RealND >::iterator it;
      std::map<size_t, RealND >::iterator end;
      //std::map<size_t, Real3D>::iterator it;
      //std::map<size_t, Real3D>::iterator end;
    };

    /** Class that stores particle positions for later analysis. */

    class ConfigurationExt {

     public:

      ConfigurationExt(); //int nParticles

      ~ConfigurationExt();

      RealND getProperties(size_t id);
      //Real3D getCoordinates(size_t id);
      //Real3D getVelocities(size_t id);

      inline size_t getSize(){return particleProperties.size();}
      //size_t getSize();

      void set(size_t id, RealND vec) ;
      //void set(size_t id, real x, real y, real z, real vx, real vy, real vz);

      //int nParticles;     // number of particles of the configuration

      static void registerPython();

      class ConfigurationExtIterator getIterator();

     private:

      std::map<size_t, RealND > particleProperties;
//      std::map<size_t, Real3D> coordinates;
//      std::map<size_t, Real3D> velocities;
    };

  }
}

#endif
