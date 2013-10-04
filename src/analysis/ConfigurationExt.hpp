// ESPP_CLASS
#ifndef _ANALYSIS_CONFIGURATIONEXT_HPP
#define _ANALYSIS_CONFIGURATIONEXT_HPP

#include "types.hpp"
#include "SystemAccess.hpp"
#include <map>

namespace espresso {
  namespace analysis {

    typedef shared_ptr<class ConfigurationExt> ConfigurationExtPtr;

    /** Iterator class for configuration to be used in Python */

    class ConfigurationExtIterator {

     public:

      ConfigurationExtIterator(std::map<size_t, std::vector<Real3D> >& coordinates);
      //ConfigurationExtIterator(std::map<size_t, Real3D>& coordinates);

      /** Get next particle id for which properties are available */

      int nextId();
      const std::vector<Real3D> nextProperties();
      //Real3D nextCoordinates();
      //Real3D nextVelocities();

     private:

      std::map<size_t, std::vector<Real3D> >::iterator it;
      std::map<size_t, std::vector<Real3D> >::iterator end;
      //std::map<size_t, Real3D>::iterator it;
      //std::map<size_t, Real3D>::iterator end;
    };

    /** Class that stores particle positions for later analysis. */

    class ConfigurationExt {

     public:

      ConfigurationExt(); //int nParticles

      ~ConfigurationExt();

      const std::vector<Real3D>& getProperties(size_t id);
      //Real3D getCoordinates(size_t id);
      //Real3D getVelocities(size_t id);

      inline size_t getSize(){return particleProperties.size();}
      //size_t getSize();

      void set(size_t id, std::vector<Real3D> vec) ;
      //void set(size_t id, real x, real y, real z, real vx, real vy, real vz);

      //int nParticles;     // number of particles of the configuration

      static void registerPython();

      class ConfigurationExtIterator getIterator();

     private:

      std::map<size_t, std::vector<Real3D> > particleProperties;
//      std::map<size_t, Real3D> coordinates;
//      std::map<size_t, Real3D> velocities;
    };

  }
}

#endif
