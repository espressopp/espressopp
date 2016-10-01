/*
  Copyright (C) 2016
      Max Planck Institute for Polymer Research & JGU Mainz
  Copyright (C) 2012-2015
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI

  This file is part of ESPResSo++.

  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// ESPP_CLASS
#ifndef _ANALYSIS_CONFIGURATIONEXT_HPP
#define _ANALYSIS_CONFIGURATIONEXT_HPP

#include "SystemAccess.hpp"
#include "RealND.hpp"
#include <map>

namespace espressopp {
  namespace analysis {

    typedef shared_ptr<class ConfigurationExt> ConfigurationExtPtr;

    /** Iterator class for configuration to be used in Python */

    class ConfigurationExtIterator {

     public:

      ConfigurationExtIterator(std::map<size_t, RealND >& coordinates);
      //ConfigurationExtIterator(std::map<size_t, Real3D>& coordinates);

      size_t currentId();
      RealND currentProperties();

      void incrementIterator();

      /** Get next particle id for which properties are available */

      size_t nextId();
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
