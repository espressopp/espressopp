/*
  Copyright (C) 2012,2013
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
#ifndef _ANALYSIS_CONFIGURATION_HPP
#define _ANALYSIS_CONFIGURATION_HPP

#include "types.hpp"
#include "SystemAccess.hpp"
#include <map>

namespace espressopp {
  namespace analysis {

    typedef shared_ptr<class Configuration> ConfigurationPtr;

    /** Iterator class for configuration to be used in Python
    class ConfigurationIterator {
    public:
      ConfigurationIterator(std::map<size_t, Real3D>& coordinates);
      // Get next particle id for which coordinates are available
      int nextId();
      Real3D nextCoordinates();
      Real3D nextVelocities();
      Real3D nextForces();
      Real3D nextRadius();
    private:
      std::map<size_t, Real3D>::iterator it;
      std::map<size_t, Real3D>::iterator end;
    };
    */

    /** Class that stores particle positions for later analysis. */
    class Configuration {
    public:
      Configuration();
      Configuration(bool _pos, bool _vel, bool _force, bool _radius);
      ~Configuration();
      Real3D getCoordinates(size_t id);
      Real3D getVelocities(size_t id);
      Real3D getForces(size_t id);
      real getRadius(size_t id);
      boost::python::list getIds();
      size_t getSize();
      void set(size_t id, real x, real y, real z);
      void setCoordinates(size_t id, Real3D _pos);
      void setVelocities(size_t id, Real3D _vel);
      void setForces(size_t id, Real3D _forces);
      void setRadius(size_t id, real _rad);
      static void registerPython();
      // class ConfigurationIterator getIterator();
    private:
      bool gatherPos, gatherVel, gatherForce, gatherRadius;
      std::map<size_t, Real3D> coordinates;
      std::map<size_t, Real3D> velocities;
      std::map<size_t, Real3D> forces;
      std::map<size_t, real> radii;
    };
  }
}

#endif
