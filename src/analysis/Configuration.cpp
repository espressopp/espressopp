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

#include "python.hpp"
#include "Configuration.hpp"
#include <cmath>
#include <stdexcept>
#include "Real3D.hpp"

using namespace espressopp;
using namespace boost::python;

#define DEFAULT_TAG 71

namespace espressopp {
  namespace analysis {

    Configuration::Configuration(bool _pos, bool _vel, bool _force, bool _radius)
                : gatherPos(_pos), gatherVel(_vel), gatherForce(_force), gatherRadius(_radius)
    {
    }

    Configuration::Configuration() {
  	  gatherPos=true;
	  gatherVel=false;
	  gatherForce=false;
	  gatherRadius=false;
    }

    Configuration::~Configuration()
    {
    }

    void Configuration::set(size_t index, real x, real y, real z) {
      if (gatherPos)
        coordinates[index] = Real3D(x, y, z);
      else {
    	  std::cout << "Error: This configuration does not store coordinates" << std::endl;
      }
    }

    void Configuration::setCoordinates(size_t index, Real3D _pos) {
      if (gatherPos)
        coordinates[index] = _pos;
      else {
    	  std::cout << "Error: This configuration does not store coordinates" << std::endl;
      }
    }

    void Configuration::setVelocities(size_t index, Real3D _vel) {
      if (gatherVel)
        velocities[index] = _vel;
      else {
    	  std::cout << "Error: This configuration does not store velocities" << std::endl;
      }
    }

    void Configuration::setForces(size_t index, Real3D _forces) {
      if (gatherForce)
        forces[index] = _forces;
      else {
    	  std::cout << "Error: This configuration does not store forces" << std::endl;
      }
    }

    void Configuration::setRadius(size_t index, real _radius) {
      if (gatherRadius)
        radii[index] = _radius;
      else {
    	  std::cout << "Error: This configuration does not store radii" << std::endl;
      }
    }

    Real3D Configuration::getCoordinates(size_t index) {
      if (gatherPos)
        return coordinates[index];
      else {
    	  std::cout << "Error: This configuration has no information about coordinates" << std::endl;
    	  return Real3D(0,0,0);
      }
    }

    Real3D Configuration::getVelocities(size_t index) {
      if (gatherVel)
        return velocities[index];
      else {
    	  std::cout << "Error: This configuration has no information about velocities" << std::endl;
    	  return Real3D(0,0,0);
      }
    }

    Real3D Configuration::getForces(size_t index) {
      if (gatherForce)
        return forces[index];
      else {
    	  std::cout << "Error: This configuration has no information about forces" << std::endl;
    	  return Real3D(0,0,0);
      }
    }

    real Configuration::getRadius(size_t index) {
      if (gatherRadius)
        return radii[index];
      else {
    	  std::cout << "Error: This configuration has no information about radii" << std::endl;
    	  return 0;
      }
    }

    size_t Configuration::getSize() {
      size_t sizePos    = coordinates.size();
      size_t sizeVel    = velocities.size();
      size_t sizeForce  = forces.size();
      size_t sizeRadius = radii.size();

      if (sizePos    > 0) return sizePos;
      if (sizeVel    > 0) return sizeVel;
      if (sizeForce  > 0) return sizeForce;
      if (sizeRadius > 0) return sizeRadius;

      return 0;
    }

/*
    ConfigurationIterator Configuration::getIterator()
    {
      return ConfigurationIterator(coordinates);
    }

    ConfigurationIterator::ConfigurationIterator(std::map<size_t, Real3D>& coordinates)
    {
      it = coordinates.begin();
      end = coordinates.end();
    }

    int ConfigurationIterator::nextId()
    {
      if (it == end) {
        PyErr_SetString(PyExc_StopIteration, "No more data.");
        boost::python::throw_error_already_set();
      }

      int id = (*it).first;
      it++;
      return id;
    }

    Real3D ConfigurationIterator::nextCoordinates()
    {
      if (it == end) {
        PyErr_SetString(PyExc_StopIteration, "No more data.");
        boost::python::throw_error_already_set();
      }

      Real3D coords = (*it).first;
      it++;
      return coords;
    }

    inline object pass_through(object const& o) { return o; }
*/

    void Configuration::registerPython() {
      using namespace espressopp::python;
/*
      class_<ConfigurationIterator>("ConfigurationIterator", no_init)
      .def("next", &ConfigurationIterator::nextId)
      .def("__iter__", pass_through)
      ;
*/
      class_<Configuration>
        ("analysis_Configuration", no_init)
      .add_property("size", &Configuration::getSize)
      .def("__getitem__", &Configuration::getCoordinates)
      // .def("__iter__", &Configuration::getIterator)
      ;
    }

  }
}
