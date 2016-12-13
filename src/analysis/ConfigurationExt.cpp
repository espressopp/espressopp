/*
   Copyright (C) 2016
      Max Planck Institute for Polymer Research & JGU Mainz
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
#include "ConfigurationExt.hpp"
#include <cmath>
#include <stdexcept>
#include "Real3D.hpp"

using namespace espressopp;
using namespace boost::python;

#define DEFAULT_TAG 71

namespace espressopp {
  namespace analysis {

    ConfigurationExt::ConfigurationExt()  //int nParticles
    {
      //this->nParticles = nParticles;
    }

    ConfigurationExt::~ConfigurationExt()
    {
    }

    //FIXME
    //passing std::vector by value is very slow and should be avoided
    //however, passing reference is dangerous, since the original object might not survive until the end of this class
    //how should I use const here?
    void ConfigurationExt::set(size_t index, RealND vec)
    {
      particleProperties[index].setDimension(vec.getDimension());
      particleProperties[index] = vec;
      //coordinates[index] = Real3D(x, y, z);
      //velocities[index] = Real3D(vx, vy, vz);
    }

    RealND ConfigurationExt::getProperties(size_t index)
    {
      return particleProperties[index];
    }
    /*
    Real3D ConfigurationExt::getCoordinates(size_t index)
    {
      return coordinates[index];
    }

    Real3D ConfigurationExt::getVelocities(size_t index)
    {
      return velocities[index];
    }

    size_t ConfigurationExt::getSize()
    {
      return coordinates.size();
    }
    */

    ConfigurationExtIterator ConfigurationExt::getIterator()
    {
      return ConfigurationExtIterator(particleProperties);
      //return ConfigurationExtIterator(coordinates);
    }

    //ConfigurationExtIterator::ConfigurationExtIterator(std::map<size_t, Real3D>& coordinates)
    ConfigurationExtIterator::ConfigurationExtIterator(std::map<size_t, RealND >& particleProperties)
    {
      it = particleProperties.begin();
      end = particleProperties.end();
      //it = coordinates.begin();
      //end = coordinates.end();
    }

    size_t ConfigurationExtIterator::currentId()
    {
      if (it == end) {
        PyErr_SetString(PyExc_StopIteration, "No more data.");
        boost::python::throw_error_already_set();
      }

      size_t id = (*it).first;
      return id;
    }

    RealND ConfigurationExtIterator::currentProperties()
    {
      if (it == end) {
        PyErr_SetString(PyExc_StopIteration, "No more data.");
        boost::python::throw_error_already_set();
      }

      RealND props = (*it).second;
      return props;
    }

    void ConfigurationExtIterator::incrementIterator()
    {
      it++;
    }

    size_t ConfigurationExtIterator::nextId()
    {
      if (it == end) {
        PyErr_SetString(PyExc_StopIteration, "No more data.");
        boost::python::throw_error_already_set();
      }

      size_t id = (*it).first;
      it++;
      return id;
    }

    //bad performance - or the compiler is smart
    const RealND ConfigurationExtIterator::nextProperties()
    {
      if (it == end) {
        PyErr_SetString(PyExc_StopIteration, "No more data.");
        boost::python::throw_error_already_set();
      }

      RealND props = (*it).second;
      it++;
      return props;
    }
    /*
    Real3D ConfigurationExtIterator::nextCoordinates()
    {
      if (it == end) {
        PyErr_SetString(PyExc_StopIteration, "No more data.");
        boost::python::throw_error_already_set();
      }

      Real3D coords = (*it).second;
      it++;
      return coords;
    }
    */

    inline object pass_through(object const& o) { return o; }

    void ConfigurationExt::registerPython() {
      using namespace espressopp::python;

      class_<ConfigurationExtIterator>("ConfigurationExtIterator", no_init)
      .def("next", &ConfigurationExtIterator::nextId)
      .def("__iter__", pass_through)
      ;

      class_<ConfigurationExt>
        ("analysis_ConfigurationExt", no_init)
      .add_property("size", &ConfigurationExt::getSize)
      .def("__getitem__", &ConfigurationExt::getProperties)
      //.def("__getitem__", &ConfigurationExt::getCoordinates)
      .def("__iter__", &ConfigurationExt::getIterator)
      ;
    }

  }
}
