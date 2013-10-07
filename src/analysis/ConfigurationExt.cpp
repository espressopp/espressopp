#include "python.hpp"
#include "ConfigurationExt.hpp"
#include <cmath>
#include <stdexcept>
#include "Real3D.hpp"

using namespace espresso;
using namespace boost::python;

#define DEFAULT_TAG 71

namespace espresso {
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

    const RealND& ConfigurationExt::getProperties(size_t index)
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

    int ConfigurationExtIterator::nextId()
    {
      if (it == end) {
        PyErr_SetString(PyExc_StopIteration, "No more data.");
        boost::python::throw_error_already_set();
      }

      int id = (*it).first;
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
      using namespace espresso::python;

      class_<ConfigurationExtIterator>("ConfigurationExtIterator", no_init)
      .def("next", &ConfigurationExtIterator::nextId)
      .def("__iter__", pass_through)
      ;

      class_<ConfigurationExt>
        ("analysis_ConfigurationExt", no_init)
      .add_property("size", &ConfigurationExt::getSize)
      //.def("__getitem__", &ConfigurationExt::getProperties) // TODO fix this!
      //.def("__getitem__", &ConfigurationExt::getCoordinates)
      .def("__iter__", &ConfigurationExt::getIterator)
      ;
    }

  }
}
