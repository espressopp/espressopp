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

    void ConfigurationExt::set(size_t index, real x, real y, real z, real vx, real vy, real vz)
    {
      coordinates[index] = Real3D(x, y, z);
      velocities[index] = Real3D(vx, vy, vz);
    }

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

    ConfigurationExtIterator ConfigurationExt::getIterator()
    {
      return ConfigurationExtIterator(coordinates);
    }

    ConfigurationExtIterator::ConfigurationExtIterator(std::map<size_t, Real3D>& coordinates)
    {
      it = coordinates.begin();
      end = coordinates.end();
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

    Real3D ConfigurationExtIterator::nextCoordinates()
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

    void ConfigurationExt::registerPython() {
      using namespace espresso::python;

      class_<ConfigurationExtIterator>("ConfigurationExtIterator", no_init)
      .def("next", &ConfigurationExtIterator::nextId)
      .def("__iter__", pass_through)
      ;

      class_<ConfigurationExt>
        ("analysis_ConfigurationExt", no_init)
      .add_property("size", &ConfigurationExt::getSize)
      .def("__getitem__", &ConfigurationExt::getCoordinates)
      .def("__iter__", &ConfigurationExt::getIterator)
      ;
    }

  }
}
