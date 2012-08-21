#include "python.hpp"
#include "Configuration.hpp"
#include <cmath>
#include <stdexcept>
#include "Real3D.hpp"

using namespace espresso;
using namespace boost::python;

#define DEFAULT_TAG 71

namespace espresso {
  namespace analysis {

    Configuration::Configuration()  //int nParticles
    {
      //this->nParticles = nParticles;
    }

    Configuration::~Configuration()
    {
    }

    void Configuration::set(size_t index, real x, real y, real z)
    {
      coordinates[index] = Real3D(x, y, z);
    }

    Real3D Configuration::getCoordinates(size_t index)
    {
      return coordinates[index];
    }

    size_t Configuration::getSize()
    {
      return coordinates.size();
    }

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

    void Configuration::registerPython() {
      using namespace espresso::python;

      class_<ConfigurationIterator>("ConfigurationIterator", no_init)
      .def("next", &ConfigurationIterator::nextId)
      .def("__iter__", pass_through)
      ;

      class_<Configuration>
        ("analysis_Configuration", no_init)
      .add_property("size", &Configuration::getSize)
      .def("__getitem__", &Configuration::getCoordinates)
      .def("__iter__", &Configuration::getIterator)
      ;
    }

  }
}
